"""多壳层卫星集群（cluster）的组合分析。

在 orbit_solver 的 Mean-J2 解析模型之上，回答"把多条轨道组装成一个集群"时
必须先算清的三件事：

  1. 面共动（构型不散架）：各壳层的 J2 节面进动速率是否一致？
     不一致 → 轨道面相对张开，几天到几个月后整个拓扑就废了。
  2. 相位复位（拓扑确定性）：各壳层是否在同一个公共周期内跑整数圈？
     是 → 全网几何严格周期复现，路由表可以预编译。
  3. 链路可达（几何连通）：星间实际间距 vs. 激光视距上限。

设计流程的关键：**高度和倾角是两个独立旋钮**——
高度由共振阶数 k 定死，倾角再由进动锁定方程解出。
两个约束互不冲突，唯一的限制是进动匹配的高度天花板。

纯标准库实现，无第三方依赖。
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, field

import orbit_solver as osv

MU = osv.MU
R_E = osv.R_E
SECONDS_PER_DAY = osv.SECONDS_PER_DAY

# 判定"面已锁定"的差分进动阈值（度/天）。1e-6 度/天 ≈ 2700 年漂 1 度。
LOCK_TOLERANCE_DEG_DAY = 1e-6
# 判定"相位已复位"的单周期残余相位阈值（度）
REPEAT_TOLERANCE_DEG = 1.0


# ----------------------------------------------------------------------
# 壳层定义
# ----------------------------------------------------------------------
@dataclass
class Shell:
    """集群中的一个轨道壳层（一组同高度同倾角的卫星）。

    planes / sats_per_plane / phasing_f 为 Walker 星座记法 i: t/p/f，
    t = planes * sats_per_plane 为该壳层总星数。
    """

    name: str
    altitude_km: float
    inclination_deg: float
    planes: int = 1
    sats_per_plane: int = 1
    phasing_f: int = 0
    raan_span_deg: float = 360.0

    def __post_init__(self) -> None:
        osv._check_altitude(self.altitude_km)
        osv._check_inclination(self.inclination_deg)
        if self.planes < 1 or self.sats_per_plane < 1:
            raise ValueError(f"{self.name}: 轨道面数与每面星数必须 >= 1")
        if not (0 <= self.phasing_f < self.planes):
            raise ValueError(f"{self.name}: 相位因子 f 需满足 0 <= f < {self.planes}")
        if not (0.0 < self.raan_span_deg <= 360.0):
            raise ValueError(f"{self.name}: 升交点跨度需在 (0, 360] 度")

    @property
    def semi_major_axis_km(self) -> float:
        return R_E + self.altitude_km

    @property
    def num_satellites(self) -> int:
        return self.planes * self.sats_per_plane

    @property
    def walker_notation(self) -> str:
        return (f"{self.inclination_deg:.2f}°: "
                f"{self.num_satellites}/{self.planes}/{self.phasing_f}")

    def features(self) -> dict:
        return osv.get_orbit_features(self.altitude_km, self.inclination_deg)


# ----------------------------------------------------------------------
# 壳层内部几何
# ----------------------------------------------------------------------
def central_angle_between_planes(inc_deg: float, delta_raan_deg: float,
                                 arg_lat_deg: float) -> float:
    """同倾角、升交点相差 delta_raan 的两个面上，同纬度幅角处两星的地心张角。

        cos(theta) = cos(dRAAN) * (cos^2 u + sin^2 u * cos^2 i) + sin^2 u * sin^2 i

    u = 0（赤道）时 theta = dRAAN 取最大；u = 90°（最高纬度）时两面收拢，间距最小。
    """
    i = math.radians(inc_deg)
    d = math.radians(delta_raan_deg)
    u = math.radians(arg_lat_deg)
    cos_theta = (math.cos(d) * (math.cos(u) ** 2 + math.sin(u) ** 2 * math.cos(i) ** 2)
                 + math.sin(u) ** 2 * math.sin(i) ** 2)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_theta))))


def chord_km(a_km: float, central_angle_deg: float) -> float:
    """地心张角对应的两星直线距离（同高度圆轨道）。"""
    return 2.0 * a_km * math.sin(math.radians(central_angle_deg) / 2.0)


def intra_shell_geometry(shell: Shell,
                         margin_km: float = osv.DEFAULT_LOS_MARGIN_KM) -> dict:
    """单壳层内部的组网几何：同面间距、邻面间距、以及是否在激光视距内。"""
    a = shell.semi_major_axis_km
    in_plane_angle = 360.0 / shell.sats_per_plane
    in_plane_km = chord_km(a, in_plane_angle)

    raan_step = shell.raan_span_deg / shell.planes
    if shell.planes > 1:
        # 邻面同相位星：赤道处最远，最高纬度处最近
        equator_angle = central_angle_between_planes(
            shell.inclination_deg, raan_step, 0.0)
        apex_angle = central_angle_between_planes(
            shell.inclination_deg, raan_step, 90.0)
        cross_max_km = chord_km(a, equator_angle)
        cross_min_km = chord_km(a, apex_angle)
    else:
        raan_step = 0.0
        cross_max_km = cross_min_km = None

    los_limit = osv.solve_los_range(
        shell.altitude_km, shell.altitude_km, margin_km)["max_range_km"]
    reach = [in_plane_km] + ([cross_max_km] if cross_max_km is not None else [])
    return {
        "name": shell.name,
        "walker": shell.walker_notation,
        "num_satellites": shell.num_satellites,
        "in_plane_spacing_deg": in_plane_angle,
        "in_plane_spacing_km": in_plane_km,
        "raan_step_deg": raan_step,
        "cross_plane_max_km": cross_max_km,
        "cross_plane_min_km": cross_min_km,
        "los_limit_km": los_limit,
        "all_links_within_los": all(r <= los_limit for r in reach),
    }


# ----------------------------------------------------------------------
# 壳层两两关系
# ----------------------------------------------------------------------
def synodic_period(shell_a: Shell, shell_b: Shell) -> float:
    """两壳层的会合周期（相对相位跑完一整圈的时间），单位秒。

        T_syn = 1 / |1/T_a - 1/T_b|

    同高度返回 inf（相对相位冻结）。
    """
    t_a = osv.period_from_a(shell_a.semi_major_axis_km)
    t_b = osv.period_from_a(shell_b.semi_major_axis_km)
    delta = abs(1.0 / t_a - 1.0 / t_b)
    return float("inf") if delta == 0.0 else 1.0 / delta


def pair_analysis(shell_a: Shell, shell_b: Shell,
                  margin_km: float = osv.DEFAULT_LOS_MARGIN_KM,
                  drift_budget_deg: float = 1.0) -> dict:
    """两个壳层之间的共动性、会合节奏与链路几何。"""
    rate_a = osv.nodal_precession_rate(shell_a.semi_major_axis_km,
                                       shell_a.inclination_deg)
    rate_b = osv.nodal_precession_rate(shell_b.semi_major_axis_km,
                                       shell_b.inclination_deg)
    d_rate = osv.rad_s_to_deg_day(rate_a - rate_b)

    if abs(d_rate) <= LOCK_TOLERANCE_DEG_DAY:
        locked = True
        days_to_budget = float("inf")
    else:
        locked = False
        days_to_budget = abs(drift_budget_deg / d_rate)

    t_syn = synodic_period(shell_a, shell_b)
    los = osv.solve_los_range(shell_a.altitude_km, shell_b.altitude_km, margin_km)
    radial_gap = abs(shell_a.semi_major_axis_km - shell_b.semi_major_axis_km)

    return {
        "pair": f"{shell_a.name} <-> {shell_b.name}",
        "differential_precession_deg_day": d_rate,
        "planes_locked": locked,
        "days_to_drift_budget": days_to_budget,
        "drift_budget_deg": drift_budget_deg,
        "synodic_period_sec": t_syn,
        "synodic_period_hours": t_syn / 3600.0,
        "conjunctions_per_day": (0.0 if math.isinf(t_syn)
                                 else SECONDS_PER_DAY / t_syn),
        "min_range_km": radial_gap,
        "los_limit_km": los["max_range_km"],
        "inclination_gap_deg": abs(shell_a.inclination_deg - shell_b.inclination_deg),
    }


# ----------------------------------------------------------------------
# 集群级：相位复位
# ----------------------------------------------------------------------
def check_common_repeat(shells: list, T_common_hours: float,
                        tol_deg: float = REPEAT_TOLERANCE_DEG) -> dict:
    """检查各壳层是否在公共周期内跑整数圈，并给出单周期残余相位误差。

    残余相位 = 360° * (实际圈数 - 最近整数圈数)，即一个公共周期后
    卫星相对"理想复位点"的沿轨相位偏差。
    """
    if T_common_hours <= 0:
        raise ValueError("公共周期必须为正数")
    if not shells:
        raise ValueError("集群至少需要一个壳层")

    t_common = T_common_hours * 3600.0
    rows = []
    for shell in shells:
        period = osv.period_from_a(shell.semi_major_axis_km)
        revs = t_common / period
        k = round(revs)
        residual_deg = 360.0 * (revs - k)
        # 残余相位累积到 360° 所需的公共周期数
        cycles = (float("inf") if abs(residual_deg) < 1e-12
                  else abs(360.0 / residual_deg))
        rows.append({
            "name": shell.name,
            "period_min": period / 60.0,
            "revs_in_common_period": revs,
            "nearest_integer_revs": k,
            "residual_phase_deg": residual_deg,
            "resonant": abs(residual_deg) <= tol_deg and k >= 1,
            "cycles_to_full_lap": cycles,
        })
    return {
        "common_period_hours": T_common_hours,
        "tolerance_deg": tol_deg,
        "cluster_repeats": all(r["resonant"] for r in rows),
        "shells": rows,
    }


# ----------------------------------------------------------------------
# 集群级：总分析
# ----------------------------------------------------------------------
def analyze_cluster(shells: list,
                    common_period_hours: float = None,
                    margin_km: float = osv.DEFAULT_LOS_MARGIN_KM,
                    drift_budget_deg: float = 1.0) -> dict:
    """对一组壳层做完整的集群级体检。"""
    if not shells:
        raise ValueError("集群至少需要一个壳层")
    names = [s.name for s in shells]
    if len(set(names)) != len(names):
        raise ValueError("壳层名称必须唯一")

    per_shell = []
    for shell in shells:
        feat = shell.features()
        geom = intra_shell_geometry(shell, margin_km)
        per_shell.append({**feat, **geom})

    pairs = [pair_analysis(shells[i], shells[j], margin_km, drift_budget_deg)
             for i in range(len(shells)) for j in range(i + 1, len(shells))]

    rates = [osv.rad_s_to_deg_day(
        osv.nodal_precession_rate(s.semi_major_axis_km, s.inclination_deg))
        for s in shells]
    spread = max(rates) - min(rates) if len(rates) > 1 else 0.0
    fully_locked = spread <= LOCK_TOLERANCE_DEG_DAY

    repeat = (check_common_repeat(shells, common_period_hours)
              if common_period_hours is not None else None)

    unreachable = [p["pair"] for p in pairs if p["min_range_km"] > p["los_limit_km"]]

    return {
        "num_shells": len(shells),
        "num_satellites": sum(s.num_satellites for s in shells),
        "precession_rates_deg_day": dict(zip(names, rates)),
        "precession_spread_deg_day": spread,
        "planes_fully_locked": fully_locked,
        "worst_case_drift_days_per_deg": (
            float("inf") if fully_locked else 1.0 / spread),
        "shells": per_shell,
        "pairs": pairs,
        "common_repeat": repeat,
        "pairs_beyond_los": unreachable,
    }


# ----------------------------------------------------------------------
# 集群设计：两个独立旋钮
# ----------------------------------------------------------------------
def design_locked_cluster(base_shell: Shell, altitudes_km: list,
                          planes: int = 1, sats_per_plane: int = 1) -> dict:
    """旋钮二：给定基准壳，为每个目标高度解出锁定进动的倾角。"""
    shells = [base_shell]
    rejected = []
    for idx, h in enumerate(altitudes_km, start=1):
        match = osv.solve_matching_shell(
            base_shell.altitude_km, base_shell.inclination_deg, h)
        if not match["feasible"]:
            rejected.append({"altitude_km": h, "reason": match["reason"],
                             "max_h_km": match["max_h2_km"]})
            continue
        shells.append(Shell(
            name=f"S{idx}@{h:.0f}km",
            altitude_km=h,
            inclination_deg=match["matched_inc2_deg"],
            planes=planes,
            sats_per_plane=sats_per_plane,
        ))
    ceiling = osv.solve_matching_shell(
        base_shell.altitude_km, base_shell.inclination_deg,
        base_shell.altitude_km)["max_h2_km"]
    return {"shells": shells, "rejected": rejected, "ceiling_km": ceiling}


def design_resonant_locked_cluster(base_inclination_deg: float,
                                   common_period_hours: float,
                                   k_values: list,
                                   planes: int = 1,
                                   sats_per_plane: int = 1) -> dict:
    """两个旋钮一起拧：高度由共振阶数定死，倾角由进动锁定解出。

    以 k_values 中的第一个（圈数最多 = 最低）壳层为基准，其倾角取
    base_inclination_deg；其余壳层高度来自共振阶梯，倾角逐个求解。
    结果同时满足"面共动"与"相位复位"两个约束。
    """
    if not k_values:
        raise ValueError("至少需要一个共振阶数 k")
    if len(set(k_values)) != len(k_values):
        raise ValueError("共振阶数 k 不能重复")

    ordered = sorted(k_values, reverse=True)  # k 越大高度越低
    ladder = {k: osv.solve_resonant_shells(common_period_hours, k, k)[0]
              for k in ordered}
    below_surface = [k for k, row in ladder.items() if not row["above_surface"]]
    if below_surface:
        raise ValueError(f"共振阶数 {below_surface} 对应高度在地表以下")

    base_k = ordered[0]
    base_h = ladder[base_k]["altitude_km"]
    base = Shell(name=f"k={base_k}", altitude_km=base_h,
                 inclination_deg=base_inclination_deg,
                 planes=planes, sats_per_plane=sats_per_plane)

    shells = [base]
    rejected = []
    for k in ordered[1:]:
        h = ladder[k]["altitude_km"]
        match = osv.solve_matching_shell(base_h, base_inclination_deg, h)
        if not match["feasible"]:
            rejected.append({"k": k, "altitude_km": h,
                             "reason": match["reason"],
                             "max_h_km": match["max_h2_km"]})
            continue
        shells.append(Shell(name=f"k={k}", altitude_km=h,
                            inclination_deg=match["matched_inc2_deg"],
                            planes=planes, sats_per_plane=sats_per_plane))

    ceiling = osv.solve_matching_shell(base_h, base_inclination_deg,
                                       base_h)["max_h2_km"]
    return {
        "shells": shells,
        "rejected": rejected,
        "ceiling_km": ceiling,
        "common_period_hours": common_period_hours,
        "base_k": base_k,
    }


# ----------------------------------------------------------------------
# 报告
# ----------------------------------------------------------------------
def _fmt(value: float, spec: str = ".2f") -> str:
    if value is None:
        return "—"
    if isinstance(value, float) and math.isinf(value):
        return "∞"
    return format(value, spec)


def format_report(analysis: dict) -> str:
    lines = []
    add = lines.append

    add(f"集群规模: {analysis['num_shells']} 个壳层 / "
        f"{analysis['num_satellites']} 颗星")
    add("")
    add("── 壳层清单 ──")
    add(f"{'壳层':<14}{'高度km':>9}{'倾角°':>9}{'周期min':>10}"
        f"{'进动°/天':>11}{'同面间距km':>12}{'邻面间距km':>12}")
    for s in analysis["shells"]:
        cross = (_fmt(s["cross_plane_max_km"], ".0f")
                 if s["cross_plane_max_km"] is not None else "—")
        add(f"{s['name']:<14}{s['altitude_km']:>9.2f}{s['inclination_deg']:>9.3f}"
            f"{s['period_min']:>10.2f}{s['precession_deg_day']:>11.4f}"
            f"{s['in_plane_spacing_km']:>12.0f}{cross:>12}")

    add("")
    add("── 面共动（构型是否散架）──")
    add(f"进动速率极差: {analysis['precession_spread_deg_day']:.3e} 度/天")
    if analysis["planes_fully_locked"]:
        add("状态: 已锁定 ✓ 各面相对方位永久保持，无需组网维持机动")
    else:
        add(f"状态: 未锁定 ✗ 相对张开 1° 仅需 "
            f"{_fmt(analysis['worst_case_drift_days_per_deg'])} 天")

    if len(analysis["pairs"]) > 0:
        add("")
        add("── 壳层间关系 ──")
        add(f"{'配对':<26}{'差分进动°/天':>15}{'会合周期h':>12}"
            f"{'最近间距km':>12}{'视距上限km':>12}")
        for p in analysis["pairs"]:
            add(f"{p['pair']:<26}{p['differential_precession_deg_day']:>15.3e}"
                f"{_fmt(p['synodic_period_hours']):>12}"
                f"{p['min_range_km']:>12.1f}{p['los_limit_km']:>12.1f}")

    add("")
    add("── 链路可达 ──")
    blocked = [s for s in analysis["shells"] if not s["all_links_within_los"]]
    if blocked:
        for s in blocked:
            worst = max(x for x in (s["in_plane_spacing_km"],
                                    s["cross_plane_max_km"] or 0.0))
            add(f"✗ {s['name']}: 最远星间距 {worst:.0f} km "
                f"> 视距上限 {s['los_limit_km']:.0f} km（需加面/加星或中继）")
    else:
        add("✓ 各壳层内部的同面与邻面链路均在视距内")
    if analysis["pairs_beyond_los"]:
        add(f"✗ 超出视距的跨层配对: {', '.join(analysis['pairs_beyond_los'])}")

    repeat = analysis["common_repeat"]
    if repeat is not None:
        add("")
        add(f"── 相位复位（{repeat['common_period_hours']} 小时公共周期）──")
        add(f"{'壳层':<14}{'圈数':>12}{'最近整数':>10}{'残余相位°':>12}"
            f"{'复位':>8}")
        for r in repeat["shells"]:
            add(f"{r['name']:<14}{r['revs_in_common_period']:>12.6f}"
                f"{r['nearest_integer_revs']:>10}{r['residual_phase_deg']:>12.4f}"
                f"{('✓' if r['resonant'] else '✗'):>8}")
        add("全网拓扑严格周期复现: "
            + ("是 ✓" if repeat["cluster_repeats"] else "否 ✗"))

    return "\n".join(lines)


# ----------------------------------------------------------------------
# 演示与命令行
# ----------------------------------------------------------------------
def demo_naive_cluster() -> None:
    """反面教材：随手挑三个高度、都用同一个倾角。"""
    print("=" * 74)
    print("算例 A：随手拼的集群（三层都用 53°）—— 看它多快散架")
    print("=" * 74)
    shells = [
        Shell("低层@550", 550, 53, planes=24, sats_per_plane=20),
        Shell("中层@800", 800, 53, planes=18, sats_per_plane=16),
        Shell("高层@1200", 1200, 53, planes=12, sats_per_plane=12),
    ]
    print(format_report(analyze_cluster(shells, common_period_hours=24)))


def demo_locked_cluster() -> None:
    """旋钮二：锁定进动，但高度是随手挑的。"""
    print("\n" + "=" * 74)
    print("算例 B：只拧倾角旋钮 —— 面锁住了，但相位不复位")
    print("=" * 74)
    base = Shell("低层@550", 550, 53, planes=24, sats_per_plane=20)
    design = design_locked_cluster(base, [800, 1200], planes=18, sats_per_plane=16)
    for s in design["shells"][1:]:
        print(f"  {s.name}: 倾角需设为 {s.inclination_deg:.3f}°")
    print(f"  高度天花板: {design['ceiling_km']:.1f} km")
    print()
    print(format_report(analyze_cluster(design["shells"], common_period_hours=24)))


def demo_resonant_locked_cluster() -> None:
    """两个旋钮一起拧：高度取共振阶梯，倾角解进动锁定。"""
    print("\n" + "=" * 74)
    print("算例 C：两个旋钮一起拧 —— 面共动 + 相位复位同时成立")
    print("=" * 74)
    design = design_resonant_locked_cluster(
        base_inclination_deg=53.0,
        common_period_hours=24.0,
        k_values=[15, 14, 13],
        planes=24, sats_per_plane=20,
    )
    print(f"  基准壳: k={design['base_k']} 圈/天, 倾角 53°")
    for s in design["shells"]:
        print(f"  {s.name}: 高度 {s.altitude_km:.2f} km, "
              f"倾角 {s.inclination_deg:.3f}°")
    if design["rejected"]:
        for r in design["rejected"]:
            print(f"  ✗ k={r['k']} ({r['altitude_km']:.1f} km) 超天花板 "
                  f"{r['max_h_km']:.1f} km")
    print(f"  高度天花板: {design['ceiling_km']:.1f} km")
    print()
    print(format_report(analyze_cluster(design["shells"], common_period_hours=24)))


def run_demo() -> None:
    demo_naive_cluster()
    demo_locked_cluster()
    demo_resonant_locked_cluster()


def _parse_shell_spec(spec: str) -> Shell:
    """解析 'name:altitude:inclination[:planes:sats_per_plane]'。"""
    parts = spec.split(":")
    if len(parts) not in (3, 5):
        raise argparse.ArgumentTypeError(
            f"壳层格式应为 name:高度:倾角[:面数:每面星数]，收到 {spec!r}")
    try:
        name, h, inc = parts[0], float(parts[1]), float(parts[2])
        planes = int(parts[3]) if len(parts) == 5 else 1
        per_plane = int(parts[4]) if len(parts) == 5 else 1
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"壳层 {spec!r} 数值无法解析: {exc}") from exc
    try:
        return Shell(name, h, inc, planes=planes, sats_per_plane=per_plane)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="多壳层卫星集群的组合分析（Mean-J2 解析模型）")
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("analyze", help="体检一组给定的壳层")
    p.add_argument("shells", nargs="+", type=_parse_shell_spec,
                   metavar="name:高度:倾角[:面数:每面星数]")
    p.add_argument("--common-period", type=float, default=None,
                   help="用于检查相位复位的公共周期（小时）")

    p = sub.add_parser("lock", help="给定基准壳，为一组高度解出锁定倾角")
    p.add_argument("base_altitude_km", type=float)
    p.add_argument("base_inclination_deg", type=float)
    p.add_argument("altitudes_km", nargs="+", type=float)
    p.add_argument("--common-period", type=float, default=None)

    p = sub.add_parser("design", help="共振高度 + 锁定倾角，一次生成完整集群")
    p.add_argument("base_inclination_deg", type=float)
    p.add_argument("common_period_hours", type=float)
    p.add_argument("k_values", nargs="+", type=int)
    p.add_argument("--planes", type=int, default=1)
    p.add_argument("--sats-per-plane", type=int, default=1)

    sub.add_parser("demo", help="运行内置的三个对比算例（默认）")
    return parser


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)

    if args.command in (None, "demo"):
        run_demo()
        return 0

    if args.command == "analyze":
        print(format_report(analyze_cluster(
            args.shells, common_period_hours=args.common_period)))
        return 0

    if args.command == "lock":
        base = Shell("base", args.base_altitude_km, args.base_inclination_deg)
        design = design_locked_cluster(base, args.altitudes_km)
        for s in design["shells"][1:]:
            print(f"{s.name}: 倾角 {s.inclination_deg:.4f}°")
        for r in design["rejected"]:
            print(f"✗ {r['altitude_km']:.1f} km 超天花板 {r['max_h_km']:.1f} km")
        print(f"天花板: {design['ceiling_km']:.1f} km\n")
        print(format_report(analyze_cluster(
            design["shells"], common_period_hours=args.common_period)))
        return 0

    if args.command == "design":
        design = design_resonant_locked_cluster(
            args.base_inclination_deg, args.common_period_hours, args.k_values,
            planes=args.planes, sats_per_plane=args.sats_per_plane)
        for s in design["shells"]:
            print(f"{s.name}: 高度 {s.altitude_km:.2f} km, "
                  f"倾角 {s.inclination_deg:.4f}°")
        for r in design["rejected"]:
            print(f"✗ k={r['k']} ({r['altitude_km']:.1f} km) 超天花板 "
                  f"{r['max_h_km']:.1f} km")
        print()
        print(format_report(analyze_cluster(
            design["shells"], common_period_hours=args.common_period_hours)))
        return 0

    raise ValueError(f"未知命令 {args.command}")  # pragma: no cover


if __name__ == "__main__":
    raise SystemExit(main())
