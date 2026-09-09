"""重复星下点轨迹与地面接触周期。

orbit_solver 的 `solve_resonant_shells` 算的是**惯性系**相位复位；
但"什么时候能和地面通信"由**星下点轨迹**决定，那需要两个额外修正：

  1. 交点周期（draconitic period）而非开普勒周期——J2 让近地点和平近点角
     都长期漂移，两次过升交点的间隔比 2π/n 长
  2. 交点日（nodal day）而非恒星日——地球要追的是**进动中的**轨道面

650 km 太阳同步轨道上这两个修正合起来把"每天圈数"从 14.735 改成 14.717，
足以选错重复周期。

对轨道算力集群来说，重复星下点周期 M 天 = **地面接触模式完全复现的周期**，
也就是 I/O 可用性的基本节拍。

纯标准库实现，无第三方依赖。
"""

from __future__ import annotations

import argparse
import math

import orbit_solver as osv

MU = osv.MU
R_E = osv.R_E
J2 = osv.J2
OMEGA_EARTH = 7.2921159e-5      # 地球恒星自转角速率, rad/s
SECONDS_PER_DAY = osv.SECONDS_PER_DAY


# ----------------------------------------------------------------------
# 交点周期与交点日
# ----------------------------------------------------------------------
def perigee_rate(a_km: float, inc_deg: float, ecc: float = 0.0) -> float:
    """近地点幅角长期进动 omega_dot，rad/s。"""
    return osv.perigee_precession_rate(a_km, inc_deg, ecc)


def mean_anomaly_rate(a_km: float, inc_deg: float, ecc: float = 0.0) -> float:
    """平近点角长期变化率 M_dot（含 J2 一阶修正），rad/s。

        M_dot = n * [1 + 1.5 J2 (R/a)^2 sqrt(1-e^2) (1 - 1.5 sin^2 i) / (1-e^2)^2]
    """
    osv._check_eccentricity(ecc)
    n = osv.mean_motion(a_km)
    i = math.radians(inc_deg)
    factor = (1.5 * J2 * (R_E / a_km) ** 2 * math.sqrt(1.0 - ecc ** 2)
              * (1.0 - 1.5 * math.sin(i) ** 2) / ((1.0 - ecc ** 2) ** 2))
    return n * (1.0 + factor)


def nodal_period(a_km: float, inc_deg: float, ecc: float = 0.0) -> float:
    """交点周期（两次过升交点的间隔），秒。

        T_nodal = 2*pi / (M_dot + omega_dot)
    """
    rate = mean_anomaly_rate(a_km, inc_deg, ecc) + perigee_rate(a_km, inc_deg, ecc)
    if rate <= 0:
        raise ValueError("交点角速率非正，参数超出模型适用范围")
    return 2.0 * math.pi / rate


def nodal_day(a_km: float, inc_deg: float, ecc: float = 0.0) -> float:
    """交点日：地球相对**进动中的**轨道面转一圈的时间，秒。

        T_nodal_day = 2*pi / (omega_earth - Omega_dot)

    太阳同步轨道上 Omega_dot 恰好等于地球公转速率，此值精确等于一个平太阳日。
    """
    omega_dot = osv.nodal_precession_rate(a_km, inc_deg, ecc)
    rate = OMEGA_EARTH - omega_dot
    if rate <= 0:
        raise ValueError("交点日角速率非正，参数超出模型适用范围")
    return 2.0 * math.pi / rate


def revs_per_nodal_day(h_km: float, inc_deg: float, ecc: float = 0.0) -> float:
    """每交点日的圈数 —— 重复星下点判据的核心量。"""
    osv._check_altitude(h_km)
    osv._check_inclination(inc_deg)
    a = R_E + h_km
    return nodal_day(a, inc_deg, ecc) / nodal_period(a, inc_deg, ecc)


# ----------------------------------------------------------------------
# 重复星下点求解
# ----------------------------------------------------------------------
def _inclination_at(h_km: float, sun_synchronous: bool, inc_deg: float) -> float:
    if not sun_synchronous:
        return inc_deg
    sso = osv.solve_sun_synchronous_inclination(h_km)
    if not sso["feasible"]:
        raise ValueError(f"{h_km:.1f} km 无太阳同步解")
    return sso["inclination_deg"]


def solve_repeat_altitude(revs_n: int, days_m: int,
                          sun_synchronous: bool = True,
                          inc_deg: float = None,
                          h_lo: float = 200.0, h_hi: float = 2000.0,
                          tol_km: float = 1e-9) -> dict:
    """反解让星下点在 M 个交点日后精确重复 N 圈的高度。

    每交点日圈数随高度单调递减，用二分法求根。
    """
    if revs_n < 1 or days_m < 1:
        raise ValueError("N 与 M 必须为正整数")
    if not sun_synchronous and inc_deg is None:
        raise ValueError("非太阳同步时必须给定倾角")

    target = revs_n / days_m

    def excess(h):
        return revs_per_nodal_day(
            h, _inclination_at(h, sun_synchronous, inc_deg)) - target

    lo, hi = h_lo, h_hi
    f_lo, f_hi = excess(lo), excess(hi)
    if f_lo * f_hi > 0:
        raise ValueError(
            f"{revs_n}/{days_m} = {target:.5f} 圈/天 在 [{h_lo}, {h_hi}] km "
            f"内无解（端点圈数 {f_lo + target:.4f} / {f_hi + target:.4f}）")

    while hi - lo > tol_km:
        mid = 0.5 * (lo + hi)
        if excess(lo) * excess(mid) <= 0:
            hi = mid
        else:
            lo = mid
    h = 0.5 * (lo + hi)
    inc = _inclination_at(h, sun_synchronous, inc_deg)
    a = R_E + h
    return {
        "revs_n": revs_n,
        "days_m": days_m,
        "revs_per_nodal_day": target,
        "altitude_km": h,
        "inclination_deg": inc,
        "nodal_period_min": nodal_period(a, inc) / 60.0,
        "keplerian_period_min": osv.period_from_a(a) / 60.0,
        "nodal_day_sec": nodal_day(a, inc),
        "repeat_cycle_days": days_m,
        "ground_track_spacing_deg": 360.0 / revs_n * days_m,
    }


def find_repeat_options(h_km: float, inc_deg: float = None,
                        max_days: int = 20,
                        tol_revs: float = 0.02,
                        sun_synchronous: bool = True) -> list:
    """列出当前高度附近可用的 N/M 重复周期，及各自需要的高度微调。"""
    inc = _inclination_at(h_km, sun_synchronous, inc_deg)
    current = revs_per_nodal_day(h_km, inc)

    seen = set()
    options = []
    for m in range(1, max_days + 1):
        n = round(current * m)
        if n < 1:
            continue
        key = (n // math.gcd(n, m), m // math.gcd(n, m))
        if key in seen:
            continue
        seen.add(key)
        ratio = n / m
        if abs(ratio - current) > tol_revs:
            continue
        try:
            sol = solve_repeat_altitude(n, m, sun_synchronous, inc_deg)
        except ValueError:
            continue
        options.append({
            "revs_n": n,
            "days_m": m,
            "revs_per_nodal_day": ratio,
            "required_altitude_km": sol["altitude_km"],
            "altitude_shift_km": sol["altitude_km"] - h_km,
            "ground_track_spacing_deg": 360.0 / n * m,
        })
    options.sort(key=lambda o: abs(o["altitude_shift_km"]))
    return {"current_revs_per_nodal_day": current,
            "inclination_deg": inc, "options": options}


# ----------------------------------------------------------------------
# 地面接触窗口
# ----------------------------------------------------------------------
def contact_window(h_km: float, elevation_mask_deg: float = 10.0,
                   inc_deg: float = None,
                   sun_synchronous: bool = True) -> dict:
    """单次过顶的几何：地心半张角与最长可见时长。

        lambda_max = arccos(R_E * cos(eps) / (R_E + h)) - eps

    最长时长对应正上方过顶，实际斜过顶更短。
    """
    osv._check_altitude(h_km)
    if not (0.0 <= elevation_mask_deg < 90.0):
        raise ValueError("仰角门限需在 [0, 90) 度")
    inc = _inclination_at(h_km, sun_synchronous, inc_deg)
    a = R_E + h_km
    eps = math.radians(elevation_mask_deg)
    ratio = R_E * math.cos(eps) / a
    if ratio > 1.0:
        raise ValueError("几何无解：高度过低")
    lam = math.acos(ratio) - eps
    lam_deg = math.degrees(lam)
    t_nodal = nodal_period(a, inc)
    return {
        "altitude_km": h_km,
        "inclination_deg": inc,
        "elevation_mask_deg": elevation_mask_deg,
        "half_angle_deg": lam_deg,
        "max_pass_sec": (2.0 * lam_deg / 360.0) * t_nodal,
        "max_pass_min": (2.0 * lam_deg / 360.0) * t_nodal / 60.0,
        "footprint_radius_km": math.radians(lam_deg) * R_E,
        "nodal_period_min": t_nodal / 60.0,
    }


def passes_per_day_estimate(h_km: float, station_latitude_deg: float,
                            elevation_mask_deg: float = 10.0,
                            inc_deg: float = None,
                            sun_synchronous: bool = True) -> dict:
    """粗估某纬度地面站每天的过顶次数与总接触时长。

    一阶几何估计：把每天的 N 条星下点轨迹看成在该纬度上等经度间隔铺开，
    站点落在任一轨迹的 ±lambda_max 经度带内即算一次过顶，升/降轨各算一次。

        passes/day ≈ 2 * (2 * lambda_max / cos(phi)) / (360 / revs_per_day)

    **仅为量级估计**：忽略轨迹倾斜、忽略高纬度轨迹收敛、
    在 |phi| 接近或超过倾角时失效。要精确值必须做真实的星历传播。
    """
    inc = _inclination_at(h_km, sun_synchronous, inc_deg)
    phi = abs(station_latitude_deg)
    if phi >= 90.0:
        raise ValueError("站点纬度需在 (-90, 90) 度内")
    window = contact_window(h_km, elevation_mask_deg, inc_deg, sun_synchronous)
    reachable = min(inc, 180.0 - inc) + window["half_angle_deg"]
    revs = revs_per_nodal_day(h_km, inc)

    if phi > reachable:
        return {"station_latitude_deg": station_latitude_deg,
                "in_coverage": False, "passes_per_day": 0.0,
                "contact_min_per_day": 0.0, "duty_cycle": 0.0,
                "note": f"纬度超出可达范围 {reachable:.1f}°"}

    lon_band = 2.0 * window["half_angle_deg"] / max(math.cos(math.radians(phi)), 1e-6)
    passes = 2.0 * lon_band / (360.0 / revs)
    passes = min(passes, 2.0 * revs)
    # 平均过顶时长取最长值的 ~2/3（斜过顶折损）
    avg_pass_min = window["max_pass_min"] * 2.0 / 3.0
    contact_min = passes * avg_pass_min
    return {
        "station_latitude_deg": station_latitude_deg,
        "in_coverage": True,
        "revs_per_nodal_day": revs,
        "max_pass_min": window["max_pass_min"],
        "avg_pass_min": avg_pass_min,
        "passes_per_day": passes,
        "contact_min_per_day": contact_min,
        "duty_cycle": contact_min / (24 * 60),
        "note": "一阶几何估计，精确值需星历传播",
    }


# ----------------------------------------------------------------------
# 周期谱
# ----------------------------------------------------------------------
def periodicity_spectrum(h_km: float, inc_deg: float = None,
                         sun_synchronous: bool = True,
                         station_latitude_deg: float = 60.0,
                         elevation_mask_deg: float = 10.0,
                         repeat_days: int = None) -> dict:
    """列出该轨道上所有与算力调度相关的时间尺度。

    这是"周期性研究"的入口：把系统里每一个重复节拍摆在同一张表上，
    再看它们之间是否可公度（commensurable）。
    """
    inc = _inclination_at(h_km, sun_synchronous, inc_deg)
    a = R_E + h_km
    t_nodal = nodal_period(a, inc)
    window = contact_window(h_km, elevation_mask_deg, inc_deg, sun_synchronous)
    access = passes_per_day_estimate(h_km, station_latitude_deg,
                                     elevation_mask_deg, inc_deg, sun_synchronous)
    repeats = find_repeat_options(h_km, inc_deg, 20, 0.02, sun_synchronous)
    chosen = None
    if repeat_days is not None:
        chosen = next((o for o in repeats["options"]
                       if o["days_m"] == repeat_days), None)
    elif repeats["options"]:
        chosen = repeats["options"][0]

    scales = [
        {"name": "星间相对椭圆", "seconds": t_nodal,
         "governs": "链路距离 → 带宽（1/d²）"},
        {"name": "交点周期（基本节拍）", "seconds": t_nodal,
         "governs": "一切几何量的基频"},
        {"name": "单次过顶窗口", "seconds": window["max_pass_sec"],
         "governs": "单次 I/O 窗口长度"},
        {"name": "交点日", "seconds": nodal_day(a, inc),
         "governs": "星下点经度推进"},
    ]
    if chosen is not None:
        scales.append({
            "name": f"重复星下点（{chosen['revs_n']}圈/{chosen['days_m']}天）",
            "seconds": chosen["days_m"] * nodal_day(a, inc),
            "governs": "**地面接触模式完全复现**"})
    scales.append({"name": "太阳同步年周期", "seconds": 365.2422 * SECONDS_PER_DAY,
                   "governs": "beta 角 / 光照 / 热"})

    return {
        "altitude_km": h_km,
        "inclination_deg": inc,
        "is_sun_synchronous": sun_synchronous,
        "nodal_period_min": t_nodal / 60.0,
        "keplerian_period_min": osv.period_from_a(a) / 60.0,
        "j2_period_correction_sec": t_nodal - osv.period_from_a(a),
        "revs_per_nodal_day": revs_per_nodal_day(h_km, inc),
        "contact": window,
        "access": access,
        "repeat_options": repeats["options"],
        "chosen_repeat": chosen,
        "timescales": scales,
    }


def format_spectrum(s: dict) -> str:
    lines = []
    add = lines.append
    add(f"轨道: {s['altitude_km']:.2f} km / {s['inclination_deg']:.3f}°"
        + ("（太阳同步）" if s["is_sun_synchronous"] else ""))
    add("")
    add("── 周期修正 ──")
    add(f"开普勒周期    : {s['keplerian_period_min']:.3f} min")
    add(f"交点周期      : {s['nodal_period_min']:.3f} min "
        f"（J2 修正 {s['j2_period_correction_sec']:+.2f} s）")
    add(f"每交点日圈数  : {s['revs_per_nodal_day']:.5f}")

    add("")
    add("── 时间尺度谱 ──")
    add(f"{'节拍':<26}{'时长':>16}   支配")
    for t in s["timescales"]:
        sec = t["seconds"]
        if sec < 7200:
            shown = f"{sec / 60:.2f} min"
        elif sec < 3 * SECONDS_PER_DAY:
            shown = f"{sec / 3600:.2f} h"
        else:
            shown = f"{sec / SECONDS_PER_DAY:.2f} 天"
        add(f"{t['name']:<26}{shown:>16}   {t['governs']}")

    add("")
    add("── 地面接触 ──")
    c, ac = s["contact"], s["access"]
    add(f"仰角门限 {c['elevation_mask_deg']:.0f}° → 地心半张角 {c['half_angle_deg']:.2f}°, "
        f"覆盖半径 {c['footprint_radius_km']:.0f} km")
    add(f"最长过顶      : {c['max_pass_min']:.2f} min")
    if ac["in_coverage"]:
        add(f"纬度 {ac['station_latitude_deg']:.0f}° 站点: "
            f"≈{ac['passes_per_day']:.1f} 次/天, "
            f"≈{ac['contact_min_per_day']:.0f} min/天, "
            f"占空比 {ac['duty_cycle'] * 100:.1f}%")
        add(f"  （{ac['note']}）")
    else:
        add(f"纬度 {ac['station_latitude_deg']:.0f}° 站点: {ac['note']}")

    if s["repeat_options"]:
        add("")
        add("── 可选重复星下点周期 ──")
        add(f"{'N圈/M天':<14}{'圈/天':>12}{'需要高度km':>14}"
            f"{'高度微调km':>13}{'轨迹间隔°':>12}")
        for o in s["repeat_options"][:6]:
            add(f"{str(o['revs_n']) + '/' + str(o['days_m']):<14}"
                f"{o['revs_per_nodal_day']:>12.5f}"
                f"{o['required_altitude_km']:>14.3f}"
                f"{o['altitude_shift_km']:>+13.3f}"
                f"{o['ground_track_spacing_deg']:>12.3f}")
    return "\n".join(lines)


# ----------------------------------------------------------------------
# 演示与命令行
# ----------------------------------------------------------------------
def run_demo() -> None:
    print("=" * 76)
    print("650 km 晨昏太阳同步轨道的周期谱")
    print("=" * 76)
    print(format_spectrum(periodicity_spectrum(650, station_latitude_deg=60)))

    print("\n" + "=" * 76)
    print("锁定到 103圈/7天 重复星下点后的精确构型")
    print("=" * 76)
    sol = solve_repeat_altitude(103, 7)
    print(f"需要高度      : {sol['altitude_km']:.4f} km "
          f"（相对 650 km 调整 {sol['altitude_km'] - 650:+.3f} km）")
    print(f"太阳同步倾角  : {sol['inclination_deg']:.4f}°")
    print(f"交点周期      : {sol['nodal_period_min']:.4f} min")
    print(f"重复周期      : {sol['days_m']} 天 = {sol['revs_n']} 圈")
    print(f"轨迹经度间隔  : {sol['ground_track_spacing_deg']:.4f}°")
    check = revs_per_nodal_day(sol["altitude_km"], sol["inclination_deg"])
    print(f"校验 圈/交点日: {check:.9f}  (目标 {103 / 7:.9f})")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="重复星下点轨迹与地面接触周期（J2 交点周期修正）")
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("spectrum", help="列出该轨道的时间尺度谱")
    p.add_argument("altitude_km", type=float)
    p.add_argument("--inclination", type=float, default=None)
    p.add_argument("--station-latitude", type=float, default=60.0)
    p.add_argument("--elevation-mask", type=float, default=10.0)
    p.add_argument("--repeat-days", type=int, default=None)

    p = sub.add_parser("repeat", help="反解 N圈/M天 重复星下点的高度")
    p.add_argument("revs_n", type=int)
    p.add_argument("days_m", type=int)
    p.add_argument("--inclination", type=float, default=None)

    p = sub.add_parser("contact", help="单次过顶窗口几何")
    p.add_argument("altitude_km", type=float)
    p.add_argument("--elevation-mask", type=float, default=10.0)
    p.add_argument("--inclination", type=float, default=None)

    sub.add_parser("demo", help="运行内置算例（默认）")
    return parser


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)

    if args.command in (None, "demo"):
        run_demo()
        return 0

    sun_sync = getattr(args, "inclination", None) is None

    if args.command == "spectrum":
        print(format_spectrum(periodicity_spectrum(
            args.altitude_km, args.inclination, sun_sync,
            args.station_latitude, args.elevation_mask, args.repeat_days)))
        return 0

    if args.command == "repeat":
        sol = solve_repeat_altitude(args.revs_n, args.days_m, sun_sync,
                                    args.inclination)
        for key, value in sol.items():
            print(f"{key:<26}: "
                  + (f"{value:.6g}" if isinstance(value, float) else str(value)))
        return 0

    if args.command == "contact":
        for key, value in contact_window(
                args.altitude_km, args.elevation_mask, args.inclination,
                sun_sync).items():
            print(f"{key:<26}: "
                  + (f"{value:.6g}" if isinstance(value, float) else str(value)))
        return 0

    raise ValueError(f"未知命令 {args.command}")  # pragma: no cover


if __name__ == "__main__":
    raise SystemExit(main())
