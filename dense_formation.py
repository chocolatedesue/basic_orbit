"""百米~公里级密集编队的可行性分析。

cluster_solver 处理的是"多个壳层组成的宏观拓扑"；本模块处理另一个尺度：
**几十上百颗星挤在 1 km 半径内、彼此相距几百米**的紧密集群。

这个尺度上，成败不取决于绝对轨道，而取决于**差分量**——
两颗星的轨道根数差多少，决定了它们会不会飘散：

  * Δa（半长轴失配）→ 沿轨长期漂移，量级 3π·Δa 每圈。**这是头号杀手。**
  * Δi（倾角失配）  → 差分 J2 节面进动 → 横向长期смеar
  * ΔΩ（升交点差）  → 因为 Ω̇ 不依赖 Ω，**差分进动恒为零**

最后一条是密集集群能成立的关键：把横向铺开做成 ΔΩ 而不是 Δi，
一阶 J2 长期漂移天然为零（J2-invariant），只剩差分气动阻力需要管。

参考算例见 `demo_suncatcher()`：81 星 / 1 km 半径 / 650 km 晨昏太阳同步轨道。

纯标准库实现，无第三方依赖。
"""

from __future__ import annotations

import argparse
import math

import orbit_solver as osv

C_KM_S = 299792.458          # 真空光速, km/s
SECONDS_PER_DAY = osv.SECONDS_PER_DAY
DAYS_PER_YEAR = 365.25


# ----------------------------------------------------------------------
# 集群堆积几何
# ----------------------------------------------------------------------
def packing_spacing_m(n_sats: int, cluster_radius_m: float,
                      geometry: str = "hex") -> float:
    """n 颗星均匀铺在半径 R 的圆盘内时的最近邻间距。

    hex  —— 六方密排：每星占据面积 (sqrt(3)/2) d^2
    grid —— 方形网格：沿直径均分 sqrt(n) - 1 段
    """
    if n_sats < 2:
        raise ValueError("集群至少 2 颗星")
    if cluster_radius_m <= 0:
        raise ValueError("集群半径必须为正数")

    if geometry == "hex":
        area_per_sat = math.pi * cluster_radius_m ** 2 / n_sats
        return math.sqrt(area_per_sat / (math.sqrt(3) / 2.0))
    if geometry == "grid":
        side = math.sqrt(n_sats)
        if abs(side - round(side)) > 1e-9:
            raise ValueError(f"方形网格需要完全平方数，收到 {n_sats}")
        return 2.0 * cluster_radius_m / (round(side) - 1)
    raise ValueError(f"未知堆积方式 {geometry!r}，可选 hex / grid")


# ----------------------------------------------------------------------
# 链路：光速时延与自由空间功率
# ----------------------------------------------------------------------
def light_time_us(distance_m: float) -> float:
    """单程光速时延，单位微秒。"""
    if distance_m < 0:
        raise ValueError("距离不能为负")
    return distance_m / 1000.0 / C_KM_S * 1e6


def link_power_gain(far_m: float, near_m: float) -> float:
    """飞得更近带来的接收功率增益倍数（自由空间发散，功率 ∝ 1/d^2）。"""
    if far_m <= 0 or near_m <= 0:
        raise ValueError("距离必须为正数")
    return (far_m / near_m) ** 2


# ----------------------------------------------------------------------
# 头号杀手：半长轴失配导致的沿轨漂移
# ----------------------------------------------------------------------
def along_track_drift(h_km: float, delta_a_m: float) -> dict:
    """半长轴差 Δa 造成的沿轨相对漂移。

    周期差 ΔT/T = 1.5 Δa/a，一圈累积的沿轨错位为 3π·Δa。
    注意这是**长期项**：每圈都在累加，不会自己回来。
    """
    osv._check_altitude(h_km)
    a = osv.R_E + h_km
    period = osv.period_from_a(a)
    per_orbit_m = 3.0 * math.pi * delta_a_m
    orbits_per_day = SECONDS_PER_DAY / period
    return {
        "delta_a_m": delta_a_m,
        "period_min": period / 60.0,
        "orbits_per_day": orbits_per_day,
        "drift_per_orbit_m": per_orbit_m,
        "drift_per_day_m": per_orbit_m * orbits_per_day,
        "drift_per_year_km": per_orbit_m * orbits_per_day * DAYS_PER_YEAR / 1000.0,
    }


def semi_major_axis_tolerance_m(h_km: float, drift_budget_m: float,
                                days: float = 1.0) -> float:
    """给定漂移预算，反解允许的半长轴失配（米）。"""
    osv._check_altitude(h_km)
    if drift_budget_m <= 0 or days <= 0:
        raise ValueError("漂移预算与时长必须为正数")
    period = osv.period_from_a(osv.R_E + h_km)
    orbits = SECONDS_PER_DAY / period * days
    return drift_budget_m / (3.0 * math.pi * orbits)


# ----------------------------------------------------------------------
# 差分 J2：Δi 会漂，ΔΩ 不会
# ----------------------------------------------------------------------
def precession_sensitivity(h_km: float, inc_deg: float) -> dict:
    """节面进动速率对 a 和 i 的偏导。

        Omega_dot ∝ a^-3.5 * cos(i)
        d(Omega_dot)/da = -3.5 * Omega_dot / a
        d(Omega_dot)/di = -Omega_dot * tan(i)

    对 Omega 的偏导恒为零 —— 这正是密集集群横向铺开要用 ΔΩ 的原因。
    """
    osv._check_altitude(h_km)
    osv._check_inclination(inc_deg)
    a = osv.R_E + h_km
    rate = osv.rad_s_to_deg_day(osv.nodal_precession_rate(a, inc_deg))
    return {
        "precession_deg_day": rate,
        "d_rate_per_km_altitude": -3.5 * rate / a,
        "d_rate_per_deg_inclination": (
            -rate * math.tan(math.radians(inc_deg)) * math.pi / 180.0),
        "d_rate_per_deg_raan": 0.0,
    }


def differential_nodal_drift(h_km: float, inc_deg: float,
                             delta_a_m: float = 0.0,
                             delta_inc_deg: float = 0.0,
                             days: float = DAYS_PER_YEAR) -> dict:
    """两星根数差引起的差分节面进动，及其横向位移后果。"""
    a = osv.R_E + h_km
    sens = precession_sensitivity(h_km, inc_deg)
    d_rate = (sens["d_rate_per_km_altitude"] * delta_a_m / 1000.0
              + sens["d_rate_per_deg_inclination"] * delta_inc_deg)
    d_angle_deg = d_rate * days
    return {
        "delta_a_m": delta_a_m,
        "delta_inc_deg": delta_inc_deg,
        "differential_rate_deg_day": d_rate,
        "days": days,
        "raan_spread_deg": d_angle_deg,
        "cross_track_km": math.radians(abs(d_angle_deg)) * a,
    }


def raan_offset_for_cross_track(h_km: float, inc_deg: float,
                                cross_track_m: float) -> float:
    """用 ΔΩ 铺开横向间距所需的升交点差（度）。

        横向振幅 ≈ a * dRAAN * sin(i)

    这条路线的差分 J2 进动恒为零。
    """
    osv._check_altitude(h_km)
    osv._check_inclination(inc_deg)
    a_m = (osv.R_E + h_km) * 1000.0
    sin_i = math.sin(math.radians(inc_deg))
    if abs(sin_i) < 1e-12:
        raise ValueError("赤道轨道无法用 ΔΩ 产生横向间距")
    return math.degrees(cross_track_m / (a_m * sin_i))


def inclination_offset_for_cross_track(h_km: float, cross_track_m: float) -> float:
    """用 Δi 铺开横向间距所需的倾角差（度）。这条路线会长期漂移。"""
    osv._check_altitude(h_km)
    a_m = (osv.R_E + h_km) * 1000.0
    return math.degrees(cross_track_m / a_m)


# ----------------------------------------------------------------------
# 综合分析
# ----------------------------------------------------------------------
def analyze_dense_cluster(h_km: float, inc_deg: float, n_sats: int,
                          cluster_radius_m: float,
                          drift_budget_m: float = 100.0,
                          reference_isl_km: float = 1000.0) -> dict:
    """密集集群的完整可行性体检。"""
    features = osv.get_orbit_features(h_km, inc_deg)
    spacing = packing_spacing_m(n_sats, cluster_radius_m, "hex")

    # 相对运动 2:1 椭圆：沿轨单边幅值取集群半径
    ellipse = osv.solve_formation_geometry(h_km, cluster_radius_m)

    tol_day = semi_major_axis_tolerance_m(h_km, drift_budget_m, 1.0)
    tol_year = semi_major_axis_tolerance_m(h_km, drift_budget_m, DAYS_PER_YEAR)
    drift_table = [along_track_drift(h_km, da) for da in (1.0, 10.0, 100.0)]

    # 两条横向铺开路线的对比
    d_inc = inclination_offset_for_cross_track(h_km, cluster_radius_m)
    d_raan = raan_offset_for_cross_track(h_km, inc_deg, cluster_radius_m)
    via_inc = differential_nodal_drift(h_km, inc_deg, delta_inc_deg=d_inc)
    via_raan = differential_nodal_drift(h_km, inc_deg, delta_inc_deg=0.0)

    sso = osv.solve_sun_synchronous_inclination(h_km)

    return {
        "altitude_km": h_km,
        "inclination_deg": inc_deg,
        "num_satellites": n_sats,
        "cluster_radius_m": cluster_radius_m,
        "period_min": features["period_min"],
        "revs_per_day": features["revs_per_day"],
        "precession_deg_day": features["precession_deg_day"],
        "sun_synchronous_inclination_deg": (
            sso["inclination_deg"] if sso["feasible"] else None),
        "is_sun_synchronous": (
            sso["feasible"]
            and abs(sso["inclination_deg"] - inc_deg) < 0.05),
        "neighbor_spacing_m": spacing,
        "relative_ellipse": ellipse,
        "relative_ellipse_period_min": features["period_min"],
        "drift_budget_m": drift_budget_m,
        "delta_a_tolerance_day_m": tol_day,
        "delta_a_tolerance_year_m": tol_year,
        "drift_table": drift_table,
        "cross_track_via_inclination": {
            "delta_inc_deg": d_inc, **via_inc},
        "cross_track_via_raan": {
            "delta_raan_deg": d_raan, **via_raan},
        "neighbor_light_time_us": light_time_us(spacing),
        "cluster_light_time_us": light_time_us(2.0 * cluster_radius_m),
        "reference_isl_light_time_us": light_time_us(reference_isl_km * 1000.0),
        "link_power_gain": link_power_gain(reference_isl_km * 1000.0, spacing),
    }


def format_dense_report(r: dict) -> str:
    lines = []
    add = lines.append

    add(f"密集集群: {r['num_satellites']} 星 / {r['cluster_radius_m']:.0f} m 半径 "
        f"/ {r['altitude_km']:.0f} km / {r['inclination_deg']:.2f}°")
    add("")
    add("── 轨道基本盘 ──")
    add(f"周期          : {r['period_min']:.2f} min（{r['revs_per_day']:.3f} 圈/天）")
    add(f"节面进动      : {r['precession_deg_day']:+.4f} °/天")
    if r["sun_synchronous_inclination_deg"] is not None:
        add(f"太阳同步倾角  : {r['sun_synchronous_inclination_deg']:.3f}° "
            f"→ 当前{'满足 ✓' if r['is_sun_synchronous'] else '不满足 ✗'}")

    add("")
    add("── 集群几何 ──")
    add(f"最近邻间距    : {r['neighbor_spacing_m']:.0f} m（六方密排）")
    ell = r["relative_ellipse"]
    add(f"2:1 相对椭圆  : e = {ell['required_eccentricity']:.3e}, "
        f"径向 ±{ell['radial_amplitude_m']:.0f} m, "
        f"沿轨 ±{ell['along_track_amplitude_m']:.0f} m")
    add(f"椭圆周期      : {r['relative_ellipse_period_min']:.2f} min（= 轨道周期，闭合）")

    add("")
    add("── 头号杀手：Δa 沿轨漂移 ──")
    add(f"{'Δa':>8}{'每圈':>12}{'每天':>12}{'每年':>12}")
    for row in r["drift_table"]:
        add(f"{row['delta_a_m']:>6.0f} m{row['drift_per_orbit_m']:>10.1f} m"
            f"{row['drift_per_day_m'] / 1000:>10.2f} km"
            f"{row['drift_per_year_km']:>10.0f} km")
    add(f"漂移预算 {r['drift_budget_m']:.0f} m/天 → Δa 容差 "
        f"{r['delta_a_tolerance_day_m'] * 100:.1f} cm")
    add(f"漂移预算 {r['drift_budget_m']:.0f} m/年 → Δa 容差 "
        f"{r['delta_a_tolerance_year_m'] * 1000:.2f} mm")

    add("")
    add("── 横向铺开：两条路线 ──")
    inc_route = r["cross_track_via_inclination"]
    raan_route = r["cross_track_via_raan"]
    add(f"走 Δi : 需 {inc_route['delta_inc_deg']:.5f}° → "
        f"差分进动 {inc_route['differential_rate_deg_day']:.3e} °/天 → "
        f"一年横向散开 {inc_route['cross_track_km']:.1f} km ✗")
    add(f"走 ΔΩ : 需 {raan_route['delta_raan_deg']:.5f}° → "
        f"差分进动 {raan_route['differential_rate_deg_day']:.3e} °/天 → "
        f"一年横向散开 {raan_route['cross_track_km']:.1f} km ✓")
    add("→ Ω̇ 不依赖 Ω，同 a 同 i 只差 Ω 的集群一阶 J2 长期漂移为零")

    add("")
    add("── 链路 ──")
    add(f"最近邻时延    : {r['neighbor_light_time_us']:.3f} μs 单程")
    add(f"集群跨端时延  : {r['cluster_light_time_us']:.3f} μs 单程")
    add(f"对比长程 ISL  : {r['reference_isl_light_time_us']:.1f} μs")
    add(f"接收功率增益  : {r['link_power_gain']:.2e} 倍（1/d² 发散）")
    return "\n".join(lines)


# ----------------------------------------------------------------------
# 参考算例
# ----------------------------------------------------------------------
def demo_suncatcher() -> None:
    """复算公开报道的 Suncatcher 构型：81 星 / 1 km 半径 / 650 km 晨昏 SSO。"""
    print("=" * 70)
    print("参考算例：81 星 / 1 km 半径 / 650 km 晨昏太阳同步轨道")
    print("=" * 70)
    inc = osv.solve_sun_synchronous_inclination(650)["inclination_deg"]
    r = analyze_dense_cluster(650, inc, 81, 1000.0)
    print(format_dense_report(r))
    print()
    print("对照公开报道的次近邻间距 100~200 m：本模型六方密排给出 "
          f"{r['neighbor_spacing_m']:.0f} m")


def run_demo() -> None:
    demo_suncatcher()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="百米~公里级密集编队的可行性分析")
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("analyze", help="体检一个密集集群")
    p.add_argument("altitude_km", type=float)
    p.add_argument("n_sats", type=int)
    p.add_argument("cluster_radius_m", type=float)
    p.add_argument("--inclination", type=float, default=None,
                   help="默认取该高度的太阳同步倾角")
    p.add_argument("--drift-budget", type=float, default=100.0,
                   help="沿轨漂移预算（米/天），默认 100")

    p = sub.add_parser("tolerance", help="反解半长轴失配容差")
    p.add_argument("altitude_km", type=float)
    p.add_argument("drift_budget_m", type=float)
    p.add_argument("--days", type=float, default=1.0)

    sub.add_parser("demo", help="运行参考算例（默认）")
    return parser


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)

    if args.command in (None, "demo"):
        run_demo()
        return 0

    if args.command == "analyze":
        inc = args.inclination
        if inc is None:
            sso = osv.solve_sun_synchronous_inclination(args.altitude_km)
            if not sso["feasible"]:
                raise ValueError(f"{args.altitude_km} km 无太阳同步解，请显式给 --inclination")
            inc = sso["inclination_deg"]
        print(format_dense_report(analyze_dense_cluster(
            args.altitude_km, inc, args.n_sats, args.cluster_radius_m,
            drift_budget_m=args.drift_budget)))
        return 0

    if args.command == "tolerance":
        tol = semi_major_axis_tolerance_m(
            args.altitude_km, args.drift_budget_m, args.days)
        print(f"漂移 {args.drift_budget_m:.0f} m / {args.days:g} 天 "
              f"→ Δa 容差 {tol:.4g} m ({tol * 100:.2f} cm)")
        return 0

    raise ValueError(f"未知命令 {args.command}")  # pragma: no cover


if __name__ == "__main__":
    raise SystemExit(main())
