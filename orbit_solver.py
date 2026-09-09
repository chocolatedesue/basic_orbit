"""Mean-J2 解析轨道模型（Mean J2 Analytic Orbit Model）。

只保留两项物理：经典二体开普勒引力 + 一级 J2 地球扁率长期摄动。
不含大气阻力、太阳光压、三体引力、高阶重力场等非保守/高阶扰动，
因此全部结果都是闭式解析解，无需数值积分。

适用范围：巨型星座宏观拓扑设计、轨道面进动匹配、整数共振壳层、
近距离编队相对运动几何、星间激光视距。

纯标准库实现，无第三方依赖。
"""

from __future__ import annotations

import argparse
import json
import math

# ================= 地球物理常数 =================
MU = 398600.4418        # 地球引力常数, km^3/s^2
R_E = 6378.137          # 地球赤道半径, km
J2 = 0.00108263         # 二阶带谐项（扁率）系数
SECONDS_PER_DAY = 86400.0
DEFAULT_LOS_MARGIN_KM = 80.0   # 星间链路的大气擦地余量

# J2 进动速率 ∝ a^-3.5，跨壳层匹配的指数
PRECESSION_EXPONENT = 3.5


# ----------------------------------------------------------------------
# 基础量
# ----------------------------------------------------------------------
def mean_motion(a_km: float) -> float:
    """平均角速率 n = sqrt(mu / a^3)，单位 rad/s。"""
    _check_semi_major_axis(a_km)
    return math.sqrt(MU / (a_km ** 3))


def period_from_a(a_km: float) -> float:
    """开普勒周期 T = 2*pi*sqrt(a^3 / mu)，单位秒。"""
    return 2.0 * math.pi / mean_motion(a_km)


def a_from_period(period_sec: float) -> float:
    """周期反解半长轴 a = (mu * (T / 2pi)^2)^(1/3)，单位 km。"""
    if period_sec <= 0:
        raise ValueError("周期必须为正数")
    return (MU * (period_sec / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)


def nodal_precession_rate(a_km: float, inc_deg: float, ecc: float = 0.0) -> float:
    """升交点赤经进动速率 Omega_dot，单位 rad/s。

        Omega_dot = -3/2 * J2 * (R_E/a)^2 * n * cos(i) / (1 - e^2)^2

    顺行轨道（i < 90°）为负（西退），逆行轨道为正（东进）。
    """
    _check_eccentricity(ecc)
    n = mean_motion(a_km)
    p_factor = (1.0 - ecc ** 2) ** 2
    return (-1.5 * J2 * ((R_E / a_km) ** 2) * n
            * math.cos(math.radians(inc_deg)) / p_factor)


def perigee_precession_rate(a_km: float, inc_deg: float, ecc: float = 0.0) -> float:
    """近地点幅角进动速率 omega_dot，单位 rad/s。

        omega_dot = 3/4 * J2 * (R_E/a)^2 * n * (5 cos^2 i - 1) / (1 - e^2)^2

    i = 63.435° / 116.565° 时为零（冻结临界倾角）。
    """
    _check_eccentricity(ecc)
    n = mean_motion(a_km)
    cos_i = math.cos(math.radians(inc_deg))
    p_factor = (1.0 - ecc ** 2) ** 2
    return 0.75 * J2 * ((R_E / a_km) ** 2) * n * (5.0 * cos_i ** 2 - 1.0) / p_factor


def rad_s_to_deg_day(rate_rad_s: float) -> float:
    """rad/s 换算为 度/天。"""
    return math.degrees(rate_rad_s) * SECONDS_PER_DAY


# ----------------------------------------------------------------------
# 1. 单星 / 单壳层基础特征（“两个时钟”）
# ----------------------------------------------------------------------
def get_orbit_features(h_km: float, inc_deg: float, ecc: float = 0.0) -> dict:
    """求解单星/单壳层的基本时钟特征。

    时钟 B = 公转周期（跑一圈多久）；时钟 A = 轨道面进动（平面每天转多少度）。
    """
    _check_altitude(h_km)
    _check_inclination(inc_deg)
    a = R_E + h_km

    n = mean_motion(a)
    period_sec = period_from_a(a)
    prec_rad_s = nodal_precession_rate(a, inc_deg, ecc)
    peri_rad_s = perigee_precession_rate(a, inc_deg, ecc)

    return {
        "altitude_km": h_km,
        "inclination_deg": inc_deg,
        "eccentricity": ecc,
        "semi_major_axis_km": a,
        "mean_motion_rad_s": n,
        "period_sec": period_sec,
        "period_min": period_sec / 60.0,
        "revs_per_day": SECONDS_PER_DAY / period_sec,
        "precession_rad_s": prec_rad_s,
        "precession_deg_day": rad_s_to_deg_day(prec_rad_s),
        "perigee_rate_deg_day": rad_s_to_deg_day(peri_rad_s),
    }


# ----------------------------------------------------------------------
# 2. 跨层组网匹配（锁定两层轨道面永不分家）
# ----------------------------------------------------------------------
def solve_matching_shell(h1_km: float, inc1_deg: float, h2_km: float) -> dict:
    """给定基准壳 1，求壳 2 的倾角，使两层的 J2 进动速率完全相等。

    进动速率 ∝ a^-3.5 * cos(i)，令两者相等即得闭式解：

        cos(i2) = cos(i1) * (a2 / a1)^3.5

    因 |cos i2| <= 1，壳 2 存在高度天花板：

        a_max = a1 * (1 / |cos i1|)^(1/3.5)
    """
    _check_altitude(h1_km)
    _check_altitude(h2_km)
    _check_inclination(inc1_deg)

    a1 = R_E + h1_km
    a2 = R_E + h2_km
    cos_i1 = math.cos(math.radians(inc1_deg))

    # 高度天花板：极轨（cos i1 = 0）不存在上限
    if abs(cos_i1) < 1e-15:
        h2_max = float("inf")
    else:
        a2_max = a1 * ((1.0 / abs(cos_i1)) ** (1.0 / PRECESSION_EXPONENT))
        h2_max = a2_max - R_E

    cos_i2 = cos_i1 * ((a2 / a1) ** PRECESSION_EXPONENT)

    if abs(cos_i2) > 1.0:
        return {
            "feasible": False,
            "h2_km": h2_km,
            "max_h2_km": h2_max,
            "reason": f"cos(i2) = {cos_i2:.4f} 超界，壳2 高度需 <= {h2_max:.1f} km",
        }

    inc2_deg = math.degrees(math.acos(cos_i2))
    return {
        "feasible": True,
        "h2_km": h2_km,
        "matched_inc2_deg": inc2_deg,
        "max_h2_km": h2_max,
        "shared_precession_deg_day": rad_s_to_deg_day(
            nodal_precession_rate(a1, inc1_deg)
        ),
    }


def solve_sun_synchronous_inclination(h_km: float, ecc: float = 0.0) -> dict:
    """太阳同步倾角：令 Omega_dot 等于地球公转平均角速率（+360°/回归年）。"""
    _check_altitude(h_km)
    a = R_E + h_km
    target_rad_s = 2.0 * math.pi / (365.2421897 * SECONDS_PER_DAY)
    n = mean_motion(a)
    cos_i = -target_rad_s * ((1.0 - ecc ** 2) ** 2) / (1.5 * J2 * ((R_E / a) ** 2) * n)
    if abs(cos_i) > 1.0:
        return {"feasible": False, "altitude_km": h_km,
                "reason": f"cos(i) = {cos_i:.4f} 超界，该高度无太阳同步解"}
    return {
        "feasible": True,
        "altitude_km": h_km,
        "inclination_deg": math.degrees(math.acos(cos_i)),
        "target_precession_deg_day": rad_s_to_deg_day(target_rad_s),
    }


# ----------------------------------------------------------------------
# 3. 确定性周期网络（整数共振壳层）
# ----------------------------------------------------------------------
def solve_resonant_shells(T_common_hours: float, k_min: int, k_max: int) -> list:
    """给定全网公共复位周期，生成一组整数圈数谐振的壳层高度。

    壳层在 T_common 内跑整整 k 圈：T_shell = T_common / k，再反解高度。

    注意：这是惯性系下的周期共振（全网相位复位），不是重复星下点轨迹
    ——后者还需计入地球自转与 J2 对交点周期的修正。
    """
    if T_common_hours <= 0:
        raise ValueError("公共周期必须为正数")
    if k_min < 1 or k_max < k_min:
        raise ValueError("圈数需满足 1 <= k_min <= k_max")

    T_common_sec = T_common_hours * 3600.0
    shells = []
    for k in range(k_min, k_max + 1):
        T_shell = T_common_sec / k
        a = a_from_period(T_shell)
        h = a - R_E
        shells.append({
            "revs_per_common_period": k,
            "altitude_km": h,
            "semi_major_axis_km": a,
            "period_min": T_shell / 60.0,
            "above_surface": h > 0.0,
        })
    return shells


# ----------------------------------------------------------------------
# 4. 百米级超紧密编队（2:1 椭圆）
# ----------------------------------------------------------------------
def solve_formation_geometry(h_km: float, along_track_range_m: float) -> dict:
    """求解 2:1 相对运动椭圆的偏心率。

    仅靠微小偏心率差产生的相对运动（CW 方程的周期解）：

        dx(t) =      a*e*cos(n t)     径向，幅值 ±a*e
        dy(t) = -2 * a*e*sin(n t)     沿轨，幅值 ±2*a*e

    输入 along_track_range_m = 沿轨方向的单边幅值 L（米），则 2*a*e = L。
    """
    _check_altitude(h_km)
    if along_track_range_m <= 0:
        raise ValueError("编队尺度必须为正数")

    a = R_E + h_km
    L_km = along_track_range_m / 1000.0
    e = (L_km / 2.0) / a
    radial_m = a * e * 1000.0
    along_m = 2.0 * a * e * 1000.0
    return {
        "altitude_km": h_km,
        "semi_major_axis_km": a,
        "required_eccentricity": e,
        "radial_amplitude_m": radial_m,
        "along_track_amplitude_m": along_m,
        "period_min": period_from_a(a) / 60.0,
    }


def formation_offset(h_km: float, ecc: float, t_sec: float) -> dict:
    """2:1 椭圆在 t 时刻的相对坐标（米）。t=0 位于径向最大偏移处。"""
    _check_altitude(h_km)
    _check_eccentricity(ecc)
    a = R_E + h_km
    n = mean_motion(a)
    phase = n * t_sec
    return {
        "t_sec": t_sec,
        "radial_dx_m": a * ecc * math.cos(phase) * 1000.0,
        "along_track_dy_m": -2.0 * a * ecc * math.sin(phase) * 1000.0,
    }


# ----------------------------------------------------------------------
# 5. 星间激光视距（LOS）
# ----------------------------------------------------------------------
def solve_los_range(h1_km: float, h2_km: float,
                    margin_km: float = DEFAULT_LOS_MARGIN_KM) -> dict:
    """两星之间不被地球及大气遮挡的最大直线通信距离。

        R_tangent = R_E + margin
        D_max = sqrt((R_E+h1)^2 - R_t^2) + sqrt((R_E+h2)^2 - R_t^2)
    """
    _check_altitude(h1_km)
    _check_altitude(h2_km)
    if margin_km < 0:
        raise ValueError("大气余量不能为负")
    r_t = R_E + margin_km
    if h1_km <= margin_km or h2_km <= margin_km:
        raise ValueError(f"卫星高度需高于擦地余量 {margin_km} km")

    r1 = R_E + h1_km
    r2 = R_E + h2_km
    leg1 = math.sqrt(r1 ** 2 - r_t ** 2)
    leg2 = math.sqrt(r2 ** 2 - r_t ** 2)
    d_max = leg1 + leg2

    # 最大链路对应的地心张角（两星连线与切点的几何）
    half_angle_deg = (math.degrees(math.acos(r_t / r1))
                      + math.degrees(math.acos(r_t / r2)))
    return {
        "h1_km": h1_km,
        "h2_km": h2_km,
        "margin_km": margin_km,
        "tangent_radius_km": r_t,
        "max_range_km": d_max,
        "sat1_horizon_leg_km": leg1,
        "sat2_horizon_leg_km": leg2,
        "max_central_angle_deg": half_angle_deg,
    }


# ----------------------------------------------------------------------
# 输入校验
# ----------------------------------------------------------------------
def _check_semi_major_axis(a_km: float) -> None:
    if a_km <= 0:
        raise ValueError("半长轴必须为正数")


def _check_altitude(h_km: float) -> None:
    if h_km <= 0:
        raise ValueError(f"高度必须大于 0（地表以上），收到 {h_km} km")


def _check_inclination(inc_deg: float) -> None:
    if not (0.0 <= inc_deg <= 180.0):
        raise ValueError(f"倾角需在 [0, 180] 度之间，收到 {inc_deg}")


def _check_eccentricity(ecc: float) -> None:
    if not (0.0 <= ecc < 1.0):
        raise ValueError(f"偏心率需在 [0, 1) 之间，收到 {ecc}")


# ----------------------------------------------------------------------
# 演示与命令行
# ----------------------------------------------------------------------
def run_demo() -> None:
    print("--- 1. 星链基准壳特征 (550 km, 53°) ---")
    shell1 = get_orbit_features(550, 53)
    print(f"公转周期  : {shell1['period_min']:.2f} 分钟 "
          f"({shell1['revs_per_day']:.3f} 圈/天)")
    print(f"面进动速度: {shell1['precession_deg_day']:.3f} 度/天")

    print("\n--- 2. 配对计算: 给定壳1(550km, 53°)，在 1000km 配壳2 ---")
    match = solve_matching_shell(550, 53, 1000)
    print(f"配对是否可行  : {match['feasible']}")
    print(f"壳2需设置倾角 : {match['matched_inc2_deg']:.2f}°")
    print(f"共同进动速度  : {match['shared_precession_deg_day']:.3f} 度/天")
    print(f"该构型天花板  : {match['max_h2_km']:.1f} km（再高倾角就补不动了）")

    print("\n--- 3. 极轨配对天花板验证 (780 km, 86.4°) ---")
    polar = solve_matching_shell(780, 86.4, 2000)
    print(f"2000km 配对倾角: {polar['matched_inc2_deg']:.2f}°")
    print(f"近极轨天花板暴涨至: {polar['max_h2_km']:.1f} km")

    print("\n--- 4. 48 小时公共周期的一组谐振壳层 (取 30~33 圈) ---")
    for s in solve_resonant_shells(48, 30, 33):
        print(f"k={s['revs_per_common_period']:>3} 圈  "
              f"高度 {s['altitude_km']:>9.2f} km  "
              f"周期 {s['period_min']:.2f} 分钟")

    print("\n--- 5. 550km 高度做 1000 米算力编队所需的偏心率 ---")
    form = solve_formation_geometry(550, 1000)
    print(f"所需偏心率: e = {form['required_eccentricity']:.6g}")
    print(f"径向摆动  : ±{form['radial_amplitude_m']:.1f} m")
    print(f"沿轨摆动  : ±{form['along_track_amplitude_m']:.1f} m（严格 2 倍）")

    print("\n--- 6. 星间激光视距 (550 km <-> 1000 km, 80 km 擦地余量) ---")
    los = solve_los_range(550, 1000)
    print(f"最大直连距离: {los['max_range_km']:.1f} km")
    print(f"对应地心张角: {los['max_central_angle_deg']:.2f}°")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Mean-J2 解析轨道模型求解器（二体开普勒 + 一级 J2）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--json", action="store_true", help="以 JSON 输出结果")
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("features", help="单壳层周期与进动速率")
    p.add_argument("altitude_km", type=float)
    p.add_argument("inclination_deg", type=float)
    p.add_argument("--ecc", type=float, default=0.0)

    p = sub.add_parser("match", help="跨层进动匹配倾角")
    p.add_argument("h1_km", type=float)
    p.add_argument("inc1_deg", type=float)
    p.add_argument("h2_km", type=float)

    p = sub.add_parser("sunsync", help="太阳同步倾角")
    p.add_argument("altitude_km", type=float)
    p.add_argument("--ecc", type=float, default=0.0)

    p = sub.add_parser("resonant", help="整数共振壳层高度")
    p.add_argument("common_period_hours", type=float)
    p.add_argument("k_min", type=int)
    p.add_argument("k_max", type=int)

    p = sub.add_parser("formation", help="2:1 编队椭圆偏心率")
    p.add_argument("altitude_km", type=float)
    p.add_argument("along_track_m", type=float)

    p = sub.add_parser("los", help="星间激光最大视距")
    p.add_argument("h1_km", type=float)
    p.add_argument("h2_km", type=float)
    p.add_argument("--margin", type=float, default=DEFAULT_LOS_MARGIN_KM)

    sub.add_parser("demo", help="运行内置算例（默认）")
    return parser


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)

    if args.command in (None, "demo"):
        run_demo()
        return 0

    if args.command == "features":
        result = get_orbit_features(args.altitude_km, args.inclination_deg, args.ecc)
    elif args.command == "match":
        result = solve_matching_shell(args.h1_km, args.inc1_deg, args.h2_km)
    elif args.command == "sunsync":
        result = solve_sun_synchronous_inclination(args.altitude_km, args.ecc)
    elif args.command == "resonant":
        result = solve_resonant_shells(args.common_period_hours, args.k_min, args.k_max)
    elif args.command == "formation":
        result = solve_formation_geometry(args.altitude_km, args.along_track_m)
    elif args.command == "los":
        result = solve_los_range(args.h1_km, args.h2_km, args.margin)
    else:  # pragma: no cover - argparse 已拦截
        raise ValueError(f"未知命令 {args.command}")

    if args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        rows = result if isinstance(result, list) else [result]
        for row in rows:
            for key, value in row.items():
                if isinstance(value, float):
                    print(f"{key:<28}: {value:.6g}")
                else:
                    print(f"{key:<28}: {value}")
            if len(rows) > 1:
                print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
