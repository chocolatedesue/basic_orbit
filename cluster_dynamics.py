"""一个轨道周期内的集群相对运动与链路物理量。

前面的模块给的都是**标量指标**（周期、进动、容差）；这个模块给**时间序列**：
把 N 颗星放进 Hill/LVLH 旋转坐标系，沿一个轨道周期传播，逐时刻算出
星间距离、光速时延、光链路耦合效率。

两个容易踩的坑，本模块都做了处理：

1. **链路预算不能套远场平方反比。** 集群尺度只有几百米，而 2.5 cm 口径
   在 1550 nm 下的瑞利距离是 1267 m —— 光束根本还没发散完。必须用高斯
   光束模型；套 1/d² 会把链路余量的波动高估两个数量级。

2. **相对运动是否"呼吸"取决于构型选择，不是物理必然。** 有界 HCW 解里
   存在一族刚性构型（GCO），整个集群像刚体一样转，所有星间距离恒定。

坐标系：Hill/LVLH，x 径向（天顶为正）、y 沿轨、z 法向，单位米。
有界（无长期漂移）解：

    x(t) = rho_x * sin(n t + alpha)
    y(t) = 2 * rho_x * cos(n t + alpha) + y_c
    z(t) = rho_z * sin(n t + beta)

    GCO（刚性）: rho_z = sqrt(3) * rho_x, beta = alpha, y_c = 0
    PCO（呼吸）: rho_z = 2 * rho_x,       beta = alpha

纯标准库实现，无第三方依赖。
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import orbit_solver as osv

C_M_S = 299792458.0
DEFAULT_WAVELENGTH_M = 1.55e-6      # 1550 nm 通信波段
DEFAULT_APERTURE_RADIUS_M = 0.025   # 5 cm 口径


# ----------------------------------------------------------------------
# 构型定义
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class RelativeOrbit:
    """一颗星在 Hill 坐标系中的有界相对轨道（米 / 度）。"""

    name: str
    rho_x_m: float
    alpha_deg: float
    rho_z_m: float
    beta_deg: float
    y_offset_m: float = 0.0

    def __post_init__(self) -> None:
        if self.rho_x_m < 0 or self.rho_z_m < 0:
            raise ValueError(f"{self.name}: 振幅不能为负")

    def position(self, n_rad_s: float, t_sec: float) -> tuple:
        """t 时刻的 (x, y, z)，单位米。"""
        theta = n_rad_s * t_sec
        a = theta + math.radians(self.alpha_deg)
        b = theta + math.radians(self.beta_deg)
        return (self.rho_x_m * math.sin(a),
                2.0 * self.rho_x_m * math.cos(a) + self.y_offset_m,
                self.rho_z_m * math.sin(b))

    def range_from_centre_m(self, n_rad_s: float, t_sec: float) -> float:
        x, y, z = self.position(n_rad_s, t_sec)
        return math.sqrt(x * x + y * y + z * z)


def build_cluster(cluster_radius_m: float, ring_counts: list = None,
                  mode: str = "gco") -> list:
    """用同心环构造集群。

    mode="gco" —— 刚性构型：rho_z = sqrt(3) rho_x，离中心距离恒为 2*rho_x，
                  且**任意两星间距均不随时间变化**（见 demo 的数值验证）。
    mode="pco" —— 投影圆构型：rho_z = 2 rho_x，沿轨/法向投影是圆，
                  但径向分量摆动，星间距离随周期呼吸。
    """
    if cluster_radius_m <= 0:
        raise ValueError("集群半径必须为正数")
    if mode not in ("gco", "pco"):
        raise ValueError(f"未知构型 {mode!r}，可选 gco / pco")
    if ring_counts is None:
        ring_counts = [1, 8, 16, 24, 32]
    if not ring_counts or ring_counts[0] < 1:
        raise ValueError("环计数需非空且首环至少 1 颗星")

    shape = math.sqrt(3.0) if mode == "gco" else 2.0
    n_rings = len(ring_counts) - 1
    sats = []
    for ring_idx, count in enumerate(ring_counts):
        if count < 1:
            raise ValueError("每环至少 1 颗星")
        # 离中心的目标距离；GCO 下等于 2*rho_x
        radius = 0.0 if n_rings == 0 else cluster_radius_m * ring_idx / n_rings
        rho_x = radius / 2.0
        for j in range(count):
            alpha = 360.0 * j / count
            sats.append(RelativeOrbit(
                name=f"R{ring_idx}-{j:02d}",
                rho_x_m=rho_x, alpha_deg=alpha,
                rho_z_m=shape * rho_x, beta_deg=alpha))
    return sats


# ----------------------------------------------------------------------
# 光链路：高斯光束，不是远场平方反比
# ----------------------------------------------------------------------
def rayleigh_range_m(waist_m: float, wavelength_m: float) -> float:
    """瑞利距离 z_R = pi * w0^2 / lambda。远小于它 → 近场，光束尚未发散。"""
    if waist_m <= 0 or wavelength_m <= 0:
        raise ValueError("束腰与波长必须为正数")
    return math.pi * waist_m ** 2 / wavelength_m


def beam_radius_m(distance_m: float, waist_m: float,
                  wavelength_m: float) -> float:
    """传播距离 z 处的高斯光束半径 w(z) = w0 * sqrt(1 + (z/z_R)^2)。"""
    if distance_m < 0:
        raise ValueError("距离不能为负")
    z_r = rayleigh_range_m(waist_m, wavelength_m)
    return waist_m * math.sqrt(1.0 + (distance_m / z_r) ** 2)


def coupling_efficiency(distance_m: float,
                        waist_m: float = DEFAULT_APERTURE_RADIUS_M,
                        wavelength_m: float = DEFAULT_WAVELENGTH_M,
                        rx_aperture_radius_m: float = None) -> float:
    """接收口径截获的功率占比。

        eta = 1 - exp(-2 a^2 / w(z)^2)

    远场（z >> z_R）时退化为 1/d^2；近场时几乎不随距离变化。
    """
    if rx_aperture_radius_m is None:
        rx_aperture_radius_m = waist_m
    if rx_aperture_radius_m <= 0:
        raise ValueError("接收口径必须为正数")
    w = beam_radius_m(distance_m, waist_m, wavelength_m)
    return 1.0 - math.exp(-2.0 * rx_aperture_radius_m ** 2 / w ** 2)


def light_time_us(distance_m: float) -> float:
    """单程光速时延，微秒。"""
    if distance_m < 0:
        raise ValueError("距离不能为负")
    return distance_m / C_M_S * 1e6


# ----------------------------------------------------------------------
# 一个周期内的传播
# ----------------------------------------------------------------------
def _distances(positions: list) -> list:
    """给定时刻所有星的位置，返回按 (i, j) 顺序展平的两两距离。"""
    n = len(positions)
    out = []
    for i in range(n):
        xi, yi, zi = positions[i]
        for j in range(i + 1, n):
            xj, yj, zj = positions[j]
            out.append(math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2
                                 + (zi - zj) ** 2))
    return out


def _pair_stats(positions: list, dists: list) -> dict:
    """某时刻的集群尺度统计量（对全体求的，会掩盖单链路波动）。"""
    n = len(positions)
    nearest = [math.inf] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = dists[k]
            k += 1
            if d < nearest[i]:
                nearest[i] = d
            if d < nearest[j]:
                nearest[j] = d
    return {
        "min_pair_m": min(dists),
        "max_pair_m": max(dists),
        "nearest_min_m": min(nearest),
        "nearest_mean_m": sum(nearest) / n,
        "nearest_max_m": max(nearest),
    }


def propagate_period(sats: list, h_km: float, steps: int = 24,
                     waist_m: float = DEFAULT_APERTURE_RADIUS_M,
                     wavelength_m: float = DEFAULT_WAVELENGTH_M) -> dict:
    """沿一个轨道周期传播集群，逐时刻给出几何与链路物理量。"""
    if len(sats) < 2:
        raise ValueError("集群至少 2 颗星")
    if steps < 2:
        raise ValueError("采样点至少 2 个")
    osv._check_altitude(h_km)

    a = osv.R_E + h_km
    n = osv.mean_motion(a)
    period = osv.period_from_a(a)

    n_pairs = len(sats) * (len(sats) - 1) // 2
    pair_min = [math.inf] * n_pairs
    pair_max = [0.0] * n_pairs

    samples = []
    for k in range(steps):
        t = period * k / steps
        positions = [s.position(n, t) for s in sats]
        dists = _distances(positions)
        for idx, d in enumerate(dists):
            if d < pair_min[idx]:
                pair_min[idx] = d
            if d > pair_max[idx]:
                pair_max[idx] = d
        stats = _pair_stats(positions, dists)
        centre_ranges = [math.sqrt(x * x + y * y + z * z)
                         for x, y, z in positions]
        near = stats["nearest_mean_m"]
        samples.append({
            "t_sec": t,
            "phase_deg": 360.0 * k / steps,
            **stats,
            "centre_range_max_m": max(centre_ranges),
            "nearest_latency_us": light_time_us(near),
            "span_latency_us": light_time_us(stats["max_pair_m"]),
            "nearest_coupling": coupling_efficiency(
                near, waist_m, wavelength_m),
            "span_coupling": coupling_efficiency(
                stats["max_pair_m"], waist_m, wavelength_m),
        })

    def spread(key):
        values = [s[key] for s in samples]
        lo, hi = min(values), max(values)
        return {"min": lo, "max": hi, "ratio": (hi / lo) if lo > 0 else math.inf,
                "swing_pct": (hi - lo) / lo * 100.0 if lo > 0 else math.inf}

    # 逐链路波动 —— 集合统计量会把它平均掉，必须单独追踪
    swings = [(hi - lo) / lo * 100.0
              for lo, hi in zip(pair_min, pair_max) if lo > 0]
    worst_idx = max(range(len(swings)), key=swings.__getitem__)
    worst_lo, worst_hi = pair_min[worst_idx], pair_max[worst_idx]
    coupling_lo = coupling_efficiency(worst_hi, waist_m, wavelength_m)
    coupling_hi = coupling_efficiency(worst_lo, waist_m, wavelength_m)

    return {
        "per_link": {
            "worst_swing_pct": max(swings),
            "median_swing_pct": sorted(swings)[len(swings) // 2],
            "worst_min_m": worst_lo,
            "worst_max_m": worst_hi,
            "worst_coupling_swing_db": 10 * math.log10(coupling_hi / coupling_lo),
            "worst_latency_swing_us": light_time_us(worst_hi - worst_lo),
        },
        "num_satellites": len(sats),
        "altitude_km": h_km,
        "period_min": period / 60.0,
        "mean_motion_rad_s": n,
        "steps": steps,
        "waist_m": waist_m,
        "wavelength_m": wavelength_m,
        "rayleigh_range_m": rayleigh_range_m(waist_m, wavelength_m),
        "samples": samples,
        "nearest_spread": spread("nearest_mean_m"),
        "span_spread": spread("max_pair_m"),
        "coupling_spread": spread("nearest_coupling"),
        "is_rigid": max(swings) < 1e-9,
    }


def format_period_report(r: dict, max_rows: int = 12) -> str:
    lines = []
    add = lines.append
    add(f"集群: {r['num_satellites']} 星 @ {r['altitude_km']:.0f} km，"
        f"周期 {r['period_min']:.2f} min，采样 {r['steps']} 点")
    z_r = r["rayleigh_range_m"]
    near = r["nearest_spread"]["min"]
    span = r["span_spread"]["max"]
    add(f"光学: 束腰/口径 {r['waist_m'] * 100:.1f} cm，"
        f"波长 {r['wavelength_m'] * 1e9:.0f} nm，瑞利距离 {z_r:.0f} m")
    add(f"      最近邻 {near:.0f} m = {near / z_r:.2f} z_R"
        f"（{'近场' if near < z_r else '远场'}），"
        f"跨端 {span:.0f} m = {span / z_r:.2f} z_R"
        f"（{'近场' if span < z_r else '远场'}）")
    add("")
    add("── 一个周期内的时间序列 ──")
    add(f"{'相位°':>7}{'t(min)':>9}{'最近邻m':>10}{'跨端m':>10}"
        f"{'最近邻μs':>11}{'跨端μs':>10}{'耦合dB':>10}")
    stride = max(1, len(r["samples"]) // max_rows)
    for s in r["samples"][::stride]:
        add(f"{s['phase_deg']:>7.0f}{s['t_sec'] / 60:>9.2f}"
            f"{s['nearest_mean_m']:>10.1f}{s['max_pair_m']:>10.1f}"
            f"{s['nearest_latency_us']:>11.3f}{s['span_latency_us']:>10.3f}"
            f"{10 * math.log10(s['nearest_coupling']):>10.2f}")

    add("")
    add("── 一周期内的波动幅度 ──")
    for label, key in (("最近邻间距", "nearest_spread"),
                       ("集群跨端距离", "span_spread"),
                       ("最近邻耦合效率", "coupling_spread")):
        sp = r[key]
        add(f"{label:<16}{sp['min']:>12.4f} → {sp['max']:>10.4f}"
            f"   波动 {sp['swing_pct']:>8.4f}%")

    add("")
    add("── 逐链路波动（集合统计量会把它平均掉）──")
    pl = r["per_link"]
    add(f"最差单链路      : {pl['worst_min_m']:.2f} → {pl['worst_max_m']:.2f} m"
        f"   波动 {pl['worst_swing_pct']:.4f}%")
    add(f"中位单链路波动  : {pl['median_swing_pct']:.4f}%")
    add(f"最差耦合摆幅    : {pl['worst_coupling_swing_db']:.3f} dB")
    add(f"最差时延摆幅    : {pl['worst_latency_swing_us']:.4f} μs")

    add("")
    if r["is_rigid"]:
        add("✓ 刚性构型：所有星间距离在整个周期内恒定，链路预算无时变")
    else:
        add(f"✗ 呼吸构型：最差单链路距离波动 {r['per_link']['worst_swing_pct']:.1f}%，"
            f"耦合摆幅 {r['per_link']['worst_coupling_swing_db']:.2f} dB")
    return "\n".join(lines)


# ----------------------------------------------------------------------
# 演示
# ----------------------------------------------------------------------
def demo_rigidity() -> None:
    print("=" * 78)
    print("GCO（刚性） vs PCO（呼吸）：81 星 / 1 km 半径 / 650 km")
    print("=" * 78)
    for mode, label in (("gco", "GCO 刚性构型"), ("pco", "PCO 投影圆构型")):
        sats = build_cluster(1000.0, mode=mode)
        r = propagate_period(sats, 650.0, steps=48)
        print(f"\n【{label}】{len(sats)} 星")
        print(format_period_report(r, max_rows=6))


def demo_aperture_sweep() -> None:
    print("\n" + "=" * 78)
    print("口径扫描：近场 vs 远场，决定链路预算是否随周期变化")
    print("=" * 78)
    sats = build_cluster(1000.0, mode="pco")
    print(f"{'口径半径':>9}{'瑞利距离m':>11}{'单链路耦合摆幅':>16}"
          f"{'最近邻耦合':>12}{'跨端耦合':>11}")
    for w_cm in (0.3, 1.0, 2.5, 5.0):
        w = w_cm / 100.0
        r = propagate_period(sats, 650.0, steps=48, waist_m=w)
        near_db = 10 * math.log10(coupling_efficiency(196.0, w))
        span_db = 10 * math.log10(coupling_efficiency(2000.0, w))
        print(f"{w_cm:>7.1f}cm{r['rayleigh_range_m']:>11.0f}"
              f"{r['per_link']['worst_coupling_swing_db']:>14.3f}dB"
              f"{near_db:>11.1f}dB{span_db:>10.1f}dB")
    print("\n同一构型（距离摆幅恒为 11.8%），远场平方反比会预测 0.97 dB 摆幅；"
          "\n口径 ≥ 2.5 cm 时实测仅 0.008 dB —— 光束还没发散，平方反比不适用。")


def run_demo() -> None:
    demo_rigidity()
    demo_aperture_sweep()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="一个轨道周期内的集群相对运动与链路物理量")
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("period", help="传播一个周期并输出时间序列")
    p.add_argument("altitude_km", type=float)
    p.add_argument("cluster_radius_m", type=float)
    p.add_argument("--mode", choices=("gco", "pco"), default="gco")
    p.add_argument("--steps", type=int, default=24)
    p.add_argument("--aperture-cm", type=float, default=2.5)
    p.add_argument("--rings", type=int, nargs="+", default=None)

    p = sub.add_parser("link", help="单条链路的高斯光束预算")
    p.add_argument("distance_m", type=float)
    p.add_argument("--aperture-cm", type=float, default=2.5)

    sub.add_parser("demo", help="运行内置对比算例（默认）")
    return parser


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)

    if args.command in (None, "demo"):
        run_demo()
        return 0

    if args.command == "period":
        sats = build_cluster(args.cluster_radius_m, args.rings, args.mode)
        print(format_period_report(propagate_period(
            sats, args.altitude_km, args.steps,
            waist_m=args.aperture_cm / 100.0)))
        return 0

    if args.command == "link":
        w = args.aperture_cm / 100.0
        d = args.distance_m
        eta = coupling_efficiency(d, w)
        print(f"距离            : {d:.1f} m")
        print(f"束腰/口径半径   : {args.aperture_cm:.2f} cm")
        print(f"瑞利距离        : {rayleigh_range_m(w, DEFAULT_WAVELENGTH_M):.1f} m")
        print(f"该处光束半径    : {beam_radius_m(d, w, DEFAULT_WAVELENGTH_M) * 100:.2f} cm")
        print(f"耦合效率        : {eta:.4f} = {10 * math.log10(eta):.2f} dB")
        print(f"单程时延        : {light_time_us(d):.3f} μs")
        return 0

    raise ValueError(f"未知命令 {args.command}")  # pragma: no cover


if __name__ == "__main__":
    raise SystemExit(main())
