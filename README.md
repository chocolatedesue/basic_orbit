# basic_orbit — Mean-J2 解析轨道模型

只保留两项物理：**经典二体开普勒引力 + 一级 $J_2$ 地球扁率长期摄动**。
不含大气阻力、太阳光压、三体引力、高阶重力场，因此全部结果都是**闭式解析解**，
无需数值积分、无需仿真软件，纯标准库几毫秒出结果。

对巨型星座的宏观拓扑、轨道面进动匹配、共振分层、编队几何、星间视距，精度完全够用。

## 快速开始

```bash
python3 orbit_solver.py                       # 运行内置算例
python3 orbit_solver.py features 550 53       # 单壳层：周期 + 进动
python3 orbit_solver.py match 550 53 1000     # 跨层进动匹配倾角
python3 orbit_solver.py sunsync 800           # 太阳同步倾角
python3 orbit_solver.py resonant 48 30 33     # 48h 公共周期的共振壳层
python3 orbit_solver.py formation 550 1000    # 1000 m 编队所需偏心率
python3 orbit_solver.py los 550 1000          # 星间激光最大视距
python3 orbit_solver.py --json match 550 53 1000   # JSON 输出

python3 -m unittest discover -p 'test_*.py'    # 150 项解析校验
```

集群级（多轨道组合成一个 cluster）：

```bash
python3 cluster_solver.py                      # 三个对比算例
python3 cluster_solver.py design 53 24 15 14 13    # 共振高度 + 锁定倾角，一次生成
python3 cluster_solver.py lock 550 53 800 1200     # 给定基准壳，解各层锁定倾角
python3 cluster_solver.py analyze 低层:550:53:24:20 中层:800:47:18:16 --common-period 24
```

密集编队（百米~公里级紧密集群）：

```bash
python3 dense_formation.py                     # 81 星 / 1 km 半径参考算例
python3 dense_formation.py analyze 650 81 1000 # 高度 / 星数 / 集群半径
python3 dense_formation.py tolerance 650 100   # 反解 Δa 容差
```

重复星下点与地面接触周期：

```bash
python3 ground_track.py                        # 650 km SSO 周期谱
python3 ground_track.py spectrum 650 --station-latitude 78
python3 ground_track.py repeat 103 7           # 反解重复星下点高度
python3 ground_track.py contact 650            # 单次过顶窗口
```

作为库使用：

```python
from orbit_solver import get_orbit_features, solve_matching_shell

get_orbit_features(550, 53)["precession_deg_day"]   # -4.489
solve_matching_shell(550, 53, 1000)["matched_inc2_deg"]   # 41.40
```

## 常数

| 量 | 值 |
| --- | --- |
| 地球引力常数 $\mu$ | 398600.4418 km³/s² |
| 地球赤道半径 $R_E$ | 6378.137 km |
| $J_2$ 扁率系数 | 0.00108263 |
| 半长轴 | $a = R_E + h$ |

## 公式速查

### 1. 单星「两个时钟」

时钟 B —— 轨道周期：

$$T = 2\pi \sqrt{a^3/\mu}$$

时钟 A —— 节面进动角速度（顺行为负 / 西退，逆行为正）：

$$\dot{\Omega} = -\frac{3}{2} J_2 \left(\frac{R_E}{a}\right)^2 \sqrt{\frac{\mu}{a^3}} \frac{\cos i}{(1-e^2)^2}$$

附带的近地点幅角进动（$i = 63.435°$ 时为零，即冻结临界倾角）：

$$\dot{\omega} = \frac{3}{4} J_2 \left(\frac{R_E}{a}\right)^2 n \frac{5\cos^2 i - 1}{(1-e^2)^2}$$

→ `get_orbit_features(h_km, inc_deg, ecc=0.0)`

### 2. 跨层组网匹配

强制两层进动速率相等（$\dot{\Omega} \propto a^{-3.5}\cos i$）：

$$\cos i_2 = \cos i_1 \cdot (a_2/a_1)^{3.5}$$

由 $|\cos i_2| \le 1$ 得高度天花板：

$$a_{\max} = a_1 \cdot (1/|\cos i_1|)^{1/3.5}$$

→ `solve_matching_shell(h1_km, inc1_deg, h2_km)`；
另附 `solve_sun_synchronous_inclination(h_km)`（令 $\dot{\Omega} = +360°/回归年$）。

### 3. 整数共振壳层

给定全网公共复位周期 $T_{\text{common}}$，壳层跑整整 $k$ 圈：

$$T_{\text{shell}} = T_{\text{common}}/k, \qquad a_k = \left(\mu (T_{\text{shell}}/2\pi)^2\right)^{1/3}$$

→ `solve_resonant_shells(T_common_hours, k_min, k_max)`

### 4. 2:1 编队椭圆

同高度上仅靠微小偏心率差产生的周期性相对运动：

$$e = \frac{L \times 10^{-3}}{2a}, \qquad \Delta x(t) = ae\cos(nt), \qquad \Delta y(t) = -2ae\sin(nt)$$

$L$ 为沿轨单边幅值（米），径向幅值恒为其一半。

→ `solve_formation_geometry(h_km, along_track_range_m)`、`formation_offset(h_km, ecc, t_sec)`

### 5. 星间激光视距

擦地余量默认 80 km：

$$R_t = R_E + 80, \qquad D_{\max} = \sqrt{(R_E+h_1)^2 - R_t^2} + \sqrt{(R_E+h_2)^2 - R_t^2}$$

→ `solve_los_range(h1_km, h2_km, margin_km=80)`

## 典型结果

| 场景 | 结果 |
| --- | --- |
| 550 km / 53° | 周期 95.65 min，15.055 圈/天，进动 −4.489 °/天 |
| 550 km / 53° ↔ 1000 km 配对 | 高层需 **41.40°**，天花板 1631.8 km |
| 780 km / 86.4° ↔ 2000 km 配对 | 高层需 83.75°，天花板 9407.4 km |
| 48 h 复位 / 30 圈 | 566.90 km，周期 96.00 min |
| 550 km 上 1000 m 编队 | $e = 7.217\times10^{-5}$，径向 ±500 m，沿轨 ±1000 m |
| 550 km ↔ 1000 km 激光 | 最大 6076.1 km，地心张角 50.14° |

> 注：550 km/53° 与 1000 km 的配对倾角闭式解为 **41.40°**。若别处看到 44.13°，
> 那是笔误——把 $a_2/a_1 = 1.06495$ 代入 $\cos i_2 = \cos 53° \times 1.06495^{3.5} = 0.75005$，
> 即 $i_2 = 41.40°$，本仓库的 `test_known_pairing_value` 锁定了这个值。

## 集群组合分析（`cluster_solver.py`）

把多条轨道装配成**一个协同 cluster** 时，必须先算清三件事：

| 问题 | 判据 | 算不清的后果 |
| --- | --- | --- |
| **面共动**：构型会不会散架 | 各壳层 $\dot\Omega$ 是否一致 | 轨道面相对张开，几天后拓扑作废 |
| **相位复位**：拓扑是否确定 | 各壳层在公共周期内是否跑整数圈 | 相位漂移，路由表无法预编译 |
| **链路可达**：能不能连上 | 星间实际间距 vs. 激光视距上限 | 设计图上的链路物理上不存在 |

### 核心结论：高度和倾角是两个独立旋钮

这是把多壳层组成确定性 cluster 的关键——两个约束**互不冲突**，可以同时满足：

1. **高度旋钮**：由共振阶数 $k$ 定死，$a_k = (\mu(T_{\text{common}}/2\pi k)^2)^{1/3}$ → 保证相位复位
2. **倾角旋钮**：再由进动锁定方程解出，$\cos i_k = \cos i_1 (a_k/a_1)^{3.5}$ → 保证面共动

唯一的限制是进动匹配的高度天花板 $a_{\max} = a_1(1/|\cos i_1|)^{1/3.5}$。

`design_resonant_locked_cluster(53, 24, [15, 14, 13])` 一次拧完两个旋钮：

| 壳层 | 高度 km | 倾角 | 周期 min | 进动 °/天 | 24h 圈数 |
| --- | --- | --- | --- | --- | --- |
| k=15 | 566.90 | 53.000° | 96.00 | −4.4511 | 15.000000 |
| k=14 | 893.80 | 45.014° | 102.86 | −4.4511 | 14.000000 |
| k=13 | 1262.09 | 32.820° | 110.77 | −4.4511 | 13.000000 |

进动极差 $1.8\times10^{-15}$ °/天（浮点噪声量级）→ 轨道面永久共动；
残余相位严格为 0 → 全网几何每 24 小时精确复现。

副产物：相邻壳层每天正好会合一次（$T_{\text{syn}} = 24$ h），k=15 与 k=13 每 12 小时会合一次
——跨层接触节奏本身也是确定的。

### 三个对比算例

`python3 cluster_solver.py` 依次跑：

- **算例 A（反面教材）**：三层都用 53°。进动极差 1.209 °/天，**相对张开 1° 只要 0.83 天**，
  半年后轨道面互相转过 200 度以上，整个拓扑失效。
- **算例 B（只拧倾角）**：高度随手挑，倾角解锁定方程。面锁住了 ✓，
  但 24h 圈数是 15.055 / 14.275 / 13.160，残余相位 19.8° / 99.1° / 57.6° ✗。
- **算例 C（两个旋钮）**：面共动 ✓ + 相位复位 ✓ 同时成立。

### 集群级 API

```python
from cluster_solver import Shell, analyze_cluster, design_resonant_locked_cluster, format_report

design = design_resonant_locked_cluster(
    base_inclination_deg=53.0, common_period_hours=24.0,
    k_values=[15, 14, 13], planes=24, sats_per_plane=20)
print(format_report(analyze_cluster(design["shells"], common_period_hours=24.0)))
```

| 函数 | 作用 |
| --- | --- |
| `Shell(name, h, i, planes, sats_per_plane, phasing_f)` | Walker 记法 `i: t/p/f` 的壳层定义 |
| `design_resonant_locked_cluster(...)` | 两个旋钮一起拧，生成完整集群 |
| `design_locked_cluster(base, altitudes)` | 只拧倾角旋钮，给定高度求锁定倾角 |
| `analyze_cluster(shells, common_period_hours)` | 集群体检：共动 / 复位 / 可达 |
| `check_common_repeat(shells, T_hours)` | 各壳层残余相位与复位判定 |
| `pair_analysis(a, b)` | 差分进动、会合周期、最近间距、视距上限 |
| `intra_shell_geometry(shell)` | 同面 / 邻面间距与视距校验 |
| `synodic_period(a, b)` | 会合周期 $1/\|1/T_a - 1/T_b\|$ |
| `central_angle_between_planes(i, dRAAN, u)` | 邻面同相位星的地心张角 |
| `format_report(analysis)` | 输出可读的集群报告 |

其中邻面几何用的是闭式解：

$$\cos\theta(u) = \cos(\Delta\Omega)\left(\cos^2 u + \sin^2 u \cos^2 i\right) + \sin^2 u \sin^2 i$$

$u=0$（赤道）时 $\theta = \Delta\Omega$ 取最大，$u=90°$ 时两面收拢——
所以邻面激光链路的最长距离出现在赤道，链路预算按赤道值算。

## 密集编队可行性（`dense_formation.py`）

`cluster_solver` 管的是"多壳层宏观拓扑"；这个模块管另一个尺度：
**几十上百颗星挤在 1 km 半径内、彼此相距几百米**。

这个尺度上成败不取决于绝对轨道，而取决于**差分量**：

| 差分量 | 后果 | 量级（650 km） |
| --- | --- | --- |
| $\Delta a$ 半长轴失配 | 沿轨长期漂移 $3\pi\Delta a$ 每圈 | **1 m 失配 → 9.4 m/圈，51 km/年** |
| $\Delta i$ 倾角失配 | 差分 J2 节面进动 → 横向散开 | 1 km 横向若走 $\Delta i$ → 44.8 km/年 |
| $\Delta\Omega$ 升交点差 | $\dot\Omega$ 不依赖 $\Omega$ → **恒为零** | 0 km/年 |

### 两条设计结论

1. **半长轴必须匹配到亚米级**。漂移预算 100 m/天 → $\Delta a$ 容差仅 **72 cm**；
   要一年只漂 100 m，容差是 **2 mm**。这是密集编队真正的工程门槛。
2. **横向铺开要走 $\Delta\Omega$ 而不是 $\Delta i$**。因为 $\dot\Omega \propto a^{-3.5}\cos i$
   与 $\Omega$ 无关，同 $a$ 同 $i$ 只差 $\Omega$ 的集群，一阶 J2 长期漂移天然为零
   （J2-invariant），只剩差分气动阻力要管——而差分阻力可以靠姿态调制无推进剂修正。

### 参考算例：81 星 / 1 km 半径 / 650 km 晨昏 SSO

```
周期          : 97.73 min（14.735 圈/天）
节面进动      : +0.9856 °/天  = 360°/365.24 天 → 轨道面跟着太阳转
最近邻间距    : 212 m（六方密排，81 星铺满 1 km 半径圆盘）
2:1 相对椭圆  : e = 7.11e-05，径向 ±500 m，沿轨 ±1000 m，周期 97.73 min（闭合）
最近邻时延    : 0.706 μs 单程
接收功率增益  : 2.2e+07 倍（相对 1000 km 长程 ISL，1/d² 发散）
```

### API

| 函数 | 作用 |
| --- | --- |
| `analyze_dense_cluster(h, i, n, radius_m)` | 密集集群完整体检 |
| `along_track_drift(h, delta_a_m)` | $\Delta a$ → 每圈/每天/每年沿轨漂移 |
| `semi_major_axis_tolerance_m(h, budget, days)` | 反解 $\Delta a$ 容差 |
| `precession_sensitivity(h, i)` | $\partial\dot\Omega/\partial a$、$\partial\dot\Omega/\partial i$（对 $\Omega$ 恒为 0） |
| `differential_nodal_drift(h, i, da, di)` | 差分进动与横向散开 |
| `raan_offset_for_cross_track(...)` / `inclination_offset_for_cross_track(...)` | 两条横向路线对比 |
| `packing_spacing_m(n, radius_m)` | 六方/网格堆积的最近邻间距 |
| `light_time_us(d)` / `link_power_gain(far, near)` | 时延与 $1/d^2$ 功率增益 |

## 重复星下点与地面接触（`ground_track.py`）

`solve_resonant_shells` 算的是**惯性系**相位复位；但"什么时候能和地面通信"
由**星下点轨迹**决定，需要两个额外修正：

1. **交点周期**而非开普勒周期。J2 让 $M$ 与 $\omega$ 都长期漂移：
   $T_{\text{nodal}} = 2\pi/(\dot M + \dot\omega)$，650 km SSO 上比开普勒周期长 **7.25 s**
2. **交点日**而非恒星日。地球要追的是**进动中的**轨道面：
   $T_{\text{nodal day}} = 2\pi/(\omega_E - \dot\Omega)$
   —— 太阳同步轨道上此值精确等于一个平太阳日（86400 s）

两个修正合起来把 650 km 的"每天圈数"从 14.735 改成 **14.717**，足以选错重复周期。

### 关键结果：+723 m 把 I/O 变成严格周期

650 km 附近的可选重复周期：

| N圈/M天 | 圈/交点日 | 需要高度 | 相对 650 km | 轨迹经度间隔 |
|---|---|---|---|---|
| **103/7** | 14.71429 | **650.723 km** | **+0.723 km** | 24.466° |
| 265/18 | 14.72222 | 648.192 km | −1.808 km | 24.453° |
| 250/17 | 14.70588 | 653.407 km | +3.407 km | 24.480° |
| 147/10 | 14.70000 | 655.286 km | +5.286 km | 24.490° |

**高度抬高 723 米，地面接触模式就每 7 天精确复现一次。** 对算力集群，
这意味着 I/O 可用性从"近似周期"变成"严格周期"——调度可以离线预编译。

### 时间尺度谱（650 km 晨昏 SSO）

| 节拍 | 时长 | 支配 |
|---|---|---|
| 交点周期 / 星间相对椭圆 | 97.85 min | 链路距离 → 带宽（$1/d^2$） |
| 单次过顶窗口 | 9.05 min（仰角 10°） | 单次 I/O 窗口 |
| 交点日 | 24.00 h | 星下点经度推进 |
| 重复星下点 | 7 天 | **地面接触模式完全复现** |
| 太阳同步年周期 | 365.24 天 | beta 角 / 光照 / 热 |

单站占空比（纬度 60°，仰角门限 10°）：≈5.4 次/天、≈33 min/天、**2.3%**。

### API

| 函数 | 作用 |
| --- | --- |
| `nodal_period(a, i)` / `nodal_day(a, i)` | J2 交点周期与交点日 |
| `revs_per_nodal_day(h, i)` | 重复星下点判据的核心量 |
| `solve_repeat_altitude(N, M)` | 反解 N圈/M天 重复的高度（二分求根） |
| `find_repeat_options(h)` | 列出附近可用周期及所需高度微调 |
| `contact_window(h, mask)` | 单次过顶半张角与最长时长 |
| `passes_per_day_estimate(h, lat)` | 过顶次数与占空比（一阶估计） |
| `periodicity_spectrum(h)` | 全部时间尺度一览 |

## 模型边界（什么时候不能用）

- **共振壳层是惯性系周期共振**（全网相位复位），不是重复星下点轨迹；
  后者还需计入地球自转与 $J_2$ 对交点周期的修正。
- **无阻力模型**：低于 ~400 km 时半长轴实际会持续衰减，需要定期抬轨维持。
- **一级长期项**：只给 $\Omega$、$\omega$、$M$ 的长期漂移，不含短周期振荡（幅值约数 km）。
- **编队几何是线性化（CW/HCW）解**，仅在相对距离 ≪ 轨道半径时成立（百米～数十公里量级）。
- **集群分析只覆盖长期项**：会合周期、残余相位是一阶估计，不含短周期振荡与相位微调。
- **密集编队的差分模型是一阶的**：$\Delta a$ 漂移与差分 J2 进动都取长期项线性化，
  不含差分气动阻力（需要大气密度模型与面质比）、不含 J2 短周期项、不含碰撞规避。
  真要定 delta-v 预算必须上数值积分。
- **进动锁定不管相位**：面锁住只保证轨道面方位不散，同一面内的卫星仍按各自周期运行。
- **`passes_per_day_estimate` 只是量级估计**：把星下点看成等经度间隔铺开，
  忽略轨迹倾斜与高纬收敛，在站点纬度接近倾角时失效。要精确的接触窗口表
  必须做真实星历传播（SGP4 + 地面站可见性）。
- 需要米级绝对定轨、机动规划、碰撞预警时，请换用完整数值积分（SGP4/高精度力模型）。

## 校验依据

`test_orbit_solver.py` 用公开工程基准做交叉验证：星链 550 km/53° 壳的周期与进动、
800 km 太阳同步倾角 98.6°、恒星日反解 GEO 高度 35786 km、临界倾角 63.435° 近地点冻结、
以及最大视距连线到地心的最近距离恰好等于 $R_E + 80$ km（几何自洽性）。

`test_cluster_solver.py` 校验集群层：极轨面在极点收拢（$\theta \to 0$）、赤道处张角等于
$\Delta\Omega$、15 与 14 圈/天的会合周期恰为 86400 s、
以及核心断言 `test_resonant_design_satisfies_both_constraints`
——同一组壳层同时满足整数共振与进动锁定。
