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

python3 -m unittest -v test_orbit_solver      # 38 项解析校验
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

## 模型边界（什么时候不能用）

- **共振壳层是惯性系周期共振**（全网相位复位），不是重复星下点轨迹；
  后者还需计入地球自转与 $J_2$ 对交点周期的修正。
- **无阻力模型**：低于 ~400 km 时半长轴实际会持续衰减，需要定期抬轨维持。
- **一级长期项**：只给 $\Omega$、$\omega$、$M$ 的长期漂移，不含短周期振荡（幅值约数 km）。
- **编队几何是线性化（CW/HCW）解**，仅在相对距离 ≪ 轨道半径时成立（百米～数十公里量级）。
- 需要米级绝对定轨、机动规划、碰撞预警时，请换用完整数值积分（SGP4/高精度力模型）。

## 校验依据

`test_orbit_solver.py` 用公开工程基准做交叉验证：星链 550 km/53° 壳的周期与进动、
800 km 太阳同步倾角 98.6°、恒星日反解 GEO 高度 35786 km、临界倾角 63.435° 近地点冻结、
以及最大视距连线到地心的最近距离恰好等于 $R_E + 80$ km（几何自洽性）。
