# 导航常用坐标系转换（C++17）

无第三方库依赖，包含可复用库、命令行程序和自动测试。

| 功能 | 接口 |
| --- | --- |
| WGS-84 大地坐标 BLH ↔ 地心地固直角坐标 ECEF XYZ | `geodetic_to_ecef` / `ecef_to_geodetic` |
| ECEF XYZ ↔ 固定参考点的当地 NEU | `LocalFrame::ecef_to_neu` / `neu_to_ecef` |
| BLH ↔ 固定参考点的当地 NEU | `LocalFrame::geodetic_to_neu` / `neu_to_geodetic` |
| 载体系向量 ↔ 地理系 NEU 向量 | `body_to_neu` / `neu_to_body` |
| 带载体原点平移的位置转换 | `body_point_to_neu` / `neu_point_to_body` |

## 1. 编译与运行

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
ctest --test-dir build --output-on-failure
./build/coordinate_demo
```

直接运行显示三类转换的示例；`--help` 查看命令。程序可用于自定义输入：

```bash
# B=39.9°，L=116.4°，H=50 m → XYZ
./build/coordinate_demo blh2xyz 39.9 116.4 50
# 输出 XYZ (m) = -2178657.082724949  4388876.233551471  4069505.747981717

./build/coordinate_demo xyz2blh -2178657.082724949 4388876.233551471 4069505.747981717

# 前三个数指定 NEU 参考点 B0 L0 H0，后三个数指定目标点
./build/coordinate_demo blh2neu 39.9 116.4 50 39.901 116.402 60
./build/coordinate_demo xyz2neu 39.9 116.4 50 -2178657.082724949 4388876.233551471 4069505.747981717
./build/coordinate_demo neu2xyz 39.9 116.4 50 100 200 10
./build/coordinate_demo neu2blh 39.9 116.4 50 100 200 10

# yaw pitch roll = 90° 0° 0°，载体前向向量 (10,0,0) → 朝东 (0,10,0)
./build/coordinate_demo body2neu 90 0 0 10 0 0
./build/coordinate_demo neu2body 90 0 0 0 10 0
```

命令行角度均为**度**，库接口角度均为**弧度**。位置和高程单位为 m。输入错误时程序返回非零退出码并输出原因。

## 2. 坐标与姿态约定

- **BLH**：B 为大地纬度（北正，范围 ±90°），L 为经度（东正，范围 ±180°），H 为椭球高。B 不是地心纬度；H 不是海拔。
- **ECEF XYZ**：原点在地球质心，X 指向赤道与零经线交点，Y 指向赤道与东经 90° 交点，Z 指向北极，构成右手系。这里的地球坐标系特指地心地固系。
- **NEU**：原点为指定参考点，三个分量依次是北、东、天；天沿参考点的椭球法线。NEU 是左手系。
- **载体系**：采用前、右、下（FRD）右手系。零姿态时前向北、右向东、下向地。
- **姿态**：航向 yaw 从北向东为正；俯仰 pitch 抬头为正；横滚 roll 右翼下沉为正。采用相对于 NED 的内禀 Z-Y-X 欧拉角定义；所有向量按列向量计算。

如果设备使用前左上、ENU、不同姿态角定义，需要先按设备接口统一轴序和符号。

## 3. WGS-84 BLH ↔ XYZ

椭球参数：

```text
a   = 6378137.0 m                 长半轴
1/f = 298.257223563               反扁率
b   = a(1-f)                     短半轴 ≈ 6356752.314245179 m
e²  = f(2-f)                     第一偏心率平方 ≈ 0.0066943799901413165
```

参数来源：[NGA WGS 84](https://earth-info.nga.mil/index.php?action=wgs84&dir=wgs84)。

正转换（Nφ 为卯酉圈曲率半径，区别于 NEU 的北向分量）：

```text
Nφ = a / sqrt(1 - e² sin²B)
X  = (Nφ + H) cosB cosL
Y  = (Nφ + H) cosB sinL
Z  = [Nφ(1 - e²) + H] sinB
```

反转换使用迭代法：

```text
p    = hypot(X, Y)
L    = atan2(Y, X)
B₀   = atan2(Z, p(1-e²))
Nₖ   = a / sqrt(1 - e² sin²Bₖ)
Bₖ₊₁ = atan2(Z + e² Nₖ sinBₖ, p)
```

当相邻纬度差小于 `1e-14 rad` 时收敛，最多迭代 100 次，未收敛抛出异常。高程用稳定的法线投影形式计算，避免极区除以 `cosB`：

```text
H = p cosB + Z sinB - a sqrt(1 - e² sin²B)
```

极轴 `X=Y=0, Z≠0`：`B=sign(Z)·90°`，`H=|Z|-b`，经度约定输出 0°。此处经度没有唯一值；构造极点 NEU 时仍可指定参考经度来定义水平轴。地心没有唯一 BLH，程序拒绝转换。

反转换面向地表附近和地表外导航；深地内部可能存在多个大地坐标解，不保证恢复原始 BLH 分支。自动测试覆盖 H=-1000 m 至 35786000 m。

## 4. XYZ / BLH ↔ NEU

给定参考点 `BLH₀=(B₀,L₀,H₀)`，先求它的 ECEF 位置 `r₀`。同一个 `LocalFrame` 始终使用这个固定参考点。

```text
          [ -sinB₀ cosL₀   -sinB₀ sinL₀    cosB₀ ]
Cₑ→ₙ =   [ -sinL₀          cosL₀          0     ]
          [  cosB₀ cosL₀    cosB₀ sinL₀    sinB₀ ]

[N E U]ᵀ = Cₑ→ₙ (r - r₀)
r         = r₀ + Cₑ→ₙᵀ [N E U]ᵀ
```

BLH 到 NEU：`BLH → XYZ → NEU`；反之：`NEU → XYZ → BLH`。这里转换的是空间直角坐标分量，N/E 不等于椭球面上的测地距离，U 也通常不等于两点椭球高之差。

**位置需要减去参考原点，自由向量只变换分量。** 因此库另外提供 `ecef_vector_to_neu` / `neu_vector_to_ecef`，用于在同一时刻表示速度、位移或加速度分量。移动参考系的速度/加速度导数关系还需要相应运动项，不由这组分量转换接口计算。

公式参考：[ESA Navipedia：ECEF 与 ENU 转换](https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates)。本程序将其中的 E/N 两行交换，得到 NEU。

## 5. 载体系 ↔ NEU

令 `ψ=yaw`、`θ=pitch`、`φ=roll`，`cψ=cosψ`、`sψ=sinψ`，其他角同理。先建立载体系到右手 NED 的姿态矩阵，再翻转向下轴：

```text
C_b→NEU = diag(1,1,-1) · Rz(ψ) · Ry(θ) · Rx(φ)

          [ cθ cψ    sφ sθ cψ - cφ sψ    cφ sθ cψ + sφ sψ ]
        = [ cθ sψ    sφ sθ sψ + cφ cψ    cφ sθ sψ - sφ cψ ]
          [ sθ      -sφ cθ               -cφ cθ             ]

v_NEU = C_b→NEU v_body
v_body = C_b→NEUᵀ v_NEU
```

由于 FRD 与 NEU 手性不同，这个矩阵的行列式为 -1，仍为正交矩阵，逆矩阵等于转置。零姿态转换为 `(N,E,U)=(前,右,-下)`。程序只通过给定欧拉角计算矩阵，即使 pitch=±90° 也能转换向量；不提供存在欧拉角奇异性的矩阵反求姿态接口。

若转换的是点的位置，且两系原点不同：

```text
p_NEU  = t_NEU + C_b→NEU p_body
p_body = C_b→NEUᵀ (p_NEU - t_NEU)
```

`t_NEU` 是载体原点在目标 NEU 中的位置，由调用者提供。姿态必须相对于**目标 NEU 的轴**定义。若姿态来自载体当前位置的当地水平系，而目标使用另一个固定参考点的 NEU，需要先转换两个当地系之间的方向，再组合变换。

