#pragma once

#include <array>

namespace navigation {

inline constexpr double kPi = 3.141592653589793238462643383279502884;
inline constexpr double kWgs84A = 6378137.0;  // 长半轴，m
inline constexpr double kWgs84InverseFlattening = 298.257223563;
inline constexpr double kWgs84F = 1.0 / kWgs84InverseFlattening;
inline constexpr double kWgs84B = kWgs84A * (1.0 - kWgs84F);
inline constexpr double kWgs84E2 = kWgs84F * (2.0 - kWgs84F);

constexpr double deg_to_rad(double degrees) { return degrees * (kPi / 180.0); }
constexpr double rad_to_deg(double radians) { return radians * (180.0 / kPi); }

// 列向量。ECEF: X/Y/Z；NEU: 北/东/天；载体系: 前/右/下。
struct Vector3 {
  double x;
  double y;
  double z;
};

struct Geodetic {
  double latitude_rad;   // B，大地纬度，[-pi/2, pi/2]
  double longitude_rad;  // L，经度，[-pi, pi]，东经为正
  double height_m;       // H，椭球高，不是海拔
};

// 载体相对于当地 NED 的内禀 Z-Y-X 姿态：航向、俯仰、横滚。
// yaw: 从北向东；pitch: 抬头；roll: 右翼下沉。全部为 rad。
struct Attitude {
  double yaw_rad;
  double pitch_rad;
  double roll_rad;
};

using Matrix3 = std::array<std::array<double, 3>, 3>;

Vector3 geodetic_to_ecef(const Geodetic& blh);
// 极轴上约定 L=0；地心抛出 invalid_argument；迭代失败抛出 runtime_error。
// 面向地表附近和地表外导航；深地内部的多值大地坐标不作唯一性保证。
Geodetic ecef_to_geodetic(const Vector3& xyz);

// 所有位置使用同一个固定参考点的切平面；构造后可重复转换。
class LocalFrame {
 public:
  explicit LocalFrame(const Geodetic& reference);

  Vector3 ecef_to_neu(const Vector3& point_ecef) const;
  Vector3 neu_to_ecef(const Vector3& point_neu) const;
  Vector3 geodetic_to_neu(const Geodetic& point_blh) const;
  Geodetic neu_to_geodetic(const Vector3& point_neu) const;

  // 自由向量只变换分量，不作原点平移；例如速度、位移、加速度。
  Vector3 ecef_vector_to_neu(const Vector3& vector_ecef) const;
  Vector3 neu_vector_to_ecef(const Vector3& vector_neu) const;

  const Vector3& reference_ecef() const { return reference_ecef_; }
  const Matrix3& ecef_to_neu_matrix() const { return ecef_to_neu_; }

 private:
  Vector3 reference_ecef_;
  Matrix3 ecef_to_neu_;
};

// NEU 是左手系；矩阵为 diag(1,1,-1) * Rz(yaw)*Ry(pitch)*Rx(roll)。
Matrix3 body_to_neu_matrix(const Attitude& attitude);
Vector3 body_to_neu(const Vector3& vector_body, const Attitude& attitude);
Vector3 neu_to_body(const Vector3& vector_neu, const Attitude& attitude);

// 位置转换另需载体原点在该 NEU 系中的坐标。
Vector3 body_point_to_neu(const Vector3& point_body, const Attitude& attitude,
                          const Vector3& body_origin_neu);
Vector3 neu_point_to_body(const Vector3& point_neu, const Attitude& attitude,
                          const Vector3& body_origin_neu);

}  // namespace navigation
