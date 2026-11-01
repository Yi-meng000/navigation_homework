#include "navigation/coordinates.hpp"

#include <cmath>
#include <stdexcept>

namespace navigation {
namespace {

void validate(const Vector3& v) {
  if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z)) {
    throw std::invalid_argument("Vector components must be finite");
  }
}

void validate(const Geodetic& g) {
  if (!std::isfinite(g.latitude_rad) || !std::isfinite(g.longitude_rad) ||
      !std::isfinite(g.height_m)) {
    throw std::invalid_argument("BLH components must be finite");
  }
  if (std::abs(g.latitude_rad) > kPi / 2.0 || std::abs(g.longitude_rad) > kPi) {
    throw std::invalid_argument("B must be in [-90,90], L in [-180,180] degrees");
  }
}

Vector3 add(const Vector3& a, const Vector3& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vector3 subtract(const Vector3& a, const Vector3& b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Matrix3 transpose(const Matrix3& m) {
  Matrix3 result{};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) result[i][j] = m[j][i];
  }
  return result;
}

Vector3 multiply(const Matrix3& m, const Vector3& v) {
  validate(v);
  return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
          m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
          m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
}

}  // namespace

Vector3 geodetic_to_ecef(const Geodetic& blh) {
  validate(blh);
  const double sb = std::sin(blh.latitude_rad);
  // 精确极点直接设置 cos(B)=0，避免生成伪水平分量。
  const double cb = std::abs(blh.latitude_rad) == kPi / 2.0
                        ? 0.0 : std::cos(blh.latitude_rad);
  const double n = kWgs84A / std::sqrt(1.0 - kWgs84E2 * sb * sb);
  return {(n + blh.height_m) * cb * std::cos(blh.longitude_rad),
          (n + blh.height_m) * cb * std::sin(blh.longitude_rad),
          (n * (1.0 - kWgs84E2) + blh.height_m) * sb};
}

Geodetic ecef_to_geodetic(const Vector3& xyz) {
  validate(xyz);
  const double p = std::hypot(xyz.x, xyz.y);
  if (!std::isfinite(std::hypot(p, xyz.z))) {
    throw std::invalid_argument("ECEF magnitude is outside double range");
  }
  if (p == 0.0) {
    if (xyz.z == 0.0) throw std::invalid_argument("BLH is undefined at Earth center");
    return {std::copysign(kPi / 2.0, xyz.z), 0.0, std::abs(xyz.z) - kWgs84B};
  }

  const double longitude = std::atan2(xyz.y, xyz.x);
  double latitude = std::atan2(xyz.z, p * (1.0 - kWgs84E2));
  // 由 z = (N*(1-e²)+H)*sin(B) 和 p = (N+H)*cos(B) 消去 H。
  for (int iteration = 0; iteration < 100; ++iteration) {
    const double sb = std::sin(latitude);
    const double n = kWgs84A / std::sqrt(1.0 - kWgs84E2 * sb * sb);
    const double next = std::atan2(xyz.z + kWgs84E2 * n * sb, p);
    if (std::abs(next - latitude) < 1e-14) {
      latitude = next;
      const double s = std::sin(latitude);
      // 投影形式在赤道和极区均稳定，避免 p/cos(B)-N 的极区除零。
      const double height = p * std::cos(latitude) + xyz.z * s -
                            kWgs84A * std::sqrt(1.0 - kWgs84E2 * s * s);
      return {latitude, longitude, height};
    }
    latitude = next;
  }
  throw std::runtime_error("ECEF to BLH did not converge within 100 iterations");
}

LocalFrame::LocalFrame(const Geodetic& reference)
    : reference_ecef_(geodetic_to_ecef(reference)) {
  const double sb = std::sin(reference.latitude_rad);
  const double cb = std::cos(reference.latitude_rad);
  const double sl = std::sin(reference.longitude_rad);
  const double cl = std::cos(reference.longitude_rad);
  ecef_to_neu_ = {{{-sb * cl, -sb * sl, cb},
                   {-sl, cl, 0.0},
                   {cb * cl, cb * sl, sb}}};
}

Vector3 LocalFrame::ecef_to_neu(const Vector3& point_ecef) const {
  validate(point_ecef);
  return ecef_vector_to_neu(subtract(point_ecef, reference_ecef_));
}

Vector3 LocalFrame::neu_to_ecef(const Vector3& point_neu) const {
  return add(reference_ecef_, neu_vector_to_ecef(point_neu));
}

Vector3 LocalFrame::geodetic_to_neu(const Geodetic& point_blh) const {
  return ecef_to_neu(geodetic_to_ecef(point_blh));
}

Geodetic LocalFrame::neu_to_geodetic(const Vector3& point_neu) const {
  return ecef_to_geodetic(neu_to_ecef(point_neu));
}

Vector3 LocalFrame::ecef_vector_to_neu(const Vector3& vector_ecef) const {
  return multiply(ecef_to_neu_, vector_ecef);
}

Vector3 LocalFrame::neu_vector_to_ecef(const Vector3& vector_neu) const {
  return multiply(transpose(ecef_to_neu_), vector_neu);
}

Matrix3 body_to_neu_matrix(const Attitude& attitude) {
  if (!std::isfinite(attitude.yaw_rad) || !std::isfinite(attitude.pitch_rad) ||
      !std::isfinite(attitude.roll_rad)) {
    throw std::invalid_argument("Attitude angles must be finite");
  }
  const double cy = std::cos(attitude.yaw_rad), sy = std::sin(attitude.yaw_rad);
  const double cp = std::cos(attitude.pitch_rad), sp = std::sin(attitude.pitch_rad);
  const double cr = std::cos(attitude.roll_rad), sr = std::sin(attitude.roll_rad);
  // 先将前右下载体系转换到 NED，再翻转第三行得到 NEU。
  return {{{cp * cy, sr * sp * cy - cr * sy, cr * sp * cy + sr * sy},
           {cp * sy, sr * sp * sy + cr * cy, cr * sp * sy - sr * cy},
           {sp, -sr * cp, -cr * cp}}};
}

Vector3 body_to_neu(const Vector3& vector_body, const Attitude& attitude) {
  return multiply(body_to_neu_matrix(attitude), vector_body);
}

Vector3 neu_to_body(const Vector3& vector_neu, const Attitude& attitude) {
  return multiply(transpose(body_to_neu_matrix(attitude)), vector_neu);
}

Vector3 body_point_to_neu(const Vector3& point_body, const Attitude& attitude,
                          const Vector3& body_origin_neu) {
  validate(body_origin_neu);
  return add(body_origin_neu, body_to_neu(point_body, attitude));
}

Vector3 neu_point_to_body(const Vector3& point_neu, const Attitude& attitude,
                          const Vector3& body_origin_neu) {
  validate(point_neu);
  validate(body_origin_neu);
  return neu_to_body(subtract(point_neu, body_origin_neu), attitude);
}

}  // namespace navigation
