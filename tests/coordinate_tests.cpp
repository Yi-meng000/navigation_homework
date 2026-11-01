#include "navigation/coordinates.hpp"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace nav = navigation;

namespace {
int checks = 0;

void near(double actual, double expected, double tolerance, const std::string& name) {
  ++checks;
  if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
    throw std::runtime_error(name + ": actual=" + std::to_string(actual) +
                             ", expected=" + std::to_string(expected));
  }
}

void near(const nav::Vector3& actual, const nav::Vector3& expected,
          double tolerance, const std::string& name) {
  near(actual.x, expected.x, tolerance, name + "[0]");
  near(actual.y, expected.y, tolerance, name + "[1]");
  near(actual.z, expected.z, tolerance, name + "[2]");
}

template <class F>
void invalid(F function, const std::string& name) {
  ++checks;
  try { function(); }
  catch (const std::invalid_argument&) { return; }
  throw std::runtime_error(name + ": expected invalid_argument");
}

void orthogonal_left_handed(const nav::Matrix3& m) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double product = 0.0;
      for (int k = 0; k < 3; ++k) product += m[i][k] * m[j][k];
      near(product, i == j ? 1.0 : 0.0, 2e-15, "matrix orthogonality");
    }
  }
  const double determinant =
      m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
      m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
      m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  near(determinant, -1.0, 2e-15, "NEU handedness");
}

void known_points() {
  near(nav::geodetic_to_ecef({0, 0, 0}), {6378137.0, 0, 0}, 1e-9, "equator Greenwich");
  near(nav::geodetic_to_ecef({0, nav::kPi / 2, 100}), {0, 6378237.0, 0}, 1e-9, "90E");
  near(nav::geodetic_to_ecef({nav::kPi / 2, 1.2, 0}),
       {0, 0, 6356752.314245179}, 1e-8, "north pole");
  // 独立参考值：PROJ cct +proj=cart +ellps=WGS84，角度输入顺序 L B H。
  const nav::Geodetic beijing{nav::deg_to_rad(39.9), nav::deg_to_rad(116.4), 50};
  const nav::Vector3 beijing_xyz{-2178657.082724949, 4388876.233551471, 4069505.747981717};
  near(nav::geodetic_to_ecef(beijing), beijing_xyz, 2e-8, "PROJ Beijing forward");
  const auto beijing_back = nav::ecef_to_geodetic(beijing_xyz);
  near(beijing_back.latitude_rad, beijing.latitude_rad, 2e-14, "PROJ Beijing latitude");
  near(beijing_back.longitude_rad, beijing.longitude_rad, 2e-14, "PROJ Beijing longitude");
  near(beijing_back.height_m, beijing.height_m, 2e-8, "PROJ Beijing height");
  near(nav::geodetic_to_ecef({nav::deg_to_rad(42), nav::deg_to_rad(-82), 200}),
       {660675.251824734, -4700948.683162266, 4245737.662222386}, 2e-8,
       "PROJ western hemisphere");
  // PROJ topocentric 输出 E/N/U，这里按 N/E/U 排列。
  const nav::LocalFrame beijing_frame(beijing);
  const nav::Geodetic target{nav::deg_to_rad(39.901), nav::deg_to_rad(116.402), 60};
  const nav::Vector3 target_neu{111.035683081, 171.035707270, 9.996740954};
  near(beijing_frame.geodetic_to_neu(target), target_neu, 2e-8, "PROJ local forward");
  const auto target_back = beijing_frame.neu_to_geodetic(target_neu);
  near(target_back.latitude_rad, target.latitude_rad, 2e-14, "PROJ local latitude");
  near(target_back.longitude_rad, target.longitude_rad, 2e-14, "PROJ local longitude");
  near(target_back.height_m, target.height_m, 2e-8, "PROJ local height");
  for (const double sign : {-1.0, 1.0}) {
    const auto pole = nav::ecef_to_geodetic({0, 0, sign * (nav::kWgs84B + 123.0)});
    near(pole.latitude_rad, sign * nav::kPi / 2, 0, "pole latitude");
    near(pole.longitude_rad, 0, 0, "pole canonical longitude");
    near(pole.height_m, 123.0, 1e-8, "pole height");
  }
  // ESA Navipedia 的轴定义：赤道零经度处 X=天、Y=东、Z=北。
  const nav::LocalFrame frame({0, 0, 0});
  near(frame.ecef_to_neu({nav::kWgs84A + 3, 2, 1}), {1, 2, 3}, 1e-9, "local axes");
  near(frame.neu_to_ecef({1, 2, 3}), {nav::kWgs84A + 3, 2, 1}, 1e-9, "local inverse");
  near(frame.ecef_vector_to_neu({3, 2, 1}), {1, 2, 3}, 1e-12, "vector has no translation");
  const nav::Attitude zero{0, 0, 0};
  near(nav::body_to_neu({1, 2, 3}, zero), {1, 2, -3}, 1e-12, "level north FRD");
  near(nav::body_to_neu({1, 0, 0}, {nav::kPi / 2, 0, 0}), {0, 1, 0}, 1e-12, "yaw east");
  near(nav::body_to_neu({0, 1, 0}, {nav::kPi / 2, 0, 0}), {-1, 0, 0}, 1e-12, "right points south");
  near(nav::body_to_neu({1, 0, 0}, {0, nav::kPi / 2, 0}), {0, 0, 1}, 1e-12, "pitch nose up");
  near(nav::body_to_neu({0, 1, 0}, {0, 0, nav::kPi / 2}), {0, 0, -1}, 1e-12, "roll right down");
  near(nav::body_point_to_neu({1, 2, 3}, zero, {10, 20, 30}), {11, 22, 27}, 1e-12,
       "body origin translation");
  near(nav::neu_point_to_body({11, 22, 27}, zero, {10, 20, 30}), {1, 2, 3}, 1e-12,
       "body origin inverse");
}

void round_trips() {
  // 覆盖南北半球、日期变更线、负高程、近极区及卫星高度。
  for (double latitude : {-90.0, -89.999999, -45.0, 0.0, 39.9, 89.999999, 90.0}) {
    for (double longitude : {-180.0, -116.4, 0.0, 116.4, 180.0}) {
      for (double height : {-1000.0, 0.0, 50.0, 1e6, 35786000.0}) {
        const nav::Geodetic g{nav::deg_to_rad(latitude), nav::deg_to_rad(longitude), height};
        const auto xyz = nav::geodetic_to_ecef(g);
        const auto back = nav::ecef_to_geodetic(xyz);
        near(back.latitude_rad, g.latitude_rad, 2e-13, "BLH round trip latitude");
        if (std::abs(latitude) != 90.0) {
          near(std::remainder(back.longitude_rad - g.longitude_rad, 2 * nav::kPi),
               0.0, 2e-14, "BLH round trip longitude");
        }
        near(back.height_m, height, 3e-8, "BLH round trip height");
        near(nav::geodetic_to_ecef(back), xyz, 5e-8, "XYZ round trip");
        const nav::LocalFrame frame(g);
        near(frame.geodetic_to_neu(g), {0, 0, 0}, 1e-9, "reference maps to origin");
        const nav::Vector3 local{1000.0, -2000.0, 30.0};
        near(frame.ecef_to_neu(frame.neu_to_ecef(local)), local, 1e-8, "NEU round trip");
        near(frame.geodetic_to_neu(frame.neu_to_geodetic(local)), local, 5e-8,
             "NEU BLH round trip");
        orthogonal_left_handed(frame.ecef_to_neu_matrix());
      }
    }
  }

  std::mt19937 engine(20261101);
  std::uniform_real_distribution<double> angle(-nav::kPi, nav::kPi);
  std::uniform_real_distribution<double> component(-1000.0, 1000.0);
  for (int i = 0; i < 1000; ++i) {
    const nav::Attitude attitude{angle(engine), angle(engine), angle(engine)};
    const nav::Vector3 body{component(engine), component(engine), component(engine)};
    near(nav::neu_to_body(nav::body_to_neu(body, attitude), attitude), body, 2e-12,
         "body vector round trip");
    orthogonal_left_handed(nav::body_to_neu_matrix(attitude));
    const nav::Vector3 origin{100, -200, 30};
    near(nav::neu_point_to_body(nav::body_point_to_neu(body, attitude, origin), attitude, origin),
         body, 2e-12, "body point round trip");
  }
}

void invalid_inputs() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  invalid([] { nav::ecef_to_geodetic({0, 0, 0}); }, "Earth center");
  invalid([] { nav::geodetic_to_ecef({nav::deg_to_rad(91), 0, 0}); }, "invalid latitude");
  invalid([] { nav::LocalFrame frame({0, nav::deg_to_rad(181), 0}); }, "invalid longitude");
  invalid([&] { nav::geodetic_to_ecef({0, 0, nan}); }, "NaN height");
  invalid([&] { nav::ecef_to_geodetic({inf, 1, 1}); }, "infinite XYZ");
  invalid([&] { nav::body_to_neu({1, 2, 3}, {nan, 0, 0}); }, "NaN attitude");
  invalid([&] { nav::neu_to_body({nan, 0, 0}, {0, 0, 0}); }, "NaN vector");
  invalid([&] { nav::body_point_to_neu({0, 0, 0}, {0, 0, 0}, {nan, 0, 0}); }, "NaN origin");
}

}  // namespace

int main() {
  try {
    known_points();
    round_trips();
    invalid_inputs();
    std::cout << "Passed " << checks << " checks\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
