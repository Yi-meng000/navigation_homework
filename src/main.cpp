#include "navigation/coordinates.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace nav = navigation;

namespace {

void usage() {
  std::cout <<
      "Navigation coordinate converter (WGS-84)\n"
      "Angles: degrees; positions: metres. BLH order: latitude longitude height.\n"
      "Usage:\n"
      "  coordinate_demo                       Run examples\n"
      "  coordinate_demo blh2xyz B L H\n"
      "  coordinate_demo xyz2blh X Y Z\n"
      "  coordinate_demo xyz2neu B0 L0 H0 X Y Z\n"
      "  coordinate_demo neu2xyz B0 L0 H0 N E U\n"
      "  coordinate_demo blh2neu B0 L0 H0 B L H\n"
      "  coordinate_demo neu2blh B0 L0 H0 N E U\n"
      "  coordinate_demo body2neu yaw pitch roll forward right down\n"
      "  coordinate_demo neu2body yaw pitch roll N E U\n"
      "Body commands transform vectors; body axes: forward/right/down.\n"
      "Yaw: north toward east; pitch: nose up; roll: right wing down.\n";
}

double number(const char* argument) {
  const std::string text(argument);
  std::size_t consumed = 0;
  const double value = std::stod(text, &consumed);
  if (consumed != text.size()) throw std::invalid_argument("Invalid number: " + text);
  return value;
}

nav::Vector3 vector(char** arguments) {
  return {number(arguments[0]), number(arguments[1]), number(arguments[2])};
}

nav::Geodetic geodetic(char** arguments) {
  return {nav::deg_to_rad(number(arguments[0])),
          nav::deg_to_rad(number(arguments[1])), number(arguments[2])};
}

void print(const char* label, const nav::Vector3& value) {
  std::cout << label << " = " << value.x << "  " << value.y << "  " << value.z << '\n';
}

void print(const nav::Geodetic& value) {
  std::cout << "BLH (deg, deg, m) = " << nav::rad_to_deg(value.latitude_rad) << "  "
            << nav::rad_to_deg(value.longitude_rad) << "  " << value.height_m << '\n';
}

void examples() {
  const nav::Geodetic reference{nav::deg_to_rad(39.9), nav::deg_to_rad(116.4), 50.0};
  const nav::LocalFrame frame(reference);
  std::cout << "1. BLH <-> ECEF\n";
  print(reference);
  const auto xyz = nav::geodetic_to_ecef(reference);
  print("XYZ (m)", xyz);
  print(nav::ecef_to_geodetic(xyz));

  std::cout << "\n2. Fixed local NEU frame, reference BLH = (39.9 deg, 116.4 deg, 50 m)\n";
  const nav::Vector3 neu{100.0, 200.0, 10.0};
  print("NEU (m)", neu);
  const auto target_xyz = frame.neu_to_ecef(neu);
  print("XYZ (m)", target_xyz);
  const auto target_blh = frame.neu_to_geodetic(neu);
  print(target_blh);
  print("NEU recovered from XYZ (m)", frame.ecef_to_neu(target_xyz));
  print("NEU recovered from BLH (m)", frame.geodetic_to_neu(target_blh));

  std::cout << "\n3. Body FRD <-> NEU, yaw=90 deg, pitch=roll=0\n";
  const nav::Attitude attitude{nav::deg_to_rad(90.0), 0.0, 0.0};
  const nav::Vector3 body{10.0, 0.0, 0.0};
  print("Body (forward, right, down)", body);
  const auto geographic = nav::body_to_neu(body, attitude);
  print("NEU", geographic);
  print("Body recovered", nav::neu_to_body(geographic, attitude));
  std::cout << "\nUse --help for custom inputs.\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::cout << std::fixed << std::setprecision(9);
  try {
    if (argc == 1) {
      examples();
      return 0;
    }
    const std::string command = argv[1];
    if (command == "--help" || command == "-h") {
      usage();
      return 0;
    }
    if (command == "blh2xyz" || command == "xyz2blh") {
      if (argc != 5) throw std::invalid_argument("Expected 3 numbers; see --help");
      if (command == "blh2xyz") print("XYZ (m)", nav::geodetic_to_ecef(geodetic(argv + 2)));
      else print(nav::ecef_to_geodetic(vector(argv + 2)));
    } else if (command == "xyz2neu" || command == "neu2xyz" ||
               command == "blh2neu" || command == "neu2blh") {
      if (argc != 8) throw std::invalid_argument("Expected 6 numbers; see --help");
      const nav::LocalFrame frame(geodetic(argv + 2));
      if (command == "xyz2neu") print("NEU (m)", frame.ecef_to_neu(vector(argv + 5)));
      else if (command == "neu2xyz") print("XYZ (m)", frame.neu_to_ecef(vector(argv + 5)));
      else if (command == "blh2neu") print("NEU (m)", frame.geodetic_to_neu(geodetic(argv + 5)));
      else print(frame.neu_to_geodetic(vector(argv + 5)));
    } else if (command == "body2neu" || command == "neu2body") {
      if (argc != 8) throw std::invalid_argument("Expected 6 numbers; see --help");
      const nav::Attitude attitude{nav::deg_to_rad(number(argv[2])),
                                   nav::deg_to_rad(number(argv[3])),
                                   nav::deg_to_rad(number(argv[4]))};
      if (command == "body2neu") print("NEU", nav::body_to_neu(vector(argv + 5), attitude));
      else print("Body (forward, right, down)", nav::neu_to_body(vector(argv + 5), attitude));
    } else {
      throw std::invalid_argument("Unknown command: " + command + "; see --help");
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << '\n';
    return 1;
  }
}
