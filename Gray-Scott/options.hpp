#pragma once

#include <string>
#include <cstring>

template <typename Real = double, typename Index = std::int64_t>
struct Options {
  private:
    template <typename T>
    [[nodiscard]] constexpr static T parse(const char * const id, int argc, const char * const argv[], T def) noexcept {
      for (int i = 0; i < argc; ++i) {
        if (std::strcmp(argv[i], id) == 0 && i + 1 < argc) {
          if constexpr (std::is_same_v<T, Real>) {
            return std::stold(argv[i + 1]);
          } else if constexpr (std::is_same_v<T, int>) {
            return std::stoll(argv[i + 1]);
          } else if constexpr (std::is_same_v<T, bool>) {
            return std::strcmp(argv[i + 1], "true") == 0;
          } else {
            return argv[i + 1];
          }
        }
      }
      return def;
    }
  
  public:
    const std::string method;
    const std::string out_file;
    const Index grid_pts_1d;
    const int order;
    const Real dt;
    const Real abs_tol;
    const Real rel_tol;
    const bool low_storage;
  
    constexpr Options(int argc, const char * const argv[]) noexcept :
      method(parse("method", argc, argv, "ERK")),
      out_file(parse("out_file", argc, argv, "solution.txt")),
      grid_pts_1d(parse("grid_pts_1d", argc, argv, 128)),
      order(parse("order", argc, argv, 2)),
      dt(parse("dt", argc, argv, Real(0))),
      abs_tol(parse("abs_tol", argc, argv, Real(1e-12L))),
      rel_tol(parse("rel_tol", argc, argv, Real(1e-4L))),
      low_storage(parse("low_storage", argc, argv, false))
      {}
  };