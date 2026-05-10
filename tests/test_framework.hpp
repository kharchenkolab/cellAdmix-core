#pragma once

#include <cmath>
#include <exception>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace celladmix::test {

using TestFunction = std::function<void()>;

struct TestCase {
  std::string name;
  TestFunction fn;
};

inline std::vector<TestCase>& registry() {
  static std::vector<TestCase> tests;
  return tests;
}

struct Registrar {
  Registrar(std::string name, TestFunction fn) {
    registry().push_back(TestCase{std::move(name), std::move(fn)});
  }
};

inline std::string make_failure(
    const std::string& file,
    int line,
    const std::string& message) {
  std::ostringstream out;
  out << file << ":" << line << ": " << message;
  return out.str();
}

}  // namespace celladmix::test

#define CELLADMIX_CONCAT_IMPL(a, b) a##b
#define CELLADMIX_CONCAT(a, b) CELLADMIX_CONCAT_IMPL(a, b)

#define TEST_CASE(name)                                                           \
  static void CELLADMIX_CONCAT(test_fn_, __LINE__)();                             \
  static ::celladmix::test::Registrar CELLADMIX_CONCAT(test_reg_, __LINE__)(      \
      name,                                                                       \
      CELLADMIX_CONCAT(test_fn_, __LINE__));                                      \
  static void CELLADMIX_CONCAT(test_fn_, __LINE__)()

#define REQUIRE(condition)                                                        \
  do {                                                                            \
    if (!(condition)) {                                                           \
      throw std::runtime_error(::celladmix::test::make_failure(                   \
          __FILE__,                                                               \
          __LINE__,                                                               \
          std::string("Requirement failed: ") + #condition));                     \
    }                                                                             \
  } while (false)

#define REQUIRE_NEAR(lhs, rhs, tol)                                               \
  do {                                                                            \
    const auto celladmix_lhs = (lhs);                                             \
    const auto celladmix_rhs = (rhs);                                             \
    const auto celladmix_tol = (tol);                                             \
    if (std::fabs(celladmix_lhs - celladmix_rhs) > celladmix_tol) {               \
      std::ostringstream celladmix_msg;                                           \
      celladmix_msg << "Expected " << #lhs << " ~= " << #rhs << " within "        \
                    << celladmix_tol << " but got " << celladmix_lhs              \
                    << " and " << celladmix_rhs;                                  \
      throw std::runtime_error(::celladmix::test::make_failure(                   \
          __FILE__, __LINE__, celladmix_msg.str()));                              \
    }                                                                             \
  } while (false)

#define REQUIRE_GT(lhs, rhs) REQUIRE((lhs) > (rhs))
#define REQUIRE_GE(lhs, rhs) REQUIRE((lhs) >= (rhs))
#define REQUIRE_LT(lhs, rhs) REQUIRE((lhs) < (rhs))
#define REQUIRE_LE(lhs, rhs) REQUIRE((lhs) <= (rhs))
#define REQUIRE_EQ(lhs, rhs) REQUIRE((lhs) == (rhs))
