#include <iostream>
#include <string>

#include "test_framework.hpp"

int main() {
  int failures = 0;
  for (const auto& test : celladmix::test::registry()) {
    try {
      test.fn();
      std::cout << "[PASS] " << test.name << "\n";
    } catch (const std::exception& error) {
      failures += 1;
      std::cerr << "[FAIL] " << test.name << "\n" << error.what() << "\n";
    } catch (...) {
      failures += 1;
      std::cerr << "[FAIL] " << test.name << "\nUnknown error\n";
    }
  }
  if (failures > 0) {
    std::cerr << failures << " test(s) failed\n";
    return 1;
  }
  return 0;
}
