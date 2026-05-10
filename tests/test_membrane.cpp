#include <vector>

#include "celladmix/membrane.hpp"
#include "test_framework.hpp"

using namespace celladmix;

TEST_CASE("Membrane score prefers points whose centroid line crosses a bright membrane") {
  Image2D image;
  image.width = 100;
  image.height = 100;
  image.pixel_size = 1.0;
  image.values.assign(100 * 100, 0.05);
  for (int y = 0; y < image.height; ++y) {
    image.values[static_cast<std::size_t>(y * image.width + 50)] = 1.0;
    image.values[static_cast<std::size_t>(y * image.width + 51)] = 1.0;
  }

  const std::vector<std::pair<double, double>> real_points = {
      {49.0, 20.0},
      {49.5, 40.0},
      {48.5, 60.0},
  };
  const std::vector<std::pair<double, double>> control_points = {
      {70.0, 20.0},
      {72.0, 40.0},
      {74.0, 60.0},
  };

  const double score = membrane_score(image, real_points, control_points, {80.0, 50.0}, 0.5);
  REQUIRE_GT(score, 1.0);
}
