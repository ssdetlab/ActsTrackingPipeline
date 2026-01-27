#include <cmath>
#include <deque>

#include "TrackingPipeline/Preprocessing/AnalysisFunctions.hpp"

void bfsClustering(std::vector<Hit> &clusterPixels, Hit &pivot,
                   HitSet &pixels) {
  std::deque<Hit> queue{pivot};
  while (!queue.empty()) {
    auto current = queue.front();
    queue.pop_front();
    auto currentPtr = pixels.find(current);

    if (currentPtr != pixels.end()) {
      // add this pixel to the cluster
      clusterPixels.push_back(current);
      // remove pixel from live pixels list
      pixels.erase(currentPtr);

      std::vector<Hit> neighbours{{current.first + 1, current.second},
                                  {current.first - 1, current.second},
                                  {current.first, current.second + 1},
                                  {current.first, current.second - 1}};
      for (const auto &n : neighbours) {
        auto nptr = pixels.find(n);
        if (nptr != pixels.end()) {
          queue.push_back(*nptr);
        }
      }
    }
  }
}

std::vector<Cluster> getClusters(HitSet &hits) {
  std::vector<Cluster> clusters;

  while (!hits.empty()) {
    auto pixel = *hits.begin();
    std::vector<Hit> clusterPixels;
    bfsClustering(clusterPixels, pixel, hits);

    std::size_t xMin = 1024;
    std::size_t xMax = 0;
    std::size_t yMin = 512;
    std::size_t yMax = 0;
    double xCenter = 0;
    double yCenter = 0;
    for (auto pix : clusterPixels) {
      xCenter += pix.first;
      yCenter += pix.second;

      if (pix.first <= xMin) {
        xMin = pix.first;
      }
      if (pix.first >= xMax) {
        xMax = pix.first;
      }
      if (pix.second <= yMin) {
        yMin = pix.second;
      }
      if (pix.second >= yMax) {
        yMax = pix.second;
      }
    }
    std::size_t dx = xMax - xMin + 1;
    std::size_t dy = yMax - yMin + 1;
    std::size_t cl_size = clusterPixels.size();
    xCenter /= clusterPixels.size();
    yCenter /= clusterPixels.size();

    clusters.push_back(Cluster{xCenter, yCenter, dx, dy, cl_size});
  }

  return clusters;
}
