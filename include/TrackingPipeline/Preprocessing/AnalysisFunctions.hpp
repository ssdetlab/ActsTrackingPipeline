#pragma once

#include <unordered_set>
#include <vector>

#include "DetectorEvent.hpp"
#include "PairHash.hpp"

using Cluster = ApollonIo::ChipEvent::Cluster;
using Hit = std::pair<std::uint16_t, std::uint16_t>;
using HitSet = std::unordered_set<Hit, PairHash>;

void bfsClustering(std::vector<Hit> &clusterPixels, Hit &pivot, HitSet &pixels);

std::vector<Cluster> getClusters(HitSet &hits);
