#include <string>
#include <vector>
#include <algorithm>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xadapt.hpp>


#include "pgmfactors/pgmfactors.hpp"

#include <fmt/core.h>

exported_class::exported_class()
    : m_name {fmt::format("{}", "pgmfactors")}
{
}

auto exported_class::name() const -> const char*
{
  return m_name.c_str();
}

namespace pgmfactors
{
// helper function to permute axes when they are specified out-of-order
// upon construction for a PFactor
std::vector<int> find_permutation_indices(const std::vector<int>& input)
{
  // Create a vector of pairs that stores the original index of each element
  std::vector<std::pair<int, int>> indexed_input;
  for (int i = 0; i < input.size(); i++) {
    indexed_input.emplace_back(input[i], i);
  }

  // Sort the vector of pairs by the element value
  std::sort(indexed_input.begin(),
            indexed_input.end(),
            [](const std::pair<int, int>& a, const std::pair<int, int>& b)
            { return a.first < b.first; });

  // Create a vector of permutation indices by mapping the sorted indices back
  // to the original indices
  std::vector<int> permutation_indices(input.size());
  for (int i = 0; i < input.size(); i++) {
    permutation_indices[indexed_input[i].second] = i;
  }

  return permutation_indices;
}


}  // namespace pgmfactors
