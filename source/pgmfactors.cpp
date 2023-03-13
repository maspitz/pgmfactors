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
std::vector<int> find_permutation_indices(factor::variable_list& rvs)
{
  // Create a vector of pairs that stores the original index of each element
  std::vector<std::pair<int, int>> indexed_input;
  for (int i = 0; i < rvs.size(); i++) {
    indexed_input.emplace_back(rvs[i].id(), i);
  }

  // Sort the vector of pairs by the element value
  std::sort(indexed_input.begin(),
            indexed_input.end(),
            [](const std::pair<int, int>& a, const std::pair<int, int>& b)
            { return a.first < b.first; });

  // Create a vector of permutation indices by mapping the sorted indices back
  // to the original indices
  std::vector<int> permutation_indices(rvs.size());
  for (int i = 0; i < rvs.size(); i++) {
    permutation_indices[indexed_input[i].second] = i;
  }

  return permutation_indices;
}



factor::factor(variable_list rand_vars, data_list data) : m_rand_vars(rand_vars),
                                                          m_data(data) {
    if (std::is_sorted(rand_vars.begin(),
                       rand_vars.end(),
                       [](auto a, auto b) -> bool { return a.id() <= b.id(); })) {
// TODO check precondition on cardinality of data
// throw exception if length(data) != product(lengths(rand_vars))
      m_data = data;
      m_rand_vars
    }

  }
    rand_vars(rand_vars), data(data) {
    auto perms = find_permutation_indices(rand_vars);

  }

}  // namespace pgmfactors
