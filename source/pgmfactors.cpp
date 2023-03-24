#include <string>
#include <vector>
#include <algorithm>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xadapt.hpp>


#include "pgmfactors/pgmfactors.hpp"

namespace pgmfactors
{
// helper function to permute axes when they are specified out-of-order
// upon construction for a PFactor
auto find_permutation_indices(const factor::rv_list& rvs) -> std::vector<unsigned int>
{
  // Create a vector of pairs that stores the original index of each element
  using index_pair = std::pair<unsigned int, unsigned int>;
  std::vector<index_pair> indexed_input;
  for (unsigned int i = 0; i < rvs.size(); i++) {
    indexed_input.emplace_back(rvs[i], i);
  }

  // Sort the vector of pairs by the element value
  std::sort(indexed_input.begin(),
            indexed_input.end(),
            [](const index_pair& a, const index_pair& b)
            { return a.first < b.first; });

  // Create a vector of permutation indices by mapping the sorted indices back
  // to the original indices
  std::vector<unsigned int> permutation_indices(rvs.size());
  for (unsigned int i = 0; i < rvs.size(); i++) {
    permutation_indices[indexed_input[i].second] = i;
  }

  return permutation_indices;
}


factor::factor(const factor::rv_list &rand_vars, const factor::data_array &data) {

// TODO check precondition on cardinality of data
// throw exception if length(data) != product(lengths(rand_vars))

  if (std::is_sorted(rand_vars.begin(),
                       rand_vars.end(),
                       [](auto a, auto b) -> bool { return a <= b; })) {
    m_data = xt::xarray<value_type>(data);
    m_rand_vars = rand_vars;
  }
  else
  {
    auto perms = find_permutation_indices(rand_vars);
    m_data = xt::transpose(data, perms);
    m_rand_vars = rand_vars;
    std::sort(m_rand_vars.begin(),
              m_rand_vars.end(),
              [](auto a, auto b) -> bool { return a < b; });
    // TODO throw exception if a rand var id is repeated.
  }


}

}  // namespace pgmfactors
