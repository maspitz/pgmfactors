#pragma once

#include <string>
#include <vector>
#include <xtensor/xarray.hpp>

namespace pgmfactors {


  // Note: Factor models only discrete factors at this time
  // Class invariants
  //   rv_id is in strictly ascending order
  //   elements of rv_id are not repeated
  class factor {
    public:
      using value_type = double;
      using rv_list = std::vector<int>;
      using data_array = xt::xarray<value_type>;
    private:
      rv_list m_rand_vars;
      xt::xarray<value_type> m_data;
    public:
/*      explicit factor(rv_list rand_vars,
                      data_array data) : rand_vars(rand_vars), data(data) { } */
      explicit factor(const rv_list rand_vars, const data_array data);
};

} // namespace pgmfactors
