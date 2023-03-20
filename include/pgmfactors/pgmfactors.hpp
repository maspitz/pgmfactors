#pragma once

#include <string>
#include <vector>


namespace pgmfactors {


  // Note: Factor models only discrete factors at this time
  // Class invariants
  //   rv_id is in strictly ascending order
  //   elements of rv_id are not repeated
  class factor {
    public:
      using value_type = double;
      using rv_list = const std::vector<int>;
      using data_array = const xt::xarray<value_type>;
    private:
      variable_list m_rand_vars;
      xt::xarray<value_type> m_data;
    public:
/*      explicit factor(variable_list rand_vars,
                      data_list data) : rand_vars(rand_vars), data(data) { } */
      explicit factor(variable_list rand_vars, data_list data);
};

} // namespace pgmfactors
