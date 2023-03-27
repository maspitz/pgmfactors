#pragma once

#include <string>
#include <vector>
#include <xtensor/xarray.hpp>

namespace pgmfactors {

  class factor;

  // Tests if two factors are identical in scope and numerically close in value.
  // rtol: relative tolerance
  // atol: absolute tolerance
  auto is_close(const factor& f_a, const factor& f_b,
                double rtol = 1e-5, double atol = 1e-8) -> bool;


  auto factor_product(const factor& f_a, const factor& f_b) -> factor;

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
      explicit factor(const rv_list &rand_vars, const data_array &data);
      auto data() const -> const data_array& { return m_data; }
      auto vars() const -> const rv_list& { return m_rand_vars; }
      auto operator==(const factor&) const -> bool = default;
  };

} // namespace pgmfactors
