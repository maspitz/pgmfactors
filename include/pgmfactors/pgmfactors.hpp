#pragma once

#include <string>
#include <vector>
#include <xtensor/xarray.hpp>

namespace pgm {

  struct var {

  };

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
      // TODO consider calling the set of variables the 'scope'. (p.104)
      rv_list m_rand_vars;
      xt::xarray<value_type> m_data;
    public:
      explicit factor(const rv_list &rand_vars, const data_array &data);
      auto data() const -> const data_array& { return m_data; }
      auto vars() const -> const rv_list& { return m_rand_vars; }
      auto operator==(const factor&) const -> bool = default;
  };


  class factor;

  // Tests if two factors are identical in scope and numerically close in value.
  // rtol: relative tolerance
  // atol: absolute tolerance
  auto is_close(const factor& f_a, const factor& f_b,
                double rtol = 1e-5, double atol = 1e-8) -> bool;


  auto factor_product(const factor& f_a, const factor& f_b) -> factor;
  auto factor_reduction(const factor& f_a, const std::map<int,int>& assignment) -> factor;
  auto factor_reduction2(const factor& f_a, const std::map<int,int>& assignment) -> factor;
  auto factor_marginalization(const factor& f_a, int rv_idx) -> factor;
  auto factor_division(const factor& f_a, const factor& f_b) -> factor;


} // namespace pgm

