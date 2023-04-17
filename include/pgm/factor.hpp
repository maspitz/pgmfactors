#ifndef PGM_FACTOR_HPP
#define PGM_FACTOR_HPP
// factor.hpp

#include <string>
#include <vector>
// TODO we could forward-declare xt::xarray<value_type> and use unique_ptr<...>
// to avoid this #include
#include <xtensor/xarray.hpp>

#include "pgm/rv.hpp"

namespace pgm
{

// Note: Factor models only discrete factors at this time
// Class invariants
//   rv_id is in strictly ascending order
//   elements of rv_id are not repeated
class factor
{
public:
  using value_type = double;
  using rv_list = std::vector<pgm::rv>;
  using data_array = xt::xarray<value_type>;

private:
  // TODO consider calling the set of variables the 'scope'. (p.104)
  rv_list m_rand_vars;
  xt::xarray<value_type> m_data;

public:
  explicit factor(const rv_list& rand_vars, const data_array& data);
  auto data() const -> const data_array& { return m_data; }
  auto vars() const -> const rv_list& { return m_rand_vars; }
  auto scope_contains(pgm::rv v) const -> bool {
    return std::find(m_rand_vars.begin(), m_rand_vars.end(), v) != m_rand_vars.end(); }
  auto operator==(const factor&) const -> bool = default;
};

// Tests if two factors are identical in scope and numerically close in value.
// rtol: relative tolerance
// atol: absolute tolerance
auto is_close(const factor& f_a,
              const factor& f_b,
              double rtol = 1e-5,
              double atol = 1e-8) -> bool;

auto factor_product(const factor& f_a, const factor& f_b) -> factor;
auto factor_reduction(const factor& input, const pgm::rv_evidence& assignments)
    -> factor;
auto factor_marginalization(const factor& input, pgm::rv summation_rv) -> factor;
auto factor_marginalization(const factor& input, const std::vector<pgm::rv>& summation_rvs) -> factor;
auto factor_normalization(const factor &f, factor::value_type norm = 1.) -> factor;
auto factor_division(const factor& f_a, const factor& f_b) -> factor;

}  // namespace pgm

#endif  // PGM_FACTOR_HPP
