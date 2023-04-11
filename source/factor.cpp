// factor.cpp
#include <algorithm>
#include <ranges>
#include <string>
#include <vector>

#include "pgm/factor.hpp"

#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xstrided_view.hpp>

#include "pgm/rv.hpp"

namespace pgm
{

auto is_close(const factor& f_a, const factor& f_b, double rtol, double atol)
    -> bool
{
  return (f_a.vars() == f_b.vars())
      && (xt::allclose(f_a.data(), f_b.data(), rtol, atol));
}

auto factor_reduction(const factor& input,
                      const pgm::rv_evidence& assignments) -> factor
{
  auto input_vars = input.vars();
  factor::rv_list output_vars;

  auto it_assignment = assignments.begin();
  xt::xstrided_slice_vector stride;

  for (auto in_var : input_vars) {
    if (in_var.id() < it_assignment->first.id()) {
      output_vars.push_back(in_var);
      stride.push_back(xt::all());
    } else if (it_assignment->first.id() < in_var.id()) {
      ++it_assignment;
    } else {
      stride.push_back(it_assignment->second);
      ++it_assignment;
    }
  }

  return pgm::factor(output_vars, xt::strided_view(input.data(), stride));
}

auto factor_marginalization(const factor& input, pgm::rv summation_rv) -> factor
{
  auto input_vars = input.vars();
  auto a_rv_it =
      std::find_if(input_vars.begin(),
                   input_vars.end(),
                   [summation_rv](auto var) { return var == summation_rv; });
  if (a_rv_it == input_vars.end()) {
    return input;
  }
  int rv_axis = a_rv_it - input_vars.begin();
  factor::rv_list output_vars;
  output_vars.reserve(input_vars.size() - 1);
  output_vars.insert(output_vars.end(), input_vars.begin(), a_rv_it);
  output_vars.insert(output_vars.end(), a_rv_it + 1, input_vars.end());

  return pgm::factor(output_vars, xt::sum(input.data(), {rv_axis}));
}

auto factor_product(const factor& f_a, const factor& f_b) -> factor
{
  auto a_vars = f_a.vars();
  auto b_vars = f_b.vars();

  std::vector<pgm::rv> product_vars;
  std::vector<int> a_shape, b_shape;

  // Merge a_vars and b_vars [invariant: strictly ascending order]
  auto it_avars = a_vars.begin(), it_bvars = b_vars.begin();

  // Step through both a_vars and b_vars,
  // appending the lesser-id variable at each step.
  while (it_avars != a_vars.end() && it_bvars != b_vars.end()) {
    if (it_avars->id() < it_bvars->id()) {
      product_vars.push_back(*it_avars);
      a_shape.push_back(it_avars->card());
      b_shape.push_back(1);
      ++it_avars;
    } else if (it_bvars->id() < it_avars->id()) {
      product_vars.push_back(*it_bvars);
      a_shape.push_back(1);
      b_shape.push_back(it_bvars->card());
      ++it_bvars;
    } else {
      // ASSERT(*it_acard == *it_bcard);
      product_vars.push_back(*it_avars);
      a_shape.push_back(it_avars->card());
      b_shape.push_back(it_avars->card());
      ++it_avars;
      ++it_bvars;
    }
  }
  while (it_avars != a_vars.end()) {
    // Add any remaining elements from a_vars
    product_vars.push_back(*it_avars);
    a_shape.push_back(it_avars->card());
    b_shape.push_back(1);
    ++it_avars;
  }
  while (it_bvars != b_vars.end()) {
    // Add any remaining elements from b_vars
    product_vars.push_back(*it_bvars);
    a_shape.push_back(1);
    b_shape.push_back(it_bvars->card());
    ++it_bvars;
  }

  // make factor data views using a_reshape and b_reshape
  // multiply the views to obtain the product_data
  // return the product factor constructed from the
  // product_vars and product_data

  auto view_a = xt::reshape_view(f_a.data(), a_shape);
  auto view_b = xt::reshape_view(f_b.data(), b_shape);
  return factor(product_vars, view_a * view_b);
}

auto factor_division(const factor& f_a, const factor& f_b) -> factor
{
  auto a_vars = f_a.vars();
  auto b_vars = f_b.vars();

  auto b_shape = std::vector<int> {};

  // Merge a_vars and b_vars [invariant: strictly ascending order]
  auto it_avars = a_vars.begin(), it_bvars = b_vars.begin();

  // Step through both a_vars and b_vars,
  // appending the lesser-id variable at each step.
  while (it_avars != a_vars.end() && it_bvars != b_vars.end()) {
    if (it_avars->id() < it_bvars->id()) {
      b_shape.push_back(1);
      ++it_avars;
    } else if (it_bvars->id() < it_avars->id()) {
      std::cerr << it_bvars->id() << " < " << it_avars->id() << std::endl;
      throw std::runtime_error(
          "Scope mismatch: some variable in f_b is not present in f_a");
    } else {
      b_shape.push_back(it_bvars->card());
      ++it_avars;
      ++it_bvars;
    }
  }
  while (it_avars != a_vars.end()) {
    b_shape.push_back(1);
    ++it_avars;
  }
  if (it_bvars != b_vars.end()) {
    throw std::runtime_error(
        "Scope mismatch: some variable in f_b is not present in f_a");
  }

  auto view_b = xt::reshape_view(f_b.data(), b_shape);
  // TODO: handle the case of 0 / 0, in which we are to take the relaxed view
  // that 0 / 0 should be taken to yield 0.
  return factor(a_vars, f_a.data() / view_b);
}

// helper function to permute axes when they are specified out-of-order
// upon construction for a PFactor
auto find_permutation_indices(const factor::rv_list& rvs)
    -> std::vector<unsigned int>
{
  // Create a vector of pairs that stores the original index of each element
  using index_pair = std::pair<unsigned int, unsigned int>;
  std::vector<index_pair> indexed_input;
  for (unsigned int i = 0; i < rvs.size(); i++) {
    indexed_input.emplace_back(rvs[i].id(), i);
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

factor::factor(const factor::rv_list& rand_vars, const factor::data_array& data)
{
  // TODO check precondition on cardinality of data
  // throw exception if length(data) != product(lengths(rand_vars))
  int dimension_product = 1;
  for (auto r : rand_vars) {
    dimension_product *= r.card();
  }
  if (data.size() != dimension_product) {
    throw std::runtime_error(
      "Cardinality of scope variables does not match size of provided factor data.");
  }

  m_rand_vars = rand_vars;
  std::vector<int> shape;
  std::transform(rand_vars.begin(),
                 rand_vars.end(),
                 std::back_inserter(shape),
                 [](auto v) { return v.card(); });
  m_data = xt::xarray<value_type>(data);
  m_data.reshape(shape);
  if (!std::is_sorted(
          m_rand_vars.begin(), m_rand_vars.end(), pgm::rv_id_comparison()))
  {
    auto perms = find_permutation_indices(rand_vars);
    m_data = xt::transpose(m_data, perms);
    std::sort(m_rand_vars.begin(), m_rand_vars.end(), pgm::rv_id_comparison());
  }

  if (std::adjacent_find(
          m_rand_vars.begin(), m_rand_vars.end(), pgm::rv_id_equality())
      != m_rand_vars.end())
  {
    throw std::runtime_error(
        "A factor's scope should not have duplicate random variables.");
  }
}

}  // namespace pgm
