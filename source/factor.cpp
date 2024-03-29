// factor.cpp
#include <algorithm>
#include <string>
#include <vector>
#include <ranges>

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

auto factor_reduction(const factor& input, const pgm::rv_evidence& assignments)
    -> factor
{
  auto input_vars = input.vars();
  factor::rv_list output_vars;
  xt::xstrided_slice_vector stride;

  auto it_assignment = assignments.begin();
  auto it_input_var = input_vars.begin();

  while (it_input_var != input_vars.end() &&
         it_assignment != assignments.end()) {
    auto assignment_id = it_assignment->first.id();
    auto in_var_id = it_input_var->id();
    if (in_var_id < assignment_id) {
      output_vars.push_back(*it_input_var);
      stride.push_back(xt::all());
      ++it_input_var;
    } else if (in_var_id > assignment_id) {
      ++it_assignment;
    } else {
      stride.push_back(it_assignment->second);
      ++it_input_var;
      ++it_assignment;
    }
  }
  while (it_input_var != input_vars.end()) {
    output_vars.push_back(*it_input_var);
    stride.push_back(xt::all());
    ++it_input_var;
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
  auto rv_axis = a_rv_it - input_vars.begin();
  factor::rv_list output_vars;
  output_vars.reserve(input_vars.size() - 1);
  output_vars.insert(output_vars.end(), input_vars.begin(), a_rv_it);
  output_vars.insert(output_vars.end(), a_rv_it + 1, input_vars.end());

  return pgm::factor(output_vars, xt::sum(input.data(), {rv_axis}));
}


// venn_action() takes two sorted ranges and acts with the given functions on their elements.
//
// The idea is to:
// call Function1 for elements of the set difference (Set1 - Set2),
// call Function2 for elements of the set difference (Set2 - Set1), and
// call Function12 for elements of the set intersection (Set1 and Set2).
//
// venn_action() returns nothing.  Its effectiveness depends entirely on the side effects
// of the given functions.
template <class InputIterator1, class InputIterator2,
          class Function1, class Function2, class Function12, class Compare>
  void venn_action (InputIterator1 first1, InputIterator1 last1,
                    InputIterator2 first2, InputIterator2 last2,
                    Function1 only_1,
                    Function2 only_2,
                    Function12 both_1_and_2,
                    Compare comp) {
  while(first1 != last1 && first2 != last2) {
    if (comp(*first1, *first2)) {
      only_1(first1);
      ++first1;
    } else if (comp(*first2, *first1)) {
      only_2(first2);
      ++first2;
    } else {
      both_1_and_2(first1);
      ++first1;
      ++first2;
    }
  }
  while(first1 != last1) {
    only_1(first1);
    ++first1;
  }
  while(first2 != last2) {
    only_2(first2);
    ++first2;
  }
}

auto factor_product(const factor& f_a, const factor& f_b) -> factor
{
  auto a_vars = f_a.vars();
  auto b_vars = f_b.vars();

  std::vector<pgm::rv> product_vars;
  std::vector<int> a_shape, b_shape;

  venn_action(a_vars.begin(), a_vars.end(),
              b_vars.begin(), b_vars.end(),
              [&](auto v) { product_vars.push_back(*v); a_shape.push_back(v->card()); b_shape.push_back(1); },
              [&](auto v) { product_vars.push_back(*v); a_shape.push_back(1); b_shape.push_back(v->card()); },
              [&](auto v) { product_vars.push_back(*v); a_shape.push_back(v->card()); b_shape.push_back(v->card()); },
              pgm::rv_id_comparison());

  auto view_a = xt::reshape_view(f_a.data(), a_shape);
  auto view_b = xt::reshape_view(f_b.data(), b_shape);
  return factor(product_vars, view_a * view_b);
}

auto factor_marginalization(const factor& input, const std::vector<pgm::rv>& summation_rvs) -> factor
{
  auto input_vars = input.vars();
  std::vector<int> summation_axes;
  factor::rv_list output_vars;

  venn_action(input_vars.begin(), input_vars.end(),
              summation_rvs.begin(), summation_rvs.end(),
              [&](auto v) { output_vars.push_back(*v); },
              [&](auto v) { },
              [&](auto v) { summation_axes.push_back(v - input_vars.begin()); },
              pgm::rv_id_comparison());

  return pgm::factor(output_vars, xt::sum(input.data(), summation_axes));
}


// factor_reduction2 is an alternative, unexported implementation of factor_reduction.
// It should behave identically, but it uses venn_action internally.
// Its main drawback is that it uses std::map<>::at() because venn_action only iterates
// over keys of the rv_evidence map, not key-value pairs.  Still, it may be useful in the future.
auto factor_reduction2(const factor& input, const pgm::rv_evidence& assignments)
    -> factor
{
  auto assignment_vars = std::views::keys(assignments);
  factor::rv_list output_vars;
  xt::xstrided_slice_vector stride;

  venn_action(input.vars().begin(), input.vars().end(),
              assignment_vars.begin(), assignment_vars.end(),
              [&](auto v) { output_vars.push_back(*v); stride.push_back(xt::all()); },
              [&](auto v) { return; },
              [&](auto v) { stride.push_back(assignments.at(*v)); },
              pgm::rv_id_comparison());

  return pgm::factor(output_vars, xt::strided_view(input.data(), stride));
}

auto factor_normalization(const factor &f, factor::value_type norm) -> factor {
  return pgm::factor(f.vars(), f.data() * norm / xt::sum(f.data()));
}


auto factor_division(const factor& f_a, const factor& f_b) -> factor
{
  auto a_vars = f_a.vars();
  auto b_vars = f_b.vars();

  auto b_shape = std::vector<int> {};

  // Merge a_vars and b_vars [invariant: strictly ascending order]
  auto it_avars = a_vars.begin();
  auto it_bvars = b_vars.begin();

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
  int dimension_product = 1;
  for (auto r : rand_vars) {
    dimension_product *= r.card();
  }
  if (data.size() != dimension_product) {
    throw std::runtime_error(
        "Cardinality of scope variables does not match size of provided factor "
        "data.");
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
