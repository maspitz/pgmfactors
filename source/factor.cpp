// factor.cpp
#include <string>
#include <vector>
#include <algorithm>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xstrided_view.hpp>

#include "pgm/factor.hpp"

namespace pgm
{

  auto is_close(const factor& f_a, const factor& f_b, double rtol, double atol) -> bool {
    return (f_a.vars() == f_b.vars()) &&
      (xt::allclose(f_a.data(), f_b.data(), rtol, atol));
  }


auto factor_reduction(const factor& f_a, const std::map<int,int>& assignment) -> factor {
  auto a_vars = f_a.vars();
  factor::rv_list b_vars;

  xt::xstrided_slice_vector sv;

  auto it_avars = a_vars.begin();
  auto it_x = assignment.begin();

  while(it_avars != a_vars.end() && it_x != assignment.end()) {
    if (*it_avars > it_x->first) {
      // assignment rv not in input factor: advance to next assignment rv.
      it_x++;
    } else if (*it_avars < it_x->first) {
      // input factor rv is unassigned: include axis in output factor.
      sv.push_back(xt::all());
      b_vars.push_back(*it_avars);
      it_avars++;
    } else {
      // input factor rv is assigned: reduce output factor according to the assignment.
      sv.push_back(it_x->second);
      it_avars++;
      it_x++;
    }
  }
  if (it_avars != a_vars.end()) {
    // add any remaining unassigned random variables and include their axes.
    sv.insert(sv.end(), (a_vars.end() - it_avars), xt::all());
    b_vars.insert(b_vars.end(), it_avars, a_vars.end());
  }
  return pgm::factor(b_vars, xt::strided_view(f_a.data(), sv));
}


auto factor_reduction2(const factor& f_a, const std::map<int,int>& assignment) -> factor {
  auto a_vars = f_a.vars();
  factor::rv_list b_vars;
  xt::xstrided_slice_vector sv;

  for(auto avar : a_vars) {
    if (auto value = assignment.find(avar); value != assignment.end()) {
      sv.push_back(value->second);
    }
    else
    {
      b_vars.push_back(avar);
      sv.push_back(xt::all());
    }
  }

  return pgm::factor(b_vars, xt::strided_view(f_a.data(), sv));
}



auto factor_marginalization(const factor& f_a, int summation_rv) -> factor {
  auto a_vars = f_a.vars();
  auto a_rv_it = std::find(a_vars.begin(), a_vars.end(), summation_rv);
  if (a_rv_it == a_vars.end()) {
    return f_a;
  }
  int rv_axis = a_rv_it - a_vars.begin();
  factor::rv_list b_vars;
  b_vars.insert(b_vars.end(),a_vars.begin(), a_rv_it);
  b_vars.insert(b_vars.end(),a_rv_it + 1, a_vars.end());

  return pgm::factor(b_vars, xt::sum(f_a.data(), {rv_axis}));
}

auto factor_product(const factor& f_a, const factor& f_b) -> factor {
	auto a_vars = f_a.vars();
	auto a_card = f_a.data().shape();

	auto b_vars = f_b.vars();
	auto b_card = f_b.data().shape();

	auto product_vars = factor::rv_list{};
	auto product_card = std::vector<int>{};
	auto a_reshape = std::vector<int>{};
	auto b_reshape = std::vector<int>{};

	// Merge a_vars and b_vars [invariant: strictly ascending order]
	auto it_avars = a_vars.begin(),	it_bvars = b_vars.begin();
	auto it_acard = a_card.begin(), it_bcard = b_card.begin();

	// Step through both a_vars and b_vars,
	// appending the lesser-id variable at each step.
	while(it_avars != a_vars.end() && it_bvars != b_vars.end()) {
		if (*it_avars < *it_bvars) {
			product_vars.push_back(*it_avars);
			product_card.push_back(*it_acard);
			a_reshape.push_back(*it_acard);
			b_reshape.push_back(1);
			++it_avars; ++it_acard;
		} else if (*it_bvars < *it_avars) {
			product_vars.push_back(*it_bvars);
			product_card.push_back(*it_bcard);
			a_reshape.push_back(1);
			b_reshape.push_back(*it_bcard);
			++it_bvars; ++it_bcard;
		} else {
			// ASSERT(*it_acard == *it_bcard);
			product_vars.push_back(*it_avars);
			product_card.push_back(*it_acard);
			a_reshape.push_back(*it_acard);
			b_reshape.push_back(*it_acard);
			++it_avars; ++it_acard;
			++it_bvars; ++it_bcard;
		}
	}
	if (it_avars != a_vars.end()) {
	// Add any remaining elements from a_vars
		product_vars.insert(product_vars.end(), it_avars, a_vars.end());
		product_card.insert(product_card.end(), it_acard, a_card.end());
		a_reshape.insert(a_reshape.end(), it_acard, a_card.end());
		b_reshape.insert(b_reshape.end(), (a_card.end() - it_acard), 1);
	} else if (it_bvars != b_vars.end()) {
	// Add any remaining elements from b_vars
		product_vars.insert(product_vars.end(), it_bvars, b_vars.end());
		product_card.insert(product_card.end(), it_bcard, b_card.end());
		b_reshape.insert(b_reshape.end(), it_bcard, b_card.end());
		a_reshape.insert(a_reshape.end(), (b_card.end() - it_bcard), 1);
	}

	// make factor data views using a_reshape and b_reshape
	// multiply the views to obtain the product_data
	// return the product factor constructed from the
	// product_vars and product_data

  auto view_a = xt::reshape_view(f_a.data(), a_reshape);
  auto view_b = xt::reshape_view(f_b.data(), b_reshape);
  return factor(product_vars, view_a * view_b);
}


auto factor_division(const factor& f_a, const factor& f_b) -> factor {
	auto a_vars = f_a.vars();
	auto a_card = f_a.data().shape();

	auto b_vars = f_b.vars();
	auto b_card = f_b.data().shape();

	auto b_reshape = std::vector<int>{};

	// Merge a_vars and b_vars [invariant: strictly ascending order]
	auto it_avars = a_vars.begin(),	it_bvars = b_vars.begin();
	auto it_bcard = b_card.begin();

	// Step through both a_vars and b_vars,
	// appending the lesser-id variable at each step.
	while(it_avars != a_vars.end()) {
		if (*it_avars < *it_bvars) {
			b_reshape.push_back(1);
			++it_avars;
		} else if (*it_bvars < *it_avars) {
      throw std::runtime_error("Scope mismatch: some variable in f_b is not present in f_a");
		} else {
      b_reshape.push_back(*it_bcard);
      ++it_avars;
      ++it_bvars;
      ++it_bcard;
		}
	}
  if (it_bvars != b_vars.end()) {
      throw std::runtime_error("Scope mismatch: some variable in f_b is not present in f_a");
  }

  auto view_b = xt::reshape_view(f_b.data(), b_reshape);
  // TODO: handle the case of 0 / 0, in which we are to take the relaxed view
  // that 0 / 0 should be taken to yield 0.
  return factor(a_vars, f_a.data() / view_b);
}



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

}  // namespace pgm