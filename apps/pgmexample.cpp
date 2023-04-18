// pgmexample.cpp
#include <iostream>
#include <iterator>
#include <map>
#include <string_view>
#include <vector>

#include <numeric>  // for std::reduce()
#include <ranges>   // for std::views::keys
#include <algorithm>  // for std::set_difference, std::remove_if

#include <xtensor/xarray.hpp>
#include <xtensor/xaxis_iterator.hpp>
#include <xtensor/xaxis_slice_iterator.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xtensor.hpp>

#include "pgm/factor.hpp"
#include "pgm/rv.hpp"

using namespace std;

template<typename T>
void print_xarray_with_coords(const xt::xarray<T>& xa)
{
  size_t flat_index {0};
  for (auto val : xa) {
    auto idx_vector = xt::unravel_index(flat_index, xa.shape());
    for (auto idx_coord : idx_vector) {
      std::cout << idx_coord << " ";
    }
    std::cout << val << "\n";
    flat_index++;
  }
  std::cout << std::flush;
}

void print_factor(const pgm::factor& f, std::string_view name = "")
{
  auto vars = f.vars();
  std::cout << "Factor " << name << ":\n";
  for (auto v : vars) {
    std::cout << v.id() << " ";
  }
  std::cout << "\n-----\n";
  print_xarray_with_coords(f.data());
  std::cout << "\n";
  std::cout << std::flush;
}


void example_misconception() {
  std::cout << "Study group misconception example:" << std::endl;

  pgm::rv A(2), B(2), C(2), D(2);
  pgm::factor phi1(pgm::factor::rv_list {A, B},
                   {{30, 5, 1, 10}});
  pgm::factor phi2(pgm::factor::rv_list {B, C},
                   {{100, 1, 1, 100}});
  pgm::factor phi3(pgm::factor::rv_list {C, D},
                   {{1, 100, 100, 1}});
  pgm::factor phi4(pgm::factor::rv_list {D, A},
                   {{100, 1, 1, 100}});

  auto phi12 = pgm::factor_product(phi1, phi2);
  auto phi34 = pgm::factor_product(phi3, phi4);
  auto phi1234 = pgm::factor_product(phi12, phi34);

  print_factor(phi1234, "phi1234");
}


void example_factor_operations() {
  pgm::rv rv_A(3), rv_B(2), rv_C(2);
  pgm::factor f_AB(pgm::factor::rv_list {rv_A, rv_B},
                   {{0.5, 0.8, 0.1, 0.0, 0.3, 0.9}});
  print_factor(f_AB, "AB");
  pgm::factor f_AB2(pgm::factor::rv_list {rv_B, rv_A},
                    {{0.5, 0.1, 0.3, 0.8, 0.0, 0.9}});
  print_factor(f_AB2, "AB2");

  pgm::factor f_BC(pgm::factor::rv_list {rv_B, rv_C}, {{0.5, 0.7, 0.1, 0.2}});

  auto f_ABC = pgm::factor_product(f_AB, f_BC);

  print_factor(f_AB, "AB");
  print_factor(f_BC, "BC");
  print_factor(f_ABC, "ABC");

  print_factor(factor_marginalization(f_ABC, rv_A), "ABC-1");
  print_factor(factor_marginalization(f_ABC, rv_B), "ABC-2");
  print_factor(factor_marginalization(f_ABC, rv_C), "ABC-3");

  pgm::factor f_d1(pgm::factor::rv_list {rv_A, rv_B},
                   {{0.5, 0.2, 0.0000001, 0.0000001, 0.3, 0.45}});
  pgm::factor f_d2(pgm::factor::rv_list {rv_A}, {{0.8, 0.001, 0.6}});
  print_factor(f_d1, "d1");
  print_factor(f_d2, "d2");
  print_factor(factor_division(f_d1, f_d2), "d1/d2");

  std::cout << "ABC reduced by evidence C = 0" << std::endl;
  pgm::rv_evidence C1_evidence;
  C1_evidence[rv_C] = 0;
  auto f_ABC_C1 = pgm::factor_reduction(f_ABC, C1_evidence);
  print_factor(f_ABC_C1);
}


// TODO: consider taking std::vector<std::reference_wrapper<pgm::factor>> as the argument.
// or for that matter, begin/end iterators that yield reference wrappers.
pgm::factor factor_joint_product(const std::vector<pgm::factor>& f) {
/*  return std::reduce(f.begin(), f.end(), [](pgm::factor const * a,
    pgm::factor const * b){return pgm::factor_product(*a, *b);});*/
  auto it = f.cbegin();
  // FIXME: this will crash if the input vector is empty
  // TODO: change factor_product to permit an identity-factor as input.
  // Then we can accumulate the product in the ordinary way.
  auto jpd = *it;
  it++;
  while(it != f.cend()) {
    jpd = pgm::factor_product(jpd, *it);
    it++;
  }
  return jpd;
}

pgm::factor naive_cpd_inference(const pgm::factor& jpd,
                                const pgm::rv_evidence& query,
                                const pgm::rv_evidence& evidence) {
  // apply evidence (i.e., project on relevant axes)
  auto f = pgm::factor_reduction(jpd, evidence);

  // marginalize (i.e., sum on all but the query var axes)

  std::vector<pgm::rv> margin_vars;
  auto query_vars = std::views::keys(query);
  std::set_difference(jpd.vars().begin(), jpd.vars().end(),
                      query_vars.begin(), query_vars.end(),
                      std::back_inserter(margin_vars),
                      pgm::rv_id_comparison());
  f = pgm::factor_marginalization(f, margin_vars);

  // normalize
  f = pgm::factor_normalization(f);

  // return query
  return pgm::factor_reduction(f, query);
}


// sum-product variable elimination

std::vector<pgm::factor> sum_product_elimination(const std::vector<pgm::factor>& factors,
                                                 const std::vector<pgm::rv>& elimination_vars) {
  std::vector<pgm::factor> flist{factors};
  for(auto v: elimination_vars) {
    std::vector<pgm::factor> do_product;
    copy_if(std::begin(flist), std::end(flist),
            std::back_inserter(do_product),
            [v](const pgm::factor &f) { return f.scope_contains(v); });
    flist.erase(remove_if(std::begin(flist), std::end(flist),
                            [v](const pgm::factor &f) { return f.scope_contains(v); }),
                  std::end(flist));
    flist.push_back(pgm::factor_marginalization(factor_joint_product(do_product), v));
  }
  return flist;
}

void example_sum_product_ve() {
  pgm::rv D(2), I(2), G(3), S(2), L(2);
  pgm::factor f_D(pgm::factor::rv_list {D}, {{0.6, 0.4}});
  pgm::factor f_I(pgm::factor::rv_list {I}, {{0.7, 0.3}});
  pgm::factor f_G(pgm::factor::rv_list {I, D, G},
                  {{0.3, 0.4, 0.3, 0.05, 0.25, 0.7,
                     0.9, 0.08, 0.02, 0.5, 0.3, 0.2}});
  pgm::factor f_S(pgm::factor::rv_list {I, S},
                  {{0.95, 0.05, 0.2, 0.8}});
  pgm::factor f_L(pgm::factor::rv_list {G, L},
                  {{0.1, 0.9, 0.4, 0.6, 0.99, 0.01}});
  ////auto rv_name = std::map<pgm::rv,std::string,pgm::rv_id_comparison>{{D,"D"}, {I,"I"}, {G,"G"}, {S,"S"}, {L,"L"}};
  // variables to be eliminated in the order given.
  // remaining variables can be queried upon.
  auto elimination_order = std::vector<pgm::rv>{D, I, G, S};
  auto factors = std::vector<pgm::factor>{f_D, f_I, f_G, f_S, f_L};

  auto f = factor_joint_product(sum_product_elimination(factors, elimination_order));
  print_factor(f, "P(L=1) ~ 0.502");

  f = factor_joint_product(sum_product_elimination(factors, {G, I, D, S}));
  print_factor(f, "P(L=1) ~ 0.502");
}


// explicitly construct jpd to perform inference
// (generally not efficient for large dimensional distributions)
void example_jpd_inference() {
  pgm::rv D(2), I(2), G(3), S(2), L(2);
  pgm::factor f_D(pgm::factor::rv_list {D}, {{0.6, 0.4}});
  pgm::factor f_I(pgm::factor::rv_list {I}, {{0.7, 0.3}});
  pgm::factor f_G(pgm::factor::rv_list {I, D, G},
                  {{0.3, 0.4, 0.3, 0.05, 0.25, 0.7,
                     0.9, 0.08, 0.02, 0.5, 0.3, 0.2}});
  pgm::factor f_S(pgm::factor::rv_list {I, S},
                  {{0.95, 0.05, 0.2, 0.8}});
  pgm::factor f_L(pgm::factor::rv_list {G, L},
                  {{0.1, 0.9, 0.4, 0.6, 0.99, 0.01}});

  auto student_jpd = factor_joint_product(std::vector<pgm::factor>{
      f_D, f_I, f_G, f_S, f_L
    });

    print_factor(
      naive_cpd_inference(student_jpd, {{I,1},{D,0},{G,1},{S,1},{L,0}}, {}),
      "P(i1,d0,g2,s1,l0) ~ 0.004608");
    print_factor(
      naive_cpd_inference(student_jpd, {{L,1}}, {}),
      "P(l1) ~ 0.502");
    print_factor(
      naive_cpd_inference(student_jpd, {{L,1}}, {{I,0}}),
      "P(l1 | i0) ~ 0.389");
    print_factor(
      naive_cpd_inference(student_jpd, {{L,1}}, {{I,0},{D,0}}),
      "P(l1 | i0, d0) ~ 0.513");
    print_factor(
      naive_cpd_inference(student_jpd, {{I,1}}, {{G,2}}),
      "P(i1 | g3) ~ 0.079");
    print_factor(
      naive_cpd_inference(student_jpd, {{D,1}}, {}),
      "P(d1) = 0.400");
    print_factor(
      naive_cpd_inference(student_jpd, {{D,1}}, {{G,2}}),
      "P(d1 | g3) ~ 0.629");
    print_factor(
      naive_cpd_inference(student_jpd, {{I,1}}, {{L,0}}),
      "P(i1 | l0) ~ 0.14");
    print_factor(
      naive_cpd_inference(student_jpd, {{I,1}}, {{L,0},{G,2}}),
      "P(i1 | l0, g3) ~ 0.079");
    print_factor(
      naive_cpd_inference(student_jpd, {{I,1}}, {{S,1},{G,2}}),
      "P(i1 | s1, g3) ~ 0.578");
    print_factor(
      naive_cpd_inference(student_jpd, {{D,1}}, {{S,1},{G,2}}),
      "P(d1 | s1, g3) ~ 0.76");
    print_factor(
      naive_cpd_inference(student_jpd, {{I,1}}, {{G,2},{D,1}}),
      "P(i1 | g3, d1) ~ 0.11");
    print_factor(
      naive_cpd_inference(student_jpd, {{I,1}}, {{G,1}}),
      "P(i1 | g2) ~ 0.175");
    print_factor(
      naive_cpd_inference(student_jpd, {{I,1}}, {{G,1},{D,1}}),
      "P(i1 | g2, d1) ~ 0.34");

}


int main()
{
  example_factor_operations();
  example_misconception();
  example_sum_product_ve();
  return 0;
}
