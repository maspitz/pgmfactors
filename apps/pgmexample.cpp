// pgmexample.cpp
#include <iostream>
#include <map>
#include <string_view>

#include <numeric>  // for std::reduce()

#include <xtensor/xarray.hpp>
#include <xtensor/xaxis_iterator.hpp>
#include <xtensor/xaxis_slice_iterator.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xtensor.hpp>

#include "pgm/factor.hpp"

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
pgm::factor factor_joint_product(const std::vector<pgm::factor const *> f) {
/*  return std::reduce(f.begin(), f.end(), [](pgm::factor const * a,
    pgm::factor const * b){return pgm::factor_product(*a, *b);});*/
  auto it = f.begin();
  auto jpd = **it;
  it++;
  while(it != f.end()) {
    jpd = pgm::factor_product(jpd, **it);
    it++;
  }
  return jpd;
}


// sum-product variable elimination
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
  // first, demonstrate inefficient inference using the full joint product:
  auto student_jpd = factor_joint_product(std::vector<pgm::factor const *>{
      &f_D, &f_I, &f_G, &f_S, &f_L
    });
  pgm::rv_evidence example_specification;
  example_specification[I] = 1;
  example_specification[D] = 0;
  example_specification[G] = 1; // note: in text, Val(G) = {1, 2, 3}.  So 1 represents 2.
  example_specification[S] = 1;
  example_specification[L] = 0;
  auto zz = factor_reduction(student_jpd, example_specification);
  print_factor(zz, "full evidence");

  auto PL = factor_marginalization(student_jpd, I);
  PL = factor_marginalization(PL, D);
  PL = factor_marginalization(PL, G);
  PL = factor_marginalization(PL, S);
  print_factor(PL, "PL");

  pgm::rv_evidence foo;
  foo[I] = 0;
  //auto PLi0 = factor_reduction(student_jpd, pgm::rv_evidence{{I, 0}});
  auto PLi0 = factor_reduction(student_jpd, foo);
  PLi0 = factor_marginalization(PLi0, D);
  PLi0 = factor_marginalization(PLi0, G);
  PLi0 = factor_marginalization(PLi0, S);
//  print_factor(student_jpd);
  print_factor(PLi0, "PL | i0");

  // TODO: marginalize student_jpd to query cpds...
}


int main()
{
  example_factor_operations();
  example_misconception();
  example_sum_product_ve();
  return 0;
}
