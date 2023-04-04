// pgmexample.cpp
#include <iostream>
#include <map>
#include <string_view>
#include <xtensor/xio.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xaxis_iterator.hpp>
#include <xtensor/xaxis_slice_iterator.hpp>
#include "pgm/factor.hpp"

template<typename T>
void print_xarray_with_coords(const xt::xarray<T> &xa) {
    size_t flat_index{0};
    for(auto val : xa) {
        auto idx_vector = xt::unravel_index(flat_index, xa.shape());
        for(auto idx_coord : idx_vector) {
            std::cout << idx_coord << " ";
        }
        std::cout << val << "\n";
        flat_index++;
    }
    std::cout << std::flush;
}

void print_factor(const pgm::factor& f, std::string_view name="") {
    auto vars = f.vars();
    std::cout << "Factor " << name << ":\n";
    for(auto v : vars) {
        std::cout << v.id() << " ";
    }
    std::cout << "\n-----\n";
    print_xarray_with_coords(f.data());
    std::cout << "\n";
    std::cout << std::flush;
}


int main() {

    auto rv_A = pgm::rv{3};
    auto rv_B = pgm::rv{2};
    auto rv_C = pgm::rv{2};
    pgm::factor f_AB(pgm::factor::rv_list{rv_A, rv_B}, {{0.5, 0.8, 0.1, 0.0, 0.3, 0.9}});
    print_factor(f_AB, "AB");
    pgm::factor f_AB2(pgm::factor::rv_list{rv_B, rv_A}, {{0.5, 0.1, 0.3, 0.8, 0.0, 0.9}});
    print_factor(f_AB2, "AB2");

    pgm::factor f_BC(pgm::factor::rv_list{rv_B, rv_C}, {{0.5, 0.7, 0.1, 0.2}});

    auto f_ABC = pgm::factor_product(f_AB, f_BC);

    print_factor(f_AB, "AB");
    print_factor(f_BC, "BC");
    print_factor(f_ABC, "ABC");


    print_factor(factor_marginalization(f_ABC, rv_A), "ABC-1");
    print_factor(factor_marginalization(f_ABC, rv_B), "ABC-2");
    print_factor(factor_marginalization(f_ABC, rv_C), "ABC-3");


    pgm::factor f_d1(pgm::factor::rv_list{rv_A, rv_B}, {{0.5, 0.2, 0.0000001, 0.0000001, 0.3, 0.45}});
    pgm::factor f_d2(pgm::factor::rv_list{rv_A}, {{0.8, 0.001, 0.6}});
    print_factor(f_d1, "d1");
    print_factor(f_d2, "d2");
    print_factor(factor_division(f_d1, f_d2), "d1/d2");

    return 0;
}
