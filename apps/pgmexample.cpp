// pgmexample.cpp
#include <iostream>
#include <map>
#include <string_view>
#include <xtensor/xio.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xaxis_iterator.hpp>
#include <xtensor/xaxis_slice_iterator.hpp>
#include "pgmfactors/pgmfactors.hpp"

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

void print_factor(const pgmfactors::factor& f, std::string_view name="") {
    auto vars = f.vars();
    std::cout << "Factor " << name << ":\n";
    for(auto rv_id : vars) {
        std::cout << rv_id << " ";
    }
    std::cout << "\n-----\n";
    print_xarray_with_coords(f.data());
    std::cout << "\n";
    std::cout << std::flush;
}


int main() {

    auto AB_data = xt::xarray<double>{0.5, 0.8, 0.1, 0.0, 0.3, 0.9};
    AB_data.reshape({3, 2});
    auto AB_vars = pgmfactors::factor::rv_list{1,2};
    pgmfactors::factor f_AB(AB_vars,
                            pgmfactors::factor::data_array(AB_data));;

    auto BC_data = xt::xarray<double>{0.5, 0.7, 0.1, 0.2};
    BC_data.reshape({2, 2});
    auto BC_vars = pgmfactors::factor::rv_list{2,3};
    pgmfactors::factor f_BC(BC_vars,
                            pgmfactors::factor::data_array(BC_data));;

    auto rvlst = pgmfactors::factor::rv_list{1,2};
    auto mydata = pgmfactors::factor::data_array({{1,2},{3,4}});
    pgmfactors::factor f(rvlst, mydata);

    auto f_ABC = pgmfactors::factor_product(f_AB, f_BC);
    print_factor(f,"f");

    print_factor(f_AB, "AB");
    print_factor(f_BC, "BC");
    print_factor(f_ABC, "ABC");


    print_factor(pgmfactors::factor_reduction(f_ABC, std::map<int,int>{{1,0}}), "rv1 = 0");
    print_factor(pgmfactors::factor_reduction(f_ABC, std::map<int,int>{{1,1}}), "rv1 = 1");
    print_factor(pgmfactors::factor_reduction(f_ABC, std::map<int,int>{{2,0},{1,1}}), "rv1 = 1, rv2 = 0");
    print_factor(pgmfactors::factor_reduction(f_ABC, std::map<int,int>{{2,1}}), "rv2 = 1");
    print_factor(pgmfactors::factor_reduction(f_ABC, std::map<int,int>{{3,0}}), "rv3 = 0");
    print_factor(pgmfactors::factor_reduction(f_ABC, std::map<int,int>{{3,1}}), "rv3 = 1");


    print_factor(pgmfactors::factor_reduction2(f_ABC, std::map<int,int>{{1,0}}), "rv1 = 0");
    print_factor(pgmfactors::factor_reduction2(f_ABC, std::map<int,int>{{1,1}}), "rv1 = 1");
    print_factor(pgmfactors::factor_reduction2(f_ABC, std::map<int,int>{{2,0},{1,1}}), "rv1 = 1, rv2 = 0");
    print_factor(pgmfactors::factor_reduction2(f_ABC, std::map<int,int>{{2,1}}), "rv2 = 1");
    print_factor(pgmfactors::factor_reduction2(f_ABC, std::map<int,int>{{3,0}}), "rv3 = 0");
    print_factor(pgmfactors::factor_reduction2(f_ABC, std::map<int,int>{{3,1}}), "rv3 = 1");

    return 0;
}
