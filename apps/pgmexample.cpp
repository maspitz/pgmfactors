// pgmexample.cpp
#include <iostream>
#include <xtensor/xstrides.hpp>
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
        std::cout << val << std::endl;
        flat_index++;
    }
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

    print_xarray_with_coords(f.data());
    return 0;
}
