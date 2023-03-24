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

    auto rvlst = pgmfactors::factor::rv_list{1,2};
    auto mydata = pgmfactors::factor::data_array({{1,2},{3,4}});
    pgmfactors::factor f(rvlst, mydata);

    print_xarray_with_coords(f.data());
    return 0;
}
