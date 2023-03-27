// pgmfactors_test.cpp

#include "pgmfactors/pgmfactors.hpp"
#include <catch2/catch_test_macros.hpp>


TEST_CASE("Factor Product", "[pgmfactors]")
{


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

    // Perform factor product and compare to expected result.
    // Use is_close() because floating point calculations
    // aren't expected to yield strict equality.

    auto f_product = pgmfactors::factor_product(f_AB, f_BC);

    auto expected_vars = pgmfactors::factor::rv_list{1,2,3};
    auto expected_data = xt::xarray<double>{
    0.25, 0.35, 0.08, 0.16,
    0.05, 0.07, 0.00, 0.00,
    0.15, 0.21, 0.09, 0.18};
    expected_data.reshape({3, 2, 2});
    auto expected_factor = pgmfactors::factor(expected_vars,
                                              pgmfactors::factor::data_array(expected_data));
    REQUIRE(is_close(f_product, expected_factor));

}
