// pgmfactors_test.cpp

#include "pgmfactors/pgmfactors.hpp"
#include <catch2/catch_test_macros.hpp>


TEST_CASE("Factor Equality Operation") {
    auto AB_data = xt::xarray<double>{0.5, 0.8, 0.1, 0.0, 0.3, 0.9};
    AB_data.reshape({3, 2});
    auto AB_vars = pgmfactors::factor::rv_list{1,2};
    pgmfactors::factor f_AB(AB_vars,
                            pgmfactors::factor::data_array(AB_data));

    REQUIRE(f_AB == f_AB);

    auto AB_data_tweaked = xt::xarray<double>{1.0, 0.8, 0.1, 0.0, 0.3, 0.9};
    AB_data_tweaked.reshape({3, 2});
    pgmfactors::factor f_AB_data_tweaked(AB_vars,
                            pgmfactors::factor::data_array(AB_data_tweaked));

    REQUIRE(f_AB != f_AB_data_tweaked);

    auto AB_vars_tweaked = pgmfactors::factor::rv_list{1,3};
    pgmfactors::factor f_AB_vars_tweaked(AB_vars_tweaked,
                            pgmfactors::factor::data_array(AB_data));

    REQUIRE(f_AB != f_AB_vars_tweaked);

}

TEST_CASE("Factor Index Sorting") {
    auto AB_data = xt::xarray<double>{0.5, 0.8, 0.1, 0.0, 0.3, 0.9};
    AB_data.reshape({3, 2});
    auto AB_vars = pgmfactors::factor::rv_list{1,2};
    pgmfactors::factor f_AB(AB_vars,
                            pgmfactors::factor::data_array(AB_data));

    // Now construct an equivalent factor, but with variables given
    // in reversed order: {2,1} rather than {1,2}.
    auto AB_permuted_data = xt::xarray<double>{0.5, 0.1, 0.3, 0.8, 0.0, 0.9};
    AB_permuted_data.reshape({2, 3});
    auto AB_permuted_vars = pgmfactors::factor::rv_list{2,1};
    pgmfactors::factor f_AB_permuted(AB_permuted_vars,
                                     pgmfactors::factor::data_array(AB_permuted_data));

    REQUIRE(f_AB == f_AB_permuted);
}

TEST_CASE("Factor Product", "[factor][operation]")
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

    SECTION("Product of factors sharing one variable") {
    // Perform factor product and compare to expected result.
    // Use is_close() because floating point calculations
    // aren't expected to yield strict equality.

      auto calculated_product = pgmfactors::factor_product(f_AB, f_BC);

      auto expected_vars = pgmfactors::factor::rv_list{1,2,3};
      auto expected_data = xt::xarray<double>{
      0.25, 0.35, 0.08, 0.16,
      0.05, 0.07, 0.00, 0.00,
      0.15, 0.21, 0.09, 0.18};
      expected_data.reshape({3, 2, 2});
      auto expected_product = pgmfactors::factor(expected_vars,
                                                 pgmfactors::factor::data_array(expected_data));
      REQUIRE(is_close(calculated_product, expected_product));
    }

    SECTION("Product of factor with itself") {
      auto calculated_self_product = pgmfactors::factor_product(f_AB, f_AB);
      auto expected_vars = f_AB.vars();
      auto expected_data = f_AB.data() * f_AB.data();  // elementwise multiplication
      auto expected_self_product = pgmfactors::factor(expected_vars,
                                                      pgmfactors::factor::data_array(expected_data));
      REQUIRE(is_close(calculated_self_product, expected_self_product));
    }


}
