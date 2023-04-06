// test_factor.cpp

#include <stdexcept>
#include <catch2/catch_test_macros.hpp>

#include "pgm/factor.hpp"

TEST_CASE("Factor Construction") {
  pgm::rv rv_A(3), rv_B(2), rv_C(2);
  SECTION("The size of a factor's data corresponds to cardinality of its scope variables")
  {
    CHECK_THROWS_AS(pgm::factor(pgm::factor::rv_list {rv_B, rv_C},
                                {{0.5, 0.8, 0.1, 0.0, 0.3, 0.9}}),
                    std::runtime_error);
    CHECK_NOTHROW(pgm::factor(pgm::factor::rv_list {rv_B, rv_C},
                              {{0.5, 0.8, 0.1, 0.0}}));
  }
  SECTION("A factor's scope does not permit repeated random variables")
  {
    CHECK_THROWS(pgm::factor(pgm::factor::rv_list {rv_B, rv_B},
                              {{0.5, 0.8, 0.1, 0.0}}));
  }

}


TEST_CASE("Factor Equality Operation")
{
  auto rv_A = pgm::rv {3};
  auto rv_B = pgm::rv {2};
  auto rv_C = pgm::rv {2};

  auto AB_data = xt::xarray<double> {0.5, 0.8, 0.1, 0.0, 0.3, 0.9};
  AB_data.reshape({3, 2});
  auto AB_vars = pgm::factor::rv_list {rv_A, rv_B};
  pgm::factor f_AB(AB_vars, AB_data);
  //                            pgm::factor::data_array(AB_data));

  REQUIRE(f_AB == f_AB);

  auto AB_data_tweaked = xt::xarray<double> {1.0, 0.8, 0.1, 0.0, 0.3, 0.9};
  AB_data_tweaked.reshape({3, 2});
  pgm::factor f_AB_data_tweaked(AB_vars,
                                pgm::factor::data_array(AB_data_tweaked));

  REQUIRE(f_AB != f_AB_data_tweaked);

  auto AB_vars_tweaked = pgm::factor::rv_list {rv_A, rv_C};
  pgm::factor f_AB_vars_tweaked(AB_vars_tweaked,
                                pgm::factor::data_array(AB_data));

  REQUIRE(f_AB != f_AB_vars_tweaked);
}


TEST_CASE("Factor Index Sorting")
{
  auto rv_A = pgm::rv {3};
  auto rv_B = pgm::rv {2};
  auto AB_data = xt::xarray<double> {0.5, 0.8, 0.1, 0.0, 0.3, 0.9};
  AB_data.reshape({3, 2});
  auto AB_vars = pgm::factor::rv_list {rv_A, rv_B};
  pgm::factor f_AB(AB_vars, pgm::factor::data_array(AB_data));

  // Now construct an equivalent factor, but with variables given
  // in reversed order: {2,1} rather than {1,2}.
  auto AB_permuted_data = xt::xarray<double> {0.5, 0.1, 0.3, 0.8, 0.0, 0.9};
  AB_permuted_data.reshape({2, 3});
  auto AB_permuted_vars = pgm::factor::rv_list {rv_B, rv_A};
  pgm::factor f_AB_permuted(AB_permuted_vars,
                            pgm::factor::data_array(AB_permuted_data));

  REQUIRE(f_AB == f_AB_permuted);
}

TEST_CASE("Factor Product", "[factor][operation]")
{
  auto rv_A = pgm::rv {3};
  auto rv_B = pgm::rv {2};
  auto rv_C = pgm::rv {2};

  auto AB_data = xt::xarray<double> {0.5, 0.8, 0.1, 0.0, 0.3, 0.9};
  AB_data.reshape({3, 2});
  auto AB_vars = pgm::factor::rv_list {rv_A, rv_B};
  pgm::factor f_AB(AB_vars, pgm::factor::data_array(AB_data));
  ;

  auto BC_data = xt::xarray<double> {0.5, 0.7, 0.1, 0.2};
  BC_data.reshape({2, 2});
  auto BC_vars = pgm::factor::rv_list {rv_B, rv_C};
  pgm::factor f_BC(BC_vars, pgm::factor::data_array(BC_data));
  ;

  SECTION("Product of factors sharing one variable")
  {
    // Perform factor product and compare to expected result.
    // Use is_close() because floating point calculations
    // aren't expected to yield strict equality.

    auto calculated_product = pgm::factor_product(f_AB, f_BC);

    auto expected_vars = pgm::factor::rv_list {rv_A, rv_B, rv_C};
    auto expected_data = xt::xarray<double> {
        0.25, 0.35, 0.08, 0.16, 0.05, 0.07, 0.00, 0.00, 0.15, 0.21, 0.09, 0.18};
    expected_data.reshape({3, 2, 2});
    auto expected_product =
        pgm::factor(expected_vars, pgm::factor::data_array(expected_data));
    REQUIRE(is_close(calculated_product, expected_product));
  }

  SECTION("Product of factor with itself")
  {
    auto calculated_self_product = pgm::factor_product(f_AB, f_AB);
    auto expected_vars = f_AB.vars();
    auto expected_data =
        f_AB.data() * f_AB.data();  // elementwise multiplication
    auto expected_self_product =
        pgm::factor(expected_vars, pgm::factor::data_array(expected_data));
    REQUIRE(is_close(calculated_self_product, expected_self_product));
  }
}
