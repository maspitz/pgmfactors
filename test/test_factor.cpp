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

  pgm::factor f_AB(pgm::factor::rv_list {rv_A, rv_B},
                   {{0.5, 0.8, 0.1, 0.0, 0.3, 0.9}});
  pgm::factor f_AB_other_data(pgm::factor::rv_list {rv_A, rv_B},
                   {{0.1, 0.8, 0.1, 0.0, 0.3, 0.9}});
  pgm::factor f_AB_other_scope(pgm::factor::rv_list {rv_A, rv_C},
                   {{0.5, 0.8, 0.1, 0.0, 0.3, 0.9}});

  REQUIRE(f_AB == f_AB);
  REQUIRE(f_AB != f_AB_other_data);
  REQUIRE(f_AB != f_AB_other_scope);
}


TEST_CASE("Factor Index Sorting")
{
  auto rv_A = pgm::rv {3};
  auto rv_B = pgm::rv {2};
  pgm::factor f_AB(pgm::factor::rv_list {rv_A, rv_B},
                   {{0.5, 0.8, 0.1, 0.0, 0.3, 0.9}});

  // Now construct an equivalent factor, but with variables given
  // in reversed order: {2,1} rather than {1,2}.
  pgm::factor f_AB_swap_vars(pgm::factor::rv_list {rv_B, rv_A},
                             {{0.5, 0.1, 0.3, 0.8, 0.0, 0.9}});

  REQUIRE(f_AB == f_AB_swap_vars);
}

TEST_CASE("Factor Product", "[factor][operation]")
{
  auto rv_A = pgm::rv {3};
  auto rv_B = pgm::rv {2};
  auto rv_C = pgm::rv {2};

  pgm::factor f_AB(pgm::factor::rv_list {rv_A, rv_B},
                   {{0.5, 0.8, 0.1, 0.0, 0.3, 0.9}});
  pgm::factor f_BC(pgm::factor::rv_list {rv_B, rv_C},
                   {{0.5, 0.7, 0.1, 0.2}});

  SECTION("Product of factors sharing one variable")
  {
    // Perform factor product and compare to expected result.
    // Use is_close() because floating point calculations
    // aren't expected to yield strict equality.

    auto calculated_product = pgm::factor_product(f_AB, f_BC);

    pgm::factor expected_product(pgm::factor::rv_list {rv_A, rv_B, rv_C},
                                 {{0.25, 0.35, 0.08, 0.16, 0.05, 0.07,
                                    0.00, 0.00, 0.15, 0.21, 0.09, 0.18}});
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

TEST_CASE("Factor Reduction", "[factor][operation]")
{
  auto rv_A = pgm::rv {3};
  auto rv_B = pgm::rv {2};
  auto rv_C = pgm::rv {2};

  pgm::factor f_AB(pgm::factor::rv_list {rv_A, rv_B},
                   {{0.5, 0.8, 0.1, 0.0, 0.3, 0.9}});
  pgm::factor f_BC(pgm::factor::rv_list {rv_B, rv_C},
                   {{0.5, 0.7, 0.1, 0.2}});

  auto f_ABC = pgm::factor_product(f_AB, f_BC);

  SECTION("Reduction by evidence C = 0")
  {
    pgm::rv_evidence C1_evidence;
    C1_evidence[rv_C] = 0;

    auto calculated_reduction = pgm::factor_reduction(f_ABC, C1_evidence);

    pgm::factor expected_reduction(pgm::factor::rv_list {rv_A, rv_B},
                                 {{0.25, 0.08, 0.05, 0.00, 0.15, 0.09}});
    REQUIRE(is_close(calculated_reduction, expected_reduction));
  }
}
