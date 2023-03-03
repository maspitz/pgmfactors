#include <string>

#include "pgmfactors/pgmfactors.hpp"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("Name is pgmfactors", "[library]")
{
  auto const exported = exported_class {};
  REQUIRE(std::string("pgmfactors") == exported.name());
}
