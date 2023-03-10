#pragma once

#include <string>
#include <vector>

#include "pgmfactors/pgmfactors_export.hpp"

/**
 * A note about the MSVC warning C4251:
 * This warning should be suppressed for private data members of the project's
 * exported classes, because there are too many ways to work around it and all
 * involve some kind of trade-off (increased code complexity requiring more
 * developer time, writing boilerplate code, longer compile times), but those
 * solutions are very situational and solve things in slightly different ways,
 * depending on the requirements of the project.
 * That is to say, there is no general solution.
 *
 * What can be done instead is understand where issues could arise where this
 * warning is spotting a legitimate bug. I will give the general description of
 * this warning's cause and break it down to make it trivial to understand.
 *
 * C4251 is emitted when an exported class has a non-static data member of a
 * non-exported class type.
 *
 * The exported class in our case is the class below (exported_class), which
 * has a non-static data member (m_name) of a non-exported class type
 * (std::string).
 *
 * The rationale here is that the user of the exported class could attempt to
 * access (directly, or via an inline member function) a static data member or
 * a non-inline member function of the data member, resulting in a linker
 * error.
 * Inline member function above means member functions that are defined (not
 * declared) in the class definition.
 *
 * Since this exported class never makes these non-exported types available to
 * the user, we can safely ignore this warning. It's fine if there are
 * non-exported class types as private member variables, because they are only
 * accessed by the members of the exported class itself.
 *
 * The name() method below returns a pointer to the stored null-terminated
 * string as a fundamental type (const char), so this is safe to use anywhere.
 * The only downside is that you can have dangling pointers if the pointer
 * outlives the class instance which stored the string.
 *
 * Shared libraries are not easy, they need some discipline to get right, but
 * they also solve some other problems that make them worth the time invested.
 */

/**
 * @brief Reports the name of the library
 *
 * Please see the note above for considerations when creating shared libraries.
 */
class PGMFACTORS_EXPORT exported_class
{
public:
  /**
   * @brief Initializes the name field to the name of the project
   */
  exported_class();

  /**
   * @brief Returns a non-owning pointer to the string stored in this class
   */
  auto name() const -> const char*;

private:
  PGMFACTORS_SUPPRESS_C4251
  std::string m_name;
};


namespace pgmfactors {

// random_variable identifies a finite valued random variable.
// It doesn't model probability.  That is handled by Factor.
class random_variable {
    const int m_cardinality;
    const int m_id;
    const std::string m_name;
    inline static int m_next_id = 0;
  public:
    explicit random_variable(int cardinality) :
      m_cardinality(cardinality), m_id(++m_next_id), m_name("unnamed rv") {}
    explicit random_variable(std::string name, int cardinality) :
      m_cardinality(cardinality), m_id(++m_next_id), m_name(std::move(name)) {}
    [[nodiscard]] auto id() const noexcept { return m_id; }
    [[nodiscard]] auto cardinality() const noexcept { return m_cardinality; }
    [[nodiscard]] auto name() const noexcept { return m_name; }
};

  // Note: Factor models only discrete factors at this time
  // Class invariants
  //   rv_id is in strictly ascending order
  //   elements of rv_id are not repeated
  class factor {
    public:
      using value_type = double;
      using variable_list = const std::vector<random_variable>;
      using data_list = const std::vector<value_type>;
    private:
      variable_list m_rand_vars;
      xt::xarray<value_type> m_data;
    public:
/*      explicit factor(variable_list rand_vars,
                      data_list data) : rand_vars(rand_vars), data(data) { } */
      explicit factor(variable_list rand_vars, data_list data);
}
