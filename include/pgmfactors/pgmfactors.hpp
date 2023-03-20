#pragma once

#include <string>
#include <vector>


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
};

} // namespace pgmfactors
