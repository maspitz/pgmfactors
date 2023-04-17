#ifndef RV_H_
#define RV_H_
// rv.hpp

#include <map>

namespace pgm
{

// rv models a Random Variable
class rv
{
private:
  int m_id;  // Unique id
  int m_card;  // Cardinality (the number of values the RV can take on)
  inline static int m_next_id = 0;

public:
  explicit rv(int card)
      : m_id(m_next_id++)
      , m_card(card)
  {
  }
  auto id() const -> int { return m_id; }
  auto card() const -> int { return m_card; }
  bool operator==(const rv& other) const { return m_id == other.id(); }
};

struct rv_id_comparison {
    bool operator()(const pgm::rv& lhs, const pgm::rv& rhs) const {
        return lhs.id() < rhs.id();
    }
};

struct rv_id_equality {
    bool operator()(const pgm::rv& lhs, const pgm::rv& rhs) const {
        return lhs.id() == rhs.id();
    }
};

using rv_evidence = std::map<pgm::rv, int, rv_id_comparison>;

}  // namespace pgm

#endif  // RV_H_
