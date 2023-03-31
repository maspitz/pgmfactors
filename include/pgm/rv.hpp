#ifndef RV_H_
#define RV_H_
// rv.hpp

namespace pgm {

// rv models a Random Variable
    class rv {
        private:
            const int m_id;    // Unique id
            const int m_card;  // Cardinality (the number of values the RV can take on)
            inline static int m_next_id = 0;
        public:
            explicit rv(int card) : m_id(m_next_id++), m_card(card) { }
            auto id() const -> const int { return m_id; }
            auto card() const -> const int { return m_card; }
    };

} // namespace pgm


#endif // RV_H_
