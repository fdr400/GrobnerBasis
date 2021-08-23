#pragma once

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <istream>
#include <iterator>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Grobner {
template <typename Field> class Monom {
  public:
    using VariableDegType = uint16_t;
    using VariableNumberType = uint16_t;

    struct Variable {
        VariableDegType deg;
        VariableNumberType number;

        Variable(VariableDegType deg_, VariableNumberType number_)
            : deg(deg_), number(number_) {}

        friend bool operator==(const Variable& v1, const Variable& v2) {
            return (v1.deg == v2.deg) && (v1.number == v2.number);
        }
        // compare Variables by given order of variables.
        friend bool operator<(const Variable& v1, const Variable& v2) {
            if (order.empty()) { // if there is default order of variables
                                 // (x_1 < x_2 < ...).
                return v1.number < v2.number;
            }
            return order[v1.number] < order[v2.number];
        }
    };
    using VariablesListType = std::vector<Variable>;

    // set order of variables.
    static void InitOrder(const std::vector<VariableNumberType>& set_order =
                              {}) { // empty vector implies default order.
        if (!order.empty()) { // clear order vector if it has already been set.
            order.clear();
            order.shrink_to_fit();
        }
        if (set_order.empty()) { // check for default order.
            return;
        }
        order.resize(set_order.size() + 1);
        for (size_t i = 0; i < set_order.size(); ++i) {
            order[set_order[i]] = i;
        }
    }

    Monom() = default;
    Monom(const Field& coef) : coefficient_(coef) {}
    Monom(std::initializer_list<Variable> l, const Field& coef = Field(1))
        : coefficient_(coef) {

        if (coefficient_ != Field(0)) { // check for zero monom.
            variables_.insert(variables_.end(), l.begin(), l.end());
            for (const auto& it : variables_) { // compute monom degree.
                monom_deg_ += it.deg;
            }
            std::sort(variables_.begin(), variables_.end()); // apply set order.
        }
    }
    Monom(std::vector<Variable> l, const Field& coef = Field(1))
        : coefficient_(coef) {

        if (coefficient_ != Field(0)) {
            variables_.insert(variables_.end(), l.begin(), l.end());
            for (const auto& it : variables_) {
                monom_deg_ += it.deg;
            }
            std::sort(variables_.begin(), variables_.end());
        }
    }
    Monom(std::string&& s) {
        std::stringstream st(s);
        st >> *this;
    }
    Monom(const Monom& other) {
        GetDeg() = other.GetDeg();
        GetCoef() = other.GetCoef();
        GetVars() = other.GetVars();
    }
    Monom(Monom&& other) {
        GetDeg() = other.GetDeg();
        GetCoef() = other.GetCoef();
        GetVars().swap(other.GetVars());
    }

    Monom& operator=(const Monom& other) {
        GetDeg() = other.GetDeg();
        GetCoef() = other.GetCoef();
        GetVars() = other.GetVars();

        return *this;
    }
    Monom& operator=(Monom&& other) {
        GetDeg() = other.GetDeg();
        GetCoef() = other.GetCoef();
        GetVars().swap(other.GetVars());

        return *this;
    }

    // make coefficient_ equal to Field(1).
    void Normalize() { coefficient_ = Field(1); }

    VariableDegType& GetDeg() { return monom_deg_; }
    const VariableDegType& GetDeg() const { return monom_deg_; }

    Field& GetCoef() { return coefficient_; }
    const Field& GetCoef() const { return coefficient_; }

    VariablesListType& GetVars() { return variables_; }
    const VariablesListType& GetVars() const { return variables_; }

    friend bool operator==(const Monom& m1, const Monom& m2) {
        return (m1.GetDeg() == m2.GetDeg()) && (m1.GetCoef() == m2.GetCoef()) &&
               (m1.GetVars() == m2.GetVars());
    }

    // multiply monom by the Field member.
    friend Monom operator*(const Monom& m, const Field& f) {
        if (f == Field(0)) { // check if f is equal to zero.
            return Monom();
        }

        Monom result(m);
        result.GetCoef() = result.GetCoef() * f;

        return result;
    }
    friend Monom operator*(const Field& f, const Monom& m) { return m * f; }

    // summarize two monoms.
    // this operator implies that given monoms have the same set of the
    // variables and degrees.
    friend Monom operator+(const Monom& m1, const Monom& m2) {
        // check that given monoms are
        // inverse in addition to each other.
        if (m1.GetCoef() + m2.GetCoef() == 0) {
            return Monom();
        }

        Monom result(m1);
        result.GetCoef() += m2.GetCoef();
        return result;
    }

    // multiply two monoms.
    friend Monom operator*(const Monom& m1, const Monom& m2) {
        Monom result(m1.GetCoef() * m2.GetCoef());

        if (result.GetCoef() == 0) {
            return result;
        }
        result.GetDeg() = m1.GetDeg() + m2.GetDeg();

        const auto& vars1 = m1.GetVars();
        const auto& vars2 = m2.GetVars();
        auto& vars_res = result.GetVars();
        VariableNumberType ptr2 = 0;
        // find for each variable variable with lower_bound number
        // (in set order).
        for (VariableNumberType ptr1 = 0; ptr1 < vars1.size(); ++ptr1) {
            while (ptr2 < vars2.size() && vars2[ptr2] < vars1[ptr1]) {
                vars_res.emplace_back(vars2[ptr2]);
                ptr2 += 1;
            }

            // if numbers of current variables are equal, their degrees sum up.
            if (ptr2 < vars2.size() &&
                vars1[ptr1].number == vars2[ptr2].number) {
                vars_res.emplace_back(vars1[ptr1].deg + vars2[ptr2].deg,
                                      vars1[ptr1].number);
                ptr2 += 1;
            } else {
                vars_res.emplace_back(vars1[ptr1]);
            }
        }
        // add omitted variables in second monom.
        while (ptr2 < vars2.size()) {
            vars_res.emplace_back(vars2[ptr2++]);
        }
        // maintain the order of variables.
        std::sort(vars_res.begin(), vars_res.end());

        return result;
    }

    // check if *this divides by other.
    bool Divides(const Monom& other) const {
        // nothing divides by zero.
        if (other.coefficient_ == Field(0)) {
            return false;
        }

        // if *this has less variables than other, obviously, *this doesn't
        // divide by other.
        if (not_null_vars() < other.not_null_vars()) {
            return false;
        }

        auto ptr = other.GetVars().cbegin();
        auto end = other.GetVars().cend();
        // 1) check if other's list of numbers of variables is the subset of the
        // *this' list of variables.
        // 2) also it is needed for variables with the same numbers in *this and
        // in other, that degree of variable in *this is not less than a degree
        // of a corresponding variable in other.
        for (auto it = GetVars().cbegin(); ptr != end; it = std::next(it)) {
            if (it == GetVars().cend() || *ptr < *it) {
                return false;
            }
            // checks 2).
            if (it->number == ptr->number && it->deg < ptr->deg) {
                return false;
            }
            if (!(*it < *ptr)) {
                ptr = std::next(ptr);
            }
        }

        return true;
    }

    // it implies that m1.Divides(m2) == true.
    // operator/ has the same idea with an operator*, but implemented using
    // iterators.
    friend Monom operator/(const Monom& m1, const Monom& m2) {
        Monom result(m1.GetCoef() / m2.GetCoef());
        if (result.GetCoef() == Field(0)) {
            return result;
        }
        result.GetDeg() = m1.GetDeg() - m2.GetDeg();
        if (result.GetDeg() == 0) {
            return result;
        }

        auto& vars_res = result.GetVars();
        auto ptr = m2.GetVars().cbegin();
        auto end = m2.GetVars().cend();
        for (const auto& var : m1.GetVars()) {
            while (ptr != end && *ptr < var) {
                ptr = std::next(ptr);
            }
            if (ptr == end) {
                vars_res.emplace_back(var);
                continue;
            }

            bool same_number = (var.number == ptr->number);
            VariableDegType new_deg = var.deg - (same_number ? ptr->deg : 0);
            if (new_deg > 0) {
                vars_res.emplace_back(new_deg, var.number);
            }
            if (same_number) {
                ptr = std::next(ptr);
            }
        }
        std::sort(vars_res.begin(), vars_res.end());

        return result;
    }

    // compute the least common multiple for two monoms, excluding coefficients.
    friend Monom GetLcm(const Monom& m1, const Monom& m2) {
        Monom result;

        // Lcm(m, 0) = Lcm(0, m) = 0.
        if (m1.GetCoef() == Field(0) || m2.GetCoef() == Field(0)) {
            return result;
        }
        result.GetCoef() = Field(1);

        auto& vars_res = result.GetVars();
        auto it1 = m1.GetVars().cbegin();
        auto end1 = m1.GetVars().cend();
        auto it2 = m2.GetVars().cbegin();
        auto end2 = m2.GetVars().cend();
        auto push_variable =
            [&vars_res](typename VariablesListType::const_iterator& it) {
                vars_res.emplace_back(*it);

                it = std::next(it);
            };
        // the degree of x_i in lcm equals the maximum degree of x_i in m1 and
        // degree of x_i in m2.
        while (it1 != end1 || it2 != end2) {
            if (it1 == end1) {
                push_variable(it2);
            } else if (it2 == end2) {
                push_variable(it1);
            } else {
                if (*it1 < *it2) {
                    push_variable(it1);
                } else if (*it2 < *it1) {
                    push_variable(it2);
                } else {
                    vars_res.emplace_back(std::max(it1->deg, it2->deg),
                                          it1->number);

                    it1 = std::next(it1);
                    it2 = std::next(it2);
                }
            }
        }
        // maintain the order of variables.
        std::sort(vars_res.begin(), vars_res.end());

        // maintain the degree of result monom.
        for (const auto& [deg, number] : result.GetVars()) {
            result.GetDeg() += deg;
        }

        return result;
    }

    // implement the input operator in accordance with LATEX rules.
    friend std::istream& operator>>(std::istream& in, Monom& m) {
        // read all string.
        std::string curr;
        in >> curr;

        // separate the coefficient from the variables.
        auto pos = curr.find('x');
        // if there are no variables all curr is the coefficient.
        std::string coef_part =
            curr.substr(0, (pos == std::string::npos ? curr.size() : pos));
        std::stringstream coef_str(coef_part);
        // check if there is no coefficient or monom starts just with a +
        // (for instance, +x_1x_2).
        if (pos == 0 || (pos == 1 && curr.front() == '+')) {
            m.GetCoef() = 1;
        } else if (pos == 1 && curr.front() == '-') { // check for -x_1x_2 case.
            m.GetCoef() = Field(-1);
        } else {
            coef_str >> m.GetCoef();
        }

        if (m.GetCoef() == Field(0)) {
            m.GetVars().clear();
            m.GetDeg() = 0;
            return in;
        }

        // read variables if there are ones.
        if (pos != std::string::npos) {
            // delete coefficient in front.
            curr = curr.substr(pos);

            // skip symbol.
            auto skip = [&curr](size_t& i, char c) {
                if (curr[i] == c) {
                    i += 1;
                }
            };
            // read degree or number from given position.
            auto read = [&curr](size_t& i) {
                uint32_t res = 0;
                while (i < curr.size() && isdigit(curr[i])) {
                    res *= 10;
                    res += (curr[i] - '0');
                    i += 1;
                }

                return res;
            };

            // divide curr in x_{number}^{degree} blocks.
            for (size_t i = 0; i < curr.size();) {
                skip(i, 'x'), skip(i, '_'); // skip x_.
                skip(i, '{');               // skip { if number is not a digit.
                VariableNumberType number = read(i); // read number;
                skip(i, '}'); // skip } if number is not a digit.
                if (i == curr.size() ||
                    curr[i] == 'x') { // check if there is no degree explicitly
                                      // written (it implies that degree is
                                      // equal to 1).
                    m.GetDeg() += 1;
                    m.GetVars().emplace_back(1, number);
                    continue;
                }
                skip(i, '^'); // skip ^
                skip(i, '{'); // skip { if number is not a digit.
                VariableDegType deg = read(i); // read degree.
                skip(i, '}');      // skip } if number is not a digit.
                m.GetDeg() += deg; // update the monom's degree.
                m.GetVars().emplace_back(deg, number);
            }
        }

        // maintain the order of variables.
        std::sort(m.GetVars().begin(), m.GetVars().end());

        return in;
    }

    // implement the input operator in accordance with LATEX rules.
    friend std::ostream& operator<<(std::ostream& out, const Monom& m) {
        const auto& coef = m.GetCoef();
        // handle 1 and -1 coefficient value cases.
        if (coef == Field(-1)) {
            if (coef < Field(0)) { // check for Z_p (-1 \equiv p-1 (mod p)).
                out << "-";
            } else {
                if (coef != Field(1)) {
                    out << coef;
                } else if (m.GetDeg() == 0) { // m == "1".
                    out << Field(1);
                }
            }
        } else if (coef != Field(1)) {
            out << coef;
        } else if (coef == Field(1) && m.GetDeg() == 0) { // m == "1".
            out << Field(1);
        }
        auto print = [&out](uint16_t num) {
            if (num >= 10) {
                out << "{" << num << "}";
            } else {
                out << num;
            }
        };
        for (const auto& [deg, number] : m.GetVars()) {
            out << "x_";
            print(number);
            if (deg > 1) {
                out << "^";
                print(deg);
            }
        }

        return out;
    }

  private:
    VariableNumberType not_null_vars() const { return GetVars().size(); }

    VariableDegType monom_deg_ = 0;
    Field coefficient_ = 0;
    VariablesListType variables_;

    inline static std::vector<VariableNumberType> order = {};
};

template <typename Field> class Order : public Monom<Field> {
#define Monom Monom<Field>

    using OrderFunction = std::function<bool(const Monom& m1, const Monom& m2)>;
    using Iterator =
        typename std::vector<typename Monom::Variable>::const_iterator;
    using ReverseIterator =
        typename std::vector<typename Monom::Variable>::const_reverse_iterator;

  public:
    // return order function by its name.
    static OrderFunction GetOrderFunction(const std::string& order) {
        if (order == "lex") {
            return get_lex();
        } else if (order == "grlex") {
            return get_grlex();
        } else if (order == "grevlex") {
            return get_grevlex();
        } else if (order == "invlex") {
            return get_invlex();
        }

        return get_lex();
    }

    // change the monom's order.
    static void InitOrder(const std::string& order) {
        order_function_ = GetOrderFunction(order);
    }

    static const OrderFunction& GetFunction() { return order_function_; }

    // this operator implement to use Order class as a comparator.
    bool operator()(const Monom& m1, const Monom& m2) const {
        return order_function_(m1, m2);
    }

  private:
    // denote x_1^a_1x_2a^2...x_n^a_n as x^a;
    // sdeg(x^a) = a_1+...+a_n.

    // x^a > x^b <=> a > b lexicographically (below lex(a, b)).
    static OrderFunction get_lex() {
        return [&](const Monom& m1, const Monom& m2) -> bool {
            auto [it1, end1, it2, end2] = get_first_not_equal(m1, m2);
            return compare_iterators<Iterator>(it1, end1, it2, end2, false,
                                               true, -1);
        };
    }
    // x^a > x^b <=> sdeg(x^a) > sdeg(x^b) or
    // (sdeg(x^a) == sdeg(x^b) and lex(a, b)).
    static OrderFunction get_grlex() {
        return [&](const Monom& m1, const Monom& m2) -> bool {
            if (m1.GetDeg() != m2.GetDeg()) {
                return m1.GetDeg() > m2.GetDeg();
            }

            auto [it1, end1, it2, end2] = get_first_not_equal(m1, m2);
            return compare_iterators<Iterator>(it1, end1, it2, end2, false,
                                               true, -1);
        };
    }
    // x^a > x^b <=> sdeg(x^a) > sdeg(x^b) or
    // (sdeg(x^a) == sdeg(x^b) and lex(reversed(b), reversed(a))).
    static OrderFunction get_grevlex() {
        return [&](const Monom& m1, const Monom& m2) -> bool {
            if (m1.GetDeg() != m2.GetDeg()) {
                return m1.GetDeg() > m2.GetDeg();
            }

            auto [it1, end1, it2, end2] = get_first_not_equal_rev(m1, m2);
            return compare_iterators<ReverseIterator>(it1, end1, it2, end2,
                                                      true, false, 1);
        };
    }
    // x^a > x^b <=> lex(reversed(a), reversed(b)).
    static OrderFunction get_invlex() {
        return [&](const Monom& m1, const Monom& m2) -> bool {
            auto [it1, end1, it2, end2] = get_first_not_equal_rev(m1, m2);
            return compare_iterators<ReverseIterator>(it1, end1, it2, end2,
                                                      false, true, -1);
        };
    }

    // find first not equal Variables.
    static std::tuple<Iterator, Iterator, Iterator, Iterator>
    get_first_not_equal(const Monom& m1, const Monom& m2) {
        Iterator it1, end1, it2, end2;
        it1 = m1.GetVars().cbegin();
        end1 = m1.GetVars().cend();
        it2 = m2.GetVars().cbegin();
        end2 = m2.GetVars().cend();
        while (it1 != end1 && it2 != end2 && *it1 == *it2) {
            it1 = std::next(it1);
            it2 = std::next(it2);
        }

        return {it1, end1, it2, end2};
    }

    // find last not equal Variables.
    static std::tuple<ReverseIterator, ReverseIterator, ReverseIterator,
                      ReverseIterator>
    get_first_not_equal_rev(const Monom& m1, const Monom& m2) {
        ReverseIterator it1, end1, it2, end2;
        it1 = m1.GetVars().crbegin();
        end1 = m1.GetVars().crend();
        it2 = m2.GetVars().crbegin();
        end2 = m2.GetVars().crend();
        while (it1 != end1 && it2 != end2 && *it1 == *it2) {
            it1 = std::next(it1);
            it2 = std::next(it2);
        }

        return {it1, end1, it2, end2};
    }

    // compare Variables depending on given order.
    template <typename Iterator>
    static bool compare_iterators(Iterator it1, Iterator end1, Iterator it2,
                                  Iterator end2, bool first, bool second,
                                  char sgn) {
        if (it1 == end1) {
            return first;
        }
        if (it2 == end2) {
            return second;
        }

        if (it1->number != it2->number) {
            return typename Monom::Variable(sgn * it1->deg, it1->number) <
                   typename Monom::Variable(sgn * it2->deg, it2->number);
        }
        return sgn * it1->deg < sgn * it2->deg;
    }

    inline static OrderFunction order_function_ = get_lex();
};

template <typename Field> class Polynomial : public Order<Field> {
#define Order Order<Field>

    using ContainerType = std::set<Monom, Order>;

  public:
    Polynomial() = default;
    Polynomial(std::initializer_list<Monom> l) {
        for (const Monom& m : l) {
            if (m.GetCoef() != Field(0)) {
                monoms_.emplace(m);
            }
        }
    }
    Polynomial(const Polynomial& other) { monoms_ = other.monoms_; }
    Polynomial(Polynomial&& other) { monoms_.swap(other.monoms_); }
    Polynomial(std::string&& s) {
        std::stringstream st(s);
        st >> *this;
    }

    Polynomial& operator=(const Polynomial& other) {
        monoms_ = other.monoms_;

        return *this;
    }
    Polynomial& operator=(Polynomial&& other) {
        monoms_.swap(other.monoms_);

        return *this;
    }

    friend bool operator==(const Polynomial& p1, const Polynomial& p2) {
        return p1.GetMonomsList() == p2.GetMonomsList();
    }

    const Monom& GetLeader() const { return *monoms_.cbegin(); }

    ContainerType& GetMonomsList() { return monoms_; }
    const ContainerType& GetMonomsList() const { return monoms_; }

    size_t Size() const { return monoms_.size(); }
    bool IsZero() const { return Size() == 0; }

    void AddMonom(const Monom& m) { GetMonomsList().emplace(m); }
    void AddMonom(Monom&& m) { GetMonomsList().emplace(m); }

    // find first monom that divides by the leader of the other and make the
    // reduction regarding it (mutate *this).
    // return true if reduction was completed.
    bool MakeElementaryReductionInPlace(const Polynomial& other) {
        for (const auto& monom : GetMonomsList()) {
            if (monom.Divides(other.GetLeader())) {
                *this = *this - (monom / other.GetLeader()) * other;

                return true;
            }
        }

        return false;
    }
    // same as previous function but *this doesn't change.
    std::pair<bool, Polynomial>
    MakeElementaryReduction(const Polynomial& other) {
        Polynomial copy(*this);
        auto flag = copy.MakeElementaryReductionInPlace(other);

        return {flag, copy};
    }

    // using the same idea as in operator* for monoms.
    friend Polynomial operator+(const Polynomial& p1, const Polynomial& p2) {
        Polynomial result;

        auto it1 = p1.GetMonomsList().cbegin();
        auto end1 = p1.GetMonomsList().cend();
        auto it2 = p2.GetMonomsList().cbegin();
        auto end2 = p2.GetMonomsList().cend();
        auto push_monom =
            [&result](typename ContainerType::const_iterator& it) {
                result.AddMonom(*it);

                it = std::next(it);
            };
        while (it1 != end1 || it2 != end2) {
            if (it1 == end1) {
                push_monom(it2);
            } else if (it2 == end2) {
                push_monom(it1);
            } else {
                if (Order::GetFunction()(*it1, *it2)) {
                    push_monom(it1);
                } else if (Order::GetFunction()(*it2, *it1)) {
                    push_monom(it2);
                } else {
                    auto sum = *it1 + *it2;
                    // check if monom sums up to zero.
                    if (sum.GetCoef() != Field(0)) {
                        result.AddMonom(sum);
                    }

                    it1 = std::next(it1);
                    it2 = std::next(it2);
                }
            }
        }

        return result;
    }

    // multiply each monom in *this by given monom.
    friend Polynomial operator*(const Polynomial& p, const Monom& m) {
        Polynomial result;

        // check if we multiplying by zero.
        if (m.GetCoef() == 0) {
            return result;
        }

        for (const auto& monom : p.GetMonomsList()) {
            result.AddMonom(monom * m);
        }

        return result;
    }

    // for ability to write polynomial * monom and monom * polynomial.
    friend Polynomial operator*(const Monom& m, const Polynomial& p) {
        return p * m;
    }

    friend Polynomial operator-(const Polynomial& p1, const Polynomial& p2) {
        return p1 + p2 * Monom(Field(-1));
    }

    // define the order of polynomials.
    // based on the strategy of a normal choice.
    friend bool operator<(const Polynomial& p1, const Polynomial& p2) {
        return Order::GetFunction()(p2.GetLeader(), p1.GetLeader());
    }

    // find the S-polynomial for two polynomials accroding to the definition
    // (S(p1, p2) = (l / L(p1)) * p1 - (l / L(p2)) * p2,
    // where l = lcm(L(p1), L(p2))).
    friend Polynomial GetSPolynomial(const Polynomial& p1,
                                     const Polynomial& p2) {
        const auto& l1 = p1.GetLeader();
        const auto& l2 = p2.GetLeader();
        auto lcm = GetLcm(l1, l2);

        return (lcm / l1) * p1 - (lcm / l2) * p2;
    }

    // read polynomial in accordance to Latex rules.
    friend std::istream& operator>>(std::istream& in, Polynomial& p) {
        // read all string.
        std::string s;
        in >> s;

        std::string curr;
        auto flush = [&]() {
            if (!curr.empty()) {
                Monom m;
                std::stringstream ss(curr);
                ss >> m;
                if (!(m == Monom(Field(0)))) {
                    p.AddMonom(m);
                }
            }
            curr.clear();
        };
        // divide string in block of monoms
        // (it should be separated by + or - signs).
        size_t i = 0;
        while (i < s.size()) {
            if (s[i] == '+' || s[i] == '-') {
                flush();
            }
            curr.push_back(s[i]);
            i += 1;
        }
        flush();

        return in;
    }

    // write polynomial in accordance to Latex rules.
    friend std::ostream& operator<<(std::ostream& out, const Polynomial& p) {
        // output 0 if given polynomial is equal to zero.
        if (p.Size() == 0) {
            return out << Field(0);
        }

        bool first = true;
        // print each monom divided by + or -.
        for (const auto& monom : p.GetMonomsList()) {
            if (!first && monom.GetCoef() > 0) {
                out << "+";
            }
            out << monom;

            first = false;
        }

        return out;
    }

    friend bool check_reciprocal(const Polynomial& p1, const Polynomial& p2) {
        auto m1 = p1.GetLeader();
        auto m2 = p2.GetLeader();
        m1.Normalize();
        m2.Normalize();

        return GetLcm(m1, m2) == m1 * m2;
    }

  private:
    ContainerType monoms_;
};

template <typename Field> class PolynomialsSet : public Polynomial<Field> {
#define Polynomial Polynomial<Field>

    using ContainerType = std::vector<Polynomial>;

  public:
    PolynomialsSet() = default;
    PolynomialsSet(std::initializer_list<Polynomial> l) {
        polynomials_.insert(polynomials_.end(), l.begin(), l.end());
    }
    PolynomialsSet(std::string&& s) {
        std::stringstream st(s);
        st >> *this;
    }

    void PushBack(const Polynomial& p) { polynomials_.emplace_back(p); }
    void PushBack(Polynomial&& p) { polynomials_.emplace_back(p); }

    ContainerType& GetPolynomialsList() { return polynomials_; }
    const ContainerType& GetPolynomialsList() const { return polynomials_; }

    // compare sets of polynomials as sets.
    friend bool operator==(const PolynomialsSet& ps1,
                           const PolynomialsSet& ps2) {
        return std::set<Polynomial>(ps1.GetPolynomialsList().begin(),
                                    ps1.GetPolynomialsList().end()) ==
               std::set<Polynomial>(ps2.GetPolynomialsList().begin(),
                                    ps2.GetPolynomialsList().end());
    }

    static void SetOptimization(const std::string& optimization) {
        optimization_ = optimization;
    }

    Polynomial MakeReduction(Polynomial p) const {
        bool reduced = true;
        // iterate till at least one reduction can be completed.
        while (reduced) {
            reduced = false;
            for (const auto& f : GetPolynomialsList()) {
                if (p.MakeElementaryReductionInPlace(f)) {
                    reduced = true;

                    // if gets zero, the remainder won't change anymore.
                    if (p.IsZero()) {
                        return p;
                    }
                }
            }
        }

        return p;
    }

    void Buchberger() {
        // it is not needed to run Buchberger algorithm on a Grobner basis.
        if (is_basis_) {
            return;
        }

        if (optimization_ == "default") {
            buchberger_default();
        } else if (optimization_ == "do_not_repeat_computation") {
            buchberger_do_not_repeat_computation();
        } else if (optimization_ == "skip_reciprocal") {
            buchberger_skip_reciprocal();
        } else if (optimization_ == "fastest_one") {
            buchberger_fastest_one();
        }

        is_basis_ = true;
    }

    // check if polynomial belongs to corresponding ideal with standart
    // algorithm: build the Grobner basis of the ideal and find the
    // remainder of the given polynomial regarding corresponding polynomials
    // set; if the resual is equal to zero than polynomial belongs to the ideal,
    // othewise -- does not belong.
    bool Belongs(const Polynomial& p) {
        Buchberger();

        return MakeReduction(p).IsZero();
    }

    // perform an autoreduction process.
    void AutoReduction() {
        bool reduced = true;
        // iterate till at least one monom of some polynomial divides by
        // some leader of other polynomials.
        while (reduced) {
            reduced = false;

            auto& pols = GetPolynomialsList();
            for (size_t i = 0; i < pols.size(); ++i) {
                for (size_t j = 0; j < pols.size(); ++j) {
                    if (i == j) {
                        continue;
                    }

                    if (pols[i].MakeElementaryReductionInPlace(pols[j])) {
                        reduced = true;
                    }
                }
            }
        }

        // delete all polynomials which was divided by others
        // (polynomials which became equal to 0).
        auto& pols = GetPolynomialsList();
        auto it =
            std::remove_if(pols.begin(), pols.end(),
                           [](const Polynomial& p) { return p.IsZero(); });
        pols.erase(it, pols.end());
    }

    // build minimum Grobner basis according to standar algorithm.
    void BuildMinimumBasis() {
        // make autoreduction and build the Grobner basis.
        AutoReduction();
        Buchberger();

        // make all the leaders' coefficients are equal to 1.
        auto& pols = GetPolynomialsList();
        for (Polynomial& p : pols) {
            p = p * Monom(Field(1) / p.GetLeader().GetCoef());
        }

        // delete all f_i from the system such as exists f_j (i != j) such
        // as L(f_i) divides L(f_j) (the ideal doesn't change in such case).
        bool deleted = true;
        while (deleted) {
            deleted = false;

            for (size_t i = 0; i < pols.size() && !deleted; ++i) {
                for (size_t j = 0; j < pols.size() && !deleted; ++j) {
                    if (i == j) {
                        continue;
                    }

                    if (pols[i].GetLeader().Divides(pols[j].GetLeader())) {
                        pols.erase(pols.begin() + i);
                        deleted = true;
                    }
                }
            }
        }

        // make autoreduction again.
        AutoReduction();
    }

    friend std::istream& operator>>(std::istream& in, PolynomialsSet& ps) {
        // read number of polynomials.
        size_t cnt_pols;
        in >> cnt_pols;

        // read each polynomial.
        ps.GetPolynomialsList().assign(cnt_pols, {});
        for (auto& pol : ps.GetPolynomialsList()) {
            in >> pol;
        }

        return in;
    }

    friend std::ostream& operator<<(std::ostream& out,
                                    const PolynomialsSet& ps) {
        const auto& pols = ps.GetPolynomialsList();
        if (pols.size() == 0) {
            return out << Field(0);
        }
        for (auto it = pols.begin(); it != pols.end(); it = std::next(it)) {
            out << *it;
            if (it != std::prev(pols.end())) {
                out << "; ";
            } else {
                out << ".";
            }
        }

        return out;
    }

  private:
    // standar buchberger algorithm implementation.
    void buchberger_default() {
        auto& pols = GetPolynomialsList();
        bool changed = true;
        // iterate till at least one polynomial adds to the set.
        while (changed) {
            changed = false;

            size_t curr_size = pols.size();
            std::vector<Polynomial> added;
            for (size_t i = 0; i + 1 < curr_size; ++i) {
                for (size_t j = i + 1; j < curr_size; ++j) {
                    auto curr_reduction_res =
                        MakeReduction(GetSPolynomial(pols[i], pols[j]));

                    if (!curr_reduction_res.IsZero()) {
                        changed = true;
                        added.emplace_back(curr_reduction_res);
                    }
                }
            }
            for (const auto& p : added) {
                pols.emplace_back(p);
            }
        }
    }

    // do not make reduction process on polynomials
    // on which it has already been done.
    void buchberger_do_not_repeat_computation() {
        auto& pols = GetPolynomialsList();
        size_t sz = pols.size();
        std::queue<std::pair<size_t, size_t>> q;
        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = i + 1; j < sz; ++j) {
                q.push({i, j});
            }
        }
        // iterate till at least one polynomial adds to the set.
        while (!q.empty()) {
            size_t curr_size = q.size();
            std::vector<Polynomial> added;
            for (size_t it = 0; it < curr_size; ++it) {
                auto [i, j] = q.front();
                q.pop();
                auto curr_reduction_res =
                    MakeReduction(GetSPolynomial(pols[i], pols[j]));

                if (!curr_reduction_res.IsZero()) {
                    for (size_t k = 0; k < curr_size; ++k) {
                        q.push({k, pols.size() + added.size()});
                    }
                    added.emplace_back(curr_reduction_res);
                }
            }
            for (const auto& p : added) {
                pols.emplace_back(p);
            }
        }
    }

    // previous optimization + do not make reduction on
    // polynomials which leader monoms are reciprocal.
    void buchberger_skip_reciprocal() {
        auto& pols = GetPolynomialsList();
        size_t sz = pols.size();
        std::queue<std::pair<size_t, size_t>> q;
        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = i + 1; j < sz; ++j) {
                if (!check_reciprocal(pols[i], pols[j])) {
                    q.push({i, j});
                }
            }
        }
        // iterate till at least one polynomial adds to the set.
        while (!q.empty()) {
            size_t curr_size = q.size();
            std::vector<Polynomial> added;
            for (size_t it = 0; it < curr_size; ++it) {
                auto [i, j] = q.front();
                q.pop();
                auto curr_reduction_res =
                    MakeReduction(GetSPolynomial(pols[i], pols[j]));

                if (!curr_reduction_res.IsZero()) {
                    for (size_t k = 0; k < curr_size; ++k) {
                        if (!check_reciprocal(pols[i], pols[j])) {
                            q.push({k, pols.size() + added.size()});
                        }
                    }
                    added.emplace_back(curr_reduction_res);
                }
            }
            for (const auto& p : added) {
                pols.emplace_back(p);
            }
        }
    }

    // previous two optimizations + criterion based on
    // syzygies of polynomial's leaders.
    void buchberger_fastest_one() {
        auto& pols = GetPolynomialsList();
        size_t sz = pols.size();
        std::queue<std::pair<size_t, size_t>> q;
        std::set<std::pair<size_t, size_t>> used;
        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = i + 1; j < sz; ++j) {
                if (!check_reciprocal(pols[i], pols[j])) {
                    q.push({i, j});
                    used.insert({i, j});
                }
            }
        }

        auto sort_pair = [](size_t i, size_t j) -> std::pair<size_t, size_t> {
            if (i < j) {
                return {i, j};
            }

            return {j, i};
        };
        auto criteria = [&used, &pols, &sort_pair](size_t i, size_t j) -> bool {
            auto lcm = GetLcm(pols[i], pols[j]);
            for (size_t l = 0; l < pols.size(); ++l) {
                if (l == i || l == j) {
                    continue;
                }

                if (used.count(sort_pair(i, l)) ||
                    used.count(sort_pair(j, l))) {
                    continue;
                }

                if (pols[l].GetLeader().Divides(lcm)) {
                    return true;
                }
            }
            return false;
        };

        // iterate till at least one polynomial adds to the set.
        while (!q.empty()) {
            size_t curr_size = q.size();
            std::vector<Polynomial> added;
            for (size_t it = 0; it < curr_size; ++it) {
                auto [i, j] = q.front();
                q.pop();
                auto curr_reduction_res =
                    MakeReduction(GetSPolynomial(pols[i], pols[j]));

                if (!curr_reduction_res.IsZero()) {
                    for (size_t k = 0; k < curr_size; ++k) {
                        if (!check_reciprocal(pols[i], pols[j]) &&
                            !criteria(i, j)) {
                            q.push({k, pols.size() + added.size()});
                            used.insert({k, pols.size() + added.size()});
                        }
                    }
                    added.emplace_back(curr_reduction_res);
                }
            }
            for (const auto& p : added) {
                pols.emplace_back(p);
            }
        }
    }

    bool is_basis_ = false;
    ContainerType polynomials_;

    inline static std::string optimization_ = "default";
};
}; // namespace Grobner