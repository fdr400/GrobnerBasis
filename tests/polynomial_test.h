#include "../include/GrobnerBasis.h"

#include <gtest/gtest.h>

namespace GrobnerTest {
using Field = double;
using namespace Grobner;

class PolynomialTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
        Monom::InitOrder();
        Order::InitOrder("lex");
    }

    // check that monoms in polynomial are in the right order.
    bool compare(const Polynomial& p, std::vector<Monom>&& v) {
        return std::equal(v.cbegin(), v.cend(), p.GetMonomsList().cbegin());
    }
};

// check if polynomial constructors work properly.
TEST_F(PolynomialTest, TestConstructors) {
    Order::InitOrder("lex");

    Polynomial p_empty;
    ASSERT_EQ(p_empty.Size(), 0);

    Monom m1("x_1");
    Monom m2("x_2");
    Monom m3("x_2^2");
    Polynomial p({m3, m1, m2});

    ASSERT_TRUE(compare(p, std::vector<Monom>{m1, m3, m2})) << p << "\n";
}

// check that Polynomial properly sorts monoms inside.
TEST_F(PolynomialTest, TestDifferentOrders) {
    Monom m1("-5x_1^3");
    Monom m2("7x_1^2x_3^2");
    Monom m3("4x_1x_2^2x_3");
    Monom m4("4x_3^2");

    Order::InitOrder("lex");
    Polynomial p1({m1, m2, m3, m4});
    Order::InitOrder("grlex");
    Polynomial p2({m1, m2, m3, m4});
    Order::InitOrder("grevlex");
    Polynomial p3({m1, m2, m3, m4});

    ASSERT_TRUE(compare(p1, std::vector<Monom>{m1, m2, m3, m4})) << p1 << "\n";
    ASSERT_TRUE(compare(p2, std::vector<Monom>{m2, m3, m1, m4})) << p2 << "\n";
    ASSERT_TRUE(compare(p3, std::vector<Monom>{m3, m2, m1, m4})) << p3 << "\n";
}

// check that GetSPolynomial function properly finds the S-polynomials for each
// pair of above polynomials.
TEST_F(PolynomialTest, TestSPolynomial) {
    Order::InitOrder("lex");
    Polynomial p1("2x_1x_2+4x_1x_3+x_2x_3^2");
    Polynomial p2("4x_1x_3^2+x_2x_3^3-4");
    Polynomial p3("x_2^2x_3^3-4x_2-8x_3");

    ASSERT_EQ(GetSPolynomial(p1, p2) * Monom(4),
              Polynomial("8x_1x_3^3-x_2^2x_3^3+2x_2x_3^4+4x_2"));
    ASSERT_EQ(GetSPolynomial(p1, p3) * Monom(2),
              Polynomial("4x_1x_2x_3^4+8x_1x_2+16x_1x_3+x_2^2x_3^5"));
    ASSERT_EQ(GetSPolynomial(p2, p3) * Monom(4),
              Polynomial("16x_1x_2+32x_1x_3+x_2^3x_3^4-4x_2^2x_3"));
}

// check if MakeElementaryReduction function properly reduces the polynomial
// regrading other polynomials.
TEST_F(PolynomialTest, TestElementaryReduction) {
    Order::InitOrder("lex");
    Polynomial p("8x_1x_3^3-x_2^2x_3^3+2x_2x_3^4+4x_2");

    Polynomial f1("4x_1x_3^2+x_2x_3^3-4");
    Polynomial f2("x_2^2x_3^3-4x_2-8x_3");

    p.MakeElementaryReductionInPlace(f1);
    ASSERT_EQ(p, Polynomial("-x_2^2x_3^3+4x_2+8x_3")) << p << "\n";
    p.MakeElementaryReductionInPlace(f2);
    ASSERT_EQ(p, Polynomial("0"));
}
}; // namespace GrobnerTest