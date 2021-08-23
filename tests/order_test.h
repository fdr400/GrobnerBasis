#include "../include/GrobnerBasis.h"

#include <gtest/gtest.h>

namespace GrobnerTest {
using Field = double;
using namespace Grobner;

#define Monom Monom<Field>

class OrderTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
        Monom::InitOrder();
        Order::InitOrder("lex");
    }
};

// check if lex order works properly.
TEST_F(OrderTest, TestLex) {
    auto lex = Order::GetOrderFunction("lex");

    Monom x1("x_1");
    Monom x2("x_2");
    ASSERT_TRUE(lex(x1, x2));

    Monom x1x2("x_1x_2");
    ASSERT_TRUE(lex(x1x2, x1));
    ASSERT_FALSE(lex(x1, x1x2));
    ASSERT_TRUE(lex(x1x2, x2));
    ASSERT_FALSE(lex(x2, x1x2));

    Monom m1("x_1x_2^2");
    Monom m2("x_2^3x_3^4");
    ASSERT_TRUE(lex(m1, m2));
    Monom m3("x_1^3x_2^2x_3^4");
    Monom m4("x_1^3x_2^2x_3");
    ASSERT_TRUE(lex(m3, m4));
}

// check if grlex order works properly.
TEST_F(OrderTest, TestGrlex) {
    auto grlex = Order::GetOrderFunction("grlex");

    Monom x1("x_1");
    Monom x2("x_2");
    ASSERT_TRUE(grlex(x1, x2));

    Monom m1("x_1x_2^2x_3^3");
    Monom m2("x_1^3x_2^2");
    ASSERT_TRUE(grlex(m1, m2));

    Monom m3("x_1x_2^2x_3^4");
    Monom m4("x_1x_2x_3^5");
    ASSERT_TRUE(grlex(m3, m4));

    Monom m5("x_1^5x_2x_3");
    Monom m6("x_1^4x_2x_3^2");
    ASSERT_TRUE(grlex(m5, m6));
}

// check if grevlex order works properly.
TEST_F(OrderTest, TestGrevlex) {
    auto grevlex = Order::GetOrderFunction("grevlex");

    Monom x1("x_1");
    Monom x2("x_2");
    ASSERT_TRUE(grevlex(x1, x2));

    Monom m1("x_1^4x_2^7x_3");
    Monom m2("x_1^4x_2^2x_3^3");
    ASSERT_TRUE(grevlex(m1, m2));

    Monom m3("x_1x_2^5x_3^2");
    Monom m4("x_1^4x_2x_3^3");
    ASSERT_TRUE(grevlex(m3, m4));

    Monom m5("x_1^5x_2x_3");
    Monom m6("x_1^4x_2x_3^2");
    ASSERT_TRUE(grevlex(m5, m6));
}

// here is no check for invlex because there is no ideological difference
// between lex and invlex.
}; // namespace GrobnerTest