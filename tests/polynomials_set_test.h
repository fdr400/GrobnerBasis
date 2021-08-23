#include "../include/GrobnerBasis.h"

#include <gtest/gtest.h>

namespace GrobnerTest {
using Field = double;
using namespace Grobner;

#define Monom Monom<Field>
#define PolynomialsSet PolynomialsSet<Field>

class PolynomialsSetTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
        Monom::InitOrder();
        Order::InitOrder("lex");
    }

    inline static std::vector<std::string> opts = {
        "default", "do_not_repeat_computation", "skip_reciprocal",
        "fastest_one"};
};

// check that PolynomialsSet constructors work properly.
TEST_F(PolynomialsSetTest, TestConstructors) {
    Monom::InitOrder({1, 2, 3});
    Order::InitOrder("lex");
    PolynomialsSet::SetOptimization("default");

    PolynomialsSet ps1;
    ASSERT_EQ(ps1.Size(), 0);

    PolynomialsSet ps2("3 x_1x_2+2x_2x_3 x_1-x_2^2 x_2x_3^2-x_2");

    ASSERT_EQ(
        ps2,
        PolynomialsSet({
            Polynomial({Monom({{1, 1}, {1, 2}}), Monom({{1, 2}, {1, 3}}, 2)}),
            Polynomial({Monom({{1, 1}}), Monom({{2, 2}}, -1)}),
            Polynomial({Monom({{1, 2}, {2, 3}}), Monom({{1, 2}}, -1)}),
        }))
        << ps2 << "\n";
}

TEST_F(PolynomialsSetTest, TestMinimalBasis1) {
    Monom::InitOrder({3, 1, 2});
    Order::InitOrder("lex");
    for (const auto& opt : opts) {
        PolynomialsSet::SetOptimization(opt);

        PolynomialsSet ps("3 x_1x_2+2x_2x_3 x_1-x_2^2 x_2x_3^2-x_2");

        ps.BuildMinimumBasis();
        ASSERT_EQ(ps, PolynomialsSet("3 x_3x_2+0.5x_2^3 x_1-x_2^2 x_2^5-4x_2"))
            << opt << "\n";
    }
}

TEST_F(PolynomialsSetTest, TestMinimalBasis2) {
    Monom::InitOrder({3, 1, 2});
    Order::InitOrder("grlex");
    for (const auto& opt : opts) {
        PolynomialsSet::SetOptimization(opt);

        PolynomialsSet ps(
            "5 x_1^3-2x_1x_2 x_1^2x_2-2x_2^2+x_1 -x_1^2 -2x_1x_2 -2x_2^2+x_1");

        ps.BuildMinimumBasis();
        ASSERT_EQ(ps, PolynomialsSet("3 x_1^2 x_1x_2 x_2^2-0.5x_1"))
            << opt << "\n";
    }
}

TEST_F(PolynomialsSetTest, TestMinimalBasis3) {
    Monom::InitOrder();
    Order::InitOrder("lex");
    for (const auto& opt : opts) {
        PolynomialsSet::SetOptimization(opt);

        PolynomialsSet ps(
            "4 x_1^2+x_2^2+x_3^2 x_1+x_2-x_3 x_2+x_3^2 x_3^4+x_3^3+x_3^2");

        ps.BuildMinimumBasis();
        ASSERT_EQ(ps,
                  PolynomialsSet("3 x_1-x_3^2-x_3 x_2+x_3^2 x_3^4+x_3^3+x_3^2"))
            << opt << "\n";
    }
}
TEST_F(PolynomialsSetTest, TestMinimalBasis4) {
    Monom::InitOrder();
    Order::InitOrder("lex");
    for (const auto& opt : opts) {
        PolynomialsSet::SetOptimization(opt);

        PolynomialsSet ps("3 3x_1-6x_2-2x_3 2x_1-4x_2+4x_4 x_1-2x_2-x_3-x_4");

        ps.BuildMinimumBasis();
        ASSERT_EQ(ps, PolynomialsSet("2 x_1-2x_2+2x_4 x_3+3x_4"))
            << opt << "\n";
    }
}

TEST_F(PolynomialsSetTest, TestBelongs1) {
    Monom::InitOrder();
    Order::InitOrder("lex");
    for (const auto& opt : opts) {
        PolynomialsSet::SetOptimization(opt);

        PolynomialsSet ps("2 x_1^2x_2+2x_3^2 x_2^2-x_2x_3");

        EXPECT_FALSE(ps.Belongs(Polynomial("x_1^3x_3^3+3x_1x_2x_3^3")))
            << opt << "\n";
        EXPECT_TRUE(ps.Belongs(Polynomial("x_1^3x_2^2x_3+2x_1x_2^2x_3^2")))
            << opt << "\n";
    }
}

TEST_F(PolynomialsSetTest, TestBelongs2) {
    Monom::InitOrder();
    Order::InitOrder("grlex");
    for (const auto& opt : opts) {
        PolynomialsSet::SetOptimization(opt);

        PolynomialsSet ps("2 x_1x_3-x_2^2 x_1^3-x_3^2");

        EXPECT_TRUE(ps.Belongs(Polynomial("-4x_1^2x_2^2x_3^2+x_2^6+3x_3^5")))
            << opt << "\n";
        EXPECT_FALSE(ps.Belongs(Polynomial("x_1x_2-5x_2^2+x_1")))
            << opt << "\n";
    }
}
} // namespace GrobnerTest