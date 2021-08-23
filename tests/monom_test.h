#include "../include/GrobnerBasis.h"

#include <gtest/gtest.h>
#include <numeric>
#include <sstream>
#include <string>

namespace GrobnerTest {
using Field = double;
using namespace Grobner;

#define Monom Monom<Field>

class MonomTest : public ::testing::Test {
    using Container = std::vector<char>;

  protected:
    virtual void SetUp() {
        Monom::InitOrder();
        Order::InitOrder("lex");
    }

    inline static const int number_of_variables = 4;
    inline static const int max_degree_of_variable = 4;
    // max_mask = (max_degree+1)^number_of_variables.
    inline static const int max_mask = 625;

    // check if corresponding to v1 monom divides by the corresponding monom to
    // v2. it implies that v1 and v2 have the same length.
    inline static std::function<bool(Container, Container)> comp =
        [](Container v1, Container v2) {
            for (size_t i = 0; i < v1.size(); ++i) {
                if (v1[i] < v2[i]) {
                    return false;
                }
            }

            return true;
        };
    // build the vector by the given mask.
    inline static std::function<Container(int)> build_vector = [](int mask) {
        Container degs;
        for (int i = 0; i < number_of_variables; ++i) {
            degs.emplace_back(mask % (max_degree_of_variable + 1));
            mask /= (max_degree_of_variable + 1);
        }

        return degs;
    };
    // build monom by the given vector.
    inline static std::function<Monom(Container)> build_monom =
        [](Container v) {
            // zero vector corresponds to Field(1) constant.
            if (std::accumulate(v.begin(), v.end(), 0) == 0) {
                return Monom((Field)1);
            }

            std::vector<Monom::Variable> vars;
            for (size_t i = 0; i < v.size(); ++i) {
                if (v[i] > 0) {
                    vars.emplace_back(v[i], i + 1);
                }
            }

            return Monom(vars);
        };
    // multiply two vector of degrees of monoms
    // (equivalent to ordinary sum of the vectors).
    inline static std::function<Container(Container, Container)>
        multiply_vectors = [](Container v1, Container v2) {
            Container result(v1.size());
            for (size_t i = 0; i < v1.size(); ++i) {
                result[i] = v1[i] + v2[i];
            }

            return result;
        };
    // divide two vector of degrees of monoms
    // (equivalent to ordinary difference of the vectors).
    inline static std::function<Container(Container, Container)>
        divide_vectors = [](Container v1, Container v2) {
            Container result(v1.size());
            for (size_t i = 0; i < v1.size(); ++i) {
                result[i] = v1[i] - v2[i];
            }

            return result;
        };
    // get lcm of two vector of degrees of monoms
    // (equivalent to ordinary coordinate maximum of the vectors).
    inline static std::function<Container(Container, Container)>
        get_lcm_vectors = [](Container v1, Container v2) {
            Container result(v1.size());
            for (size_t i = 0; i < v1.size(); ++i) {
                result[i] = std::max(v1[i], v2[i]);
            }

            return result;
        };
};

// check if all possible constructors work properly.
TEST_F(MonomTest, TestConstructors) {
    Monom m1;
    ASSERT_EQ(m1.GetDeg(), 0);
    ASSERT_EQ(m1.GetCoef(), 0);
    ASSERT_EQ(m1.GetVars().size(), 0);

    Monom m2("7");
    ASSERT_EQ(m2.GetDeg(), 0);
    ASSERT_EQ(m2.GetCoef(), 7);
    ASSERT_EQ(m2.GetVars().size(), 0);

    std::vector<Monom::Variable> tmp{{1, 1}, {1, 2}};

    Monom m3("x_1x_2");
    ASSERT_EQ(m3.GetDeg(), 2);
    ASSERT_EQ(m3.GetCoef(), 1);
    ASSERT_EQ(m3.GetVars(), tmp);

    Monom m4({});
    ASSERT_EQ(m4.GetDeg(), 0);
    ASSERT_EQ(m4.GetCoef(), 1);
    ASSERT_EQ(m4.GetVars().size(), 0);

    Monom m5("7x_1x_2");
    ASSERT_EQ(m5.GetDeg(), 2);
    ASSERT_EQ(m5.GetCoef(), 7);
    ASSERT_EQ(m5.GetVars(), tmp);

    // check if monom with 0 coefficient automatically clear its
    // variables vector.
    Monom m6({{1, 1}, {1, 2}}, 0);
    ASSERT_EQ(m6.GetDeg(), 0);
    ASSERT_EQ(m6.GetCoef(), 0);
    ASSERT_EQ(m6.GetVars().size(), 0);
}

// brute force all existing monoms with {number_of_variables} with maximum
// degree of each variable equal to {max_degree_of_variable}.
TEST_F(MonomTest, TestDividesAndDivideOperator) {
    for (int mask1 = 0; mask1 < max_mask; ++mask1) {
        auto degs1 = build_vector(mask1);
        for (int mask2 = 0; mask2 < max_mask; ++mask2) {
            auto degs2 = build_vector(mask2);

            auto m1 = build_monom(degs1);
            auto m2 = build_monom(degs2);
            EXPECT_EQ(comp(degs1, degs2), m1.Divides(m2))
                << m1 << " " << m2 << "\n";
            if (!(m2 == Monom()) && m1.Divides(m2)) {
                EXPECT_EQ(build_monom(divide_vectors(degs1, degs2)), m1 / m2)
                    << m1 << " " << m2 << "\n";
            }
        }
    }
}

// same as TestDividesAndDivideOperator but for operator*.
TEST_F(MonomTest, TestOperatorsBruteForce) {
    for (int mask1 = 0; mask1 < max_mask; ++mask1) {
        auto degs1 = build_vector(mask1);
        for (int mask2 = 0; mask2 < max_mask; ++mask2) {
            auto degs2 = build_vector(mask2);

            auto m1 = build_monom(degs1);
            auto m2 = build_monom(degs2);
            EXPECT_EQ(build_monom(multiply_vectors(degs1, degs2)), m1 * m2)
                << m1 << " " << m2 << "\n";
        }
    }
}

// test if input and output of monom works properly.
TEST_F(MonomTest, TestInOutOperators) {
    auto test_helper = [&](std::string read, std::string need) {
        std::stringstream s1(read);
        Monom m1;
        s1 >> m1;
        std::stringstream s2;
        s2 << m1;
        ASSERT_EQ(s2.str(), need);
    };

    test_helper("7", "7");
    test_helper("7x_1", "7x_1");
    test_helper("x_2", "x_2");
    test_helper("x_1x_2", "x_1x_2");
    test_helper("-x_1^2x_2", "-x_1^2x_2");
    test_helper("-3x_{10}^{10}x_{47}^{23}x_{1000}^{189}",
                "-3x_{10}^{10}x_{47}^{23}x_{1000}^{189}");

    // for Z_2.
    /*
    test_helper("7", "1");
    test_helper("7x_1", "x_1");
    test_helper("x_2", "x_2");
    test_helper("x_1x_2", "x_1x_2");
    test_helper("-x_1^2x_2", "x_1^2x_2");
    test_helper("-3x_{10}^{10}x_{47}^{23}x_{1000}^{189}",
                "x_{10}^{10}x_{47}^{23}x_{1000}^{189}");
    */
}

// check if that order for variables works.
TEST_F(MonomTest, TestOrderOfVariables) {
    {
        Monom m("x_2x_1x_3");
        auto tmp = std::vector<Monom::Variable>{{1, 1}, {1, 2}, {1, 3}};
        ASSERT_TRUE(m.GetVars() == tmp);
    }
    {
        Monom::InitOrder({3, 1, 2});
        Monom m("x_1x_2x_3");
        auto tmp = std::vector<Monom::Variable>{{1, 3}, {1, 1}, {1, 2}};
        ASSERT_TRUE(m.GetVars() == tmp);
    }
}

// same as TestDividesAndDivideOperator but for GetLcm.
TEST_F(MonomTest, TestLcmBruteForce) {
    Monom::InitOrder();
    for (int mask1 = 0; mask1 < max_mask; ++mask1) {
        auto degs1 = build_vector(mask1);
        for (int mask2 = 0; mask2 < max_mask; ++mask2) {
            auto degs2 = build_vector(mask2);

            auto m1 = build_monom(degs1);
            auto m2 = build_monom(degs2);
            ASSERT_EQ(build_monom(get_lcm_vectors(degs1, degs2)),
                      GetLcm(m1, m2))
                << m1 << " " << m2 << "\n";
        }
    }
}
} // namespace GrobnerTest