#include "../include/modular.h"

#include <gtest/gtest.h>

namespace ModularTest {
TEST(ModularTest, TestConstructorsModulo2) {
    Modular::SetModulo(2);
    for (int i = -10, parity = 0; i <= 10; ++i, parity = 1 - parity) {
        Modular m(i);
        ASSERT_EQ(m.GetValue(), parity) << m.GetValue() << " " << i << "\n";
    }
}

TEST(ModularTest, TestConstructorsModulo3) {
    Modular::SetModulo(3);
    for (int i = -30, remainder = 0; i <= 30;
         ++i, remainder = (remainder + 1 == 3 ? 0 : remainder + 1)) {
        Modular m(i);
        ASSERT_EQ(m.GetValue(), remainder) << m.GetValue() << " " << i << "\n";
    }
}

TEST(ModularTest, TestAdditionOperators) {
    Modular::SetModulo(2);
    for (int i = -10, par_i = 0; i <= 10; ++i, par_i = 1 - par_i) {
        for (int j = -10, par_j = 0; j <= 10; ++j, par_j = 1 - par_j) {
            Modular m1(i), m2(j);

            ASSERT_EQ(m1 + m2, (int)(par_i != par_j))
                << m1 << " " << m2 << " " << i << " " << j << "\n";
            ASSERT_EQ(m1 - m2, (int)(par_i != par_j))
                << m1 << " " << m2 << " " << i << " " << j << "\n";
        }
    }
}

TEST(ModularTest, TestProductOperator) {
    Modular::SetModulo(2);
    for (int i = -10, par_i = 0; i <= 10; ++i, par_i = 1 - par_i) {
        for (int j = -10, par_j = 0; j <= 10; ++j, par_j = 1 - par_j) {
            Modular m1(i), m2(j);

            ASSERT_EQ(m1 * m2, (int)(par_i + par_j == 2))
                << m1 << " " << m2 << " " << i << " " << j << "\n";
        }
    }
}

TEST(ModularTest, TestGetInverse) {
    std::vector<int> prime_modulos = {2, 3, 5, 7, 11, 13, 17};

    for (const auto& mod : prime_modulos) {
        Modular::SetModulo(mod);

        for (int i = 1; i <= 100; ++i) {
            Modular m(i);
            if (m == Modular(0)) {
                continue;
            }

            ASSERT_EQ(m.GetInverse() * m, 1) << m << " " << mod << "\n";
        }
    }
}

TEST(ModularTest, TestDivideOperator) {
    std::vector<int> prime_modulos = {2, 3, 5, 7, 11, 13, 17};

    for (const auto& mod : prime_modulos) {
        Modular::SetModulo(mod);

        for (int i = 0; i <= 100; ++i) {
            for (int j = 1; j <= 100; ++j) {
                Modular m1(i), m2(j);
                if (m2 == Modular(0)) {
                    continue;
                }

                auto res = m1 / m2;
                ASSERT_EQ(m1, res * m2) << i << " " << j << " " << mod << "\n";
            }
        }
    }
}

TEST(ModularTest, TestInOutOperators) {
    auto test_helper = [&](std::string read, std::string need) {
        std::stringstream s1(read);
        Modular m1;
        s1 >> m1;
        std::stringstream s2;
        s2 << m1;
        ASSERT_EQ(s2.str(), need);
    };

    Modular::SetModulo(2);
    test_helper("0", "0");
    test_helper("1", "1");
    test_helper("-1", "1");
    test_helper("-2", "0");
    test_helper("2", "0");

    Modular::SetModulo(3);
    test_helper("0", "0");
    test_helper("1", "1");
    test_helper("-1", "2");
    test_helper("-2", "1");
    test_helper("2", "2");
    test_helper("3", "0");
}
} // namespace ModularTest