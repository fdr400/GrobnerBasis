#include <gtest/gtest.h>

#include "monom_test.h"
#include "order_test.h"
#include "polynomial_test.h"
#include "polynomials_set_test.h"

#include "modular_test.h"

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}