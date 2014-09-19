#include "src/Satsuma/SatsumaProt.h"
#include "gtest/gtest.h"

// Tests factorial of negative numbers.
TEST(SatsumaProtTest, Basic) {
  // This test is named "Negative", and belongs to the "FactorialTest"
  // test case.

  SatsumaProtParams params(false, false, false, 4, 1, 50, 0, 0, 50, false, false, 0, 2);
//  SatsumaProt protAligner(REF, DB, QUERY, params);
//  protAligner.alignAll();

  EXPECT_EQ(1, 1);
  EXPECT_EQ(-1, -1);
  EXPECT_TRUE(10 > 0);

}

