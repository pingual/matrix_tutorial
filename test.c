#include "utest.h"
#include "matrix.h"

UTEST(matrix_create, simple) {
  matrix m = matrix_create(3, 3, 0.);
  ASSERT_TRUE(m.ok);
}

UTEST_MAIN()  
