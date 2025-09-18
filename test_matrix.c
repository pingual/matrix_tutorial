#include "utest.h"
#include "matrix.h"

UTEST(matrix_create, simple) {
  matrix m = matrix_create(3, 5, 0.);
  ASSERT_TRUE(m.ok && m.n1==3 && m.n2==5);
}

UTEST(matrix_create, too_large) {
  matrix m = matrix_create((unsigned)-1, (unsigned)-1, 0.);
  ASSERT_FALSE(m.ok);
}

UTEST(matrix_create, zero) {
  matrix m = matrix_create(0, 0, 0.);
  ASSERT_TRUE(m.ok && m.n1==0 && m.n2==0);
}

UTEST(matrix_identity, simple) {
  matrix m = matrix_identity(42);
  ASSERT_TRUE(m.ok && m.n1==42 && m.n2==42);
}

UTEST(matrix_identity, too_large) {
  matrix m = matrix_identity((unsigned) -1);
  ASSERT_TRUE(!m.ok);
}

UTEST(matrix_identity, zero) {
  matrix m = matrix_identity(0);
  ASSERT_TRUE(m.ok && m.n1==0 && m.n2==0);
}

UTEST(matrix_destroy, simple) {
  matrix m = matrix_create(2, 5, 0.);
  matrix_destroy(m);
  ASSERT_TRUE(1);
}

UTEST(matrix_destroy, wrong) {
  matrix m;
  m.ok=false;
  m.data=(scalar*)0xdeadbeef;
  matrix_destroy(m);
  ASSERT_TRUE(1);
}

UTEST(matrix_get, simple) {
  matrix m = matrix_identity(12);
  ASSERT_TRUE(*matrix_get(m,3,3)==1.);
}

UTEST(matrix_get, wrong) {
  matrix m;
  m.ok=false;
  ASSERT_TRUE(matrix_get(m, 0, 0)==NULL);
}

UTEST(add_matrix, simple) {
  matrix m1 = matrix_create(2, 2, 1.);
  matrix m2 = matrix_create(2, 2, 1.);
  matrix m3 = matrix_add(m1, m2);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      ASSERT_TRUE(*matrix_get(m3, i, j) == 2.);
    }
  }
}

UTEST(matrix_multiply, simple) {
  matrix m1 = matrix_create(2, 2, 2.);
  matrix m2 = matrix_create(2, 2, 2.);
  matrix m3 = matrix_multiply(m1, m2);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      ASSERT_TRUE(*matrix_get(m3, i, j) == 8.);
    }
  }
}

UTEST(matrix_scalar, simple) {
  matrix m1 = matrix_create(2, 2, 2.);
  matrix m2 = matrix_scalar(m1, 3.);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      ASSERT_TRUE(*matrix_get(m2, i, j) == 6.);
    }
  }
}

UTEST(matrix_power, simple) {
  matrix m1 = matrix_create(2, 2, 2.);
  matrix m2 = matrix_power(m1, 2);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      ASSERT_TRUE(*matrix_get(m2, i, j) == 8.);
    }
  }
}

UTEST(matrix_det, simple) {
  matrix m = matrix_create(3, 3, 0.);
  *matrix_get(m, 0, 0) = 1.;
  *matrix_get(m, 0, 1) = 1.;
  *matrix_get(m, 0, 2) = 1.;
  *matrix_get(m, 1, 0) = 1.;
  *matrix_get(m, 1, 1) = 2.;
  *matrix_get(m, 1, 2) = 3.;
  *matrix_get(m, 2, 0) = 3.;
  *matrix_get(m, 2, 1) = 4.;
  *matrix_get(m, 2, 2) = 6.;
  ASSERT_NEAR(1., matrix_det(m), 0.00001);
}


UTEST_MAIN()  
