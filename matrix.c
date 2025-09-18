#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "matrix.h"

#define PRINT_PRECISION "6"
#define PRINT_DECIMAL_PRECISION "2"

matrix matrix_create(unsigned n1, unsigned n2, scalar v)
{
  matrix m={n1,n2,true,NULL};
  if(n1==0 || n2==0)
    return m;
  else if(SIZE_MAX / n1 <= n2)
    m.ok = false;
  else if (!(m.data=calloc(((size_t) n1)*n2, sizeof(scalar))))
    m.ok = false;
  else {
    for(unsigned i=0; i<n1; ++i)
      for(unsigned j=0; j<n2; ++j)
        *matrix_get(m,i,j) = v;
  }
 
  return m;
}

matrix matrix_identity(unsigned n)
{
  matrix m = matrix_create(n, n, 0.);
  if(!m.ok)
    return m;

  for(unsigned i=0; i<n; ++i)
    *matrix_get(m, i, i) = 1.;
  return m;
}

void matrix_destroy(matrix m)
{
  if(m.ok) {
    m.ok = false;
    free(m.data);
  }
}

scalar *matrix_get(matrix m, unsigned i, unsigned j)
{
  if(!m.ok || i>m.n1 || i>m.n2)
    return NULL;

  return &m.data[i*m.n2+j];
}

matrix matrix_add(matrix m, matrix n)
{
  matrix res={0,0,false,NULL};

  if(m.n1!=n.n1 || m.n2!=n.n2 || !m.ok || !n.ok)
    return res;

  res=matrix_create(m.n1, m.n2, 0.);
  for(unsigned i=0; i<m.n1; ++i)
    for(unsigned j=0; j<m.n2; ++j)
      *matrix_get(res, i, j) = *matrix_get(m, i, j) + *matrix_get(n, i, j);

  return res;
}

void matrix_print(FILE *f, matrix m)
{
  if(!m.ok)
    fprintf(f, "Invalid matrix\n");
  else {
    for(unsigned i=0; i<m.n1; ++i) {
      for(unsigned j=0; j<m.n2; ++j)
        fprintf(
            f, 
            "%"PRINT_PRECISION"."PRINT_DECIMAL_PRECISION"f ",
            *matrix_get(m, i, j));
      fprintf(f, "\n");
    }
  }
}

matrix matrix_multiply(matrix m, matrix n) {
    if (m.n2 != n.n1) return matrix_create(0, 0, 0.);

    matrix res = matrix_create(m.n1, n.n2, 0.);
    for (int i = 0; i < m.n1; i++) {
        for (int j = 0; j < n.n2; j++) {
            for (int k = 0; k < m.n2; k++) {
                *matrix_get(res, i, j) += (*matrix_get(m, i, k))*(*matrix_get(n, k, j));
            }
        }
    }

    return res;
}

matrix matrix_scalar(matrix m, double lambda) {
    matrix res = matrix_create(m.n1, m.n2, 0.);
    for (int i = 0; i < m.n1; i++) {
        for (int j = 0; j < m.n2; j++) {
            *matrix_get(res, i, j) = *matrix_get(m, i, j)*lambda;
        }
    }
    return res;
}

matrix matrix_power(matrix m, int n) {
    if (n == 0) return matrix_identity(m.n1);

    int maxbit = 0;
    for (int val = 1; val <= n; val *= 2) maxbit++;
    
    matrix* pows = malloc(maxbit*sizeof(matrix));
    pows[0] = m;
    for (int i = 1; i < maxbit; i++) {
        pows[i] = matrix_multiply(pows[i - 1], pows[i - 1]);
    }

    matrix res = matrix_identity(m.n1);
    int decomp = n;
    int pos = 0;
    while (decomp > 0) {
        if (decomp % 2 == 1) {
            matrix temp = matrix_multiply(res, pows[pos]);
            matrix_destroy(res);
            res = temp;
        }
        decomp = decomp / 2;
        pos++;
    }

    for (int i = 1; i < maxbit; i++) {
        matrix_destroy(pows[i]);
    }
    free(pows);

    return res;
}

matrix matrix_copy(matrix m) {
    matrix res = matrix_create(m.n1, m.n2, 0.);
    for (int i = 0; i < m.n1; i++) {
        for (int j = 0; j < m.n2; j++) {
            *matrix_get(res, i, j) = *matrix_get(m, i, j);
        }
    }
    return res;
}

void echange_lignes(matrix m, int l1, int l2) {
    for (int i = 0; i < m.n2; i++) {
        double temp = *matrix_get(m, l1, i);
        *matrix_get(m, l1, i) = *matrix_get(m, l2, i);
        *matrix_get(m, l2, i) = temp;
    }
}

void transvection(matrix m, int ligne_source, int ligne_dest, double lambda) {
    for (int i = 0; i < m.n2; i++) {
        *matrix_get(m, ligne_dest, i) += lambda*(*matrix_get(m, ligne_source, i));
    }
}

void dilatation_ligne(matrix m, int ligne, double lambda) {
    for (int i = 0; i < m.n2; i++) {
        *matrix_get(m, ligne, i) *= lambda;
    }
}



matrix matrix_inv(matrix m) {
    if (m.n2 != m.n1) return matrix_create(0, 0, 0.);

    matrix a_red = matrix_copy(m);
    matrix res = matrix_identity(m.n1);

    //descente
    for (int i = 0; i < m.n1; i++) {
        int j = i;
        while ((j < m.n2) && (*matrix_get(a_red, j, i) == 0.)) {
            j++;
        }
        if (j >= m.n2) return matrix_create(0, 0, 0.);
        echange_lignes(a_red, i, j);
        echange_lignes(res, i, j);
        dilatation_ligne(res, i, 1/(*matrix_get(a_red, i, i)));
        dilatation_ligne(a_red, i, 1/(*matrix_get(a_red, i, i)));
        *matrix_get(a_red, i, i) = 1.;

        for (int j = i+1; j < m.n2; j++) {
            transvection(res, i, j, -(*matrix_get(a_red, j, i)));
            transvection(a_red, i, j, -(*matrix_get(a_red, j, i)));
            *matrix_get(a_red, j, i) = 0.;
        }
    }

    //montée
    for (int i = m.n1 - 1; i >= 0; i--) {
        for (int j = 0; j < i; j++) {
            transvection(res, i, j, -(*matrix_get(a_red, j, i)));
            transvection(a_red, i, j, -(*matrix_get(a_red, j, i)));
            *matrix_get(a_red, j, i) = 0.;
        }
    }

    matrix_destroy(a_red);
    return res;
}