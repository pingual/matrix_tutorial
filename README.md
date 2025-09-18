# matrix_tutorial
Simple implementation of matrices in C for educational purposes about C programming under Linux.

## Documentation
Note that all functions that return a matrix allocate some space for the matrix that needs to be freed with *matrix_destroy* after usage.

### matrix matrix_create(unsigned n1, unsigned n2, scalar v)
Create a matrix of size n1*n2 (n1 lines and n2 columns) and fills it with v.

### void matrix_destroy(matrix m)
Free the space occupied by a matrix. THE MATRIX CANNOT BE USED AFTER.

### matrix matrix_identity(unsigned n)
Returns the idendity matrix of size n.

### scalar *matrix_get(matrix m, unsigned i, unsigned j)
Return a pointer to the element in poisition (i, j) from the matrix m.

### matrix matrix_add(matrix m, matrix n)
Adds the matrix m and the matrix m

### void matrix_print(FILE *f, matrix m)
Prints the matrix m in f

### matrix matrix_multiply(matrix m, matrix n)
Multiplies the matrix m and the matrix n

### matrix matrix_scalar(matrix m, scalar lambda)
Multiplies the matrix m by the scalar lambda

matrix matrix_power(matrix m, int n)
matrix matrix_copy(matrix m)
matrix matrix_inv(matrix m)