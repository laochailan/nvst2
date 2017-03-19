#ifndef CG_H
#define CG_H

// This is an implementation of the conjugate gradient method
// for solving the equation op(x) = b
//
// res, x, x0, tmp1, tmp2 are n-sized vectors.
//
int conjugate_gradient(double *x0, void (*op)(double *res, double *x, int n, void *arg), double *b, int n, void *arg, double tol, int maxiter, double *tmp1, double *tmp2);

// laplace2D is a conjugate_gradient compatible 1st order laplace operator in
// two dimensions with unit grid spacing. The outermost cells on the grid are
// treated as constant border conditions.
//
void laplace2D(double *res, double *x, int n, void *arg);
struct laplace2Dargs {
	int w;
};

void test_conjugate_gradient(void);
void test_conjugate_gradient2(void);

#endif
