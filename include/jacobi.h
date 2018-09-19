#ifndef JACOBI_H
#define JACOBI_H

#include <stdexcept>
#include <armadillo>

using namespace arma;

/*
 Linear equation solution by Jacobi iterative method.
 Input: Number of equations; a matrix A; a vector b; initial guess;
        maximum number of iterations; error tolerance
 Output: Approximate solution of the system or an error message

 Ref: Burden, R., Faires, D., Burden, A. 2016. Numerical Analysis, 10rd ed.
      Cengage: Somewhere awful.
 */
vec jacobi(const mat &A, const mat &b, vec &XO, const int N, const double TOL)
{
  /* First draft: Too much copying of vectors */
  int k = 0;

  int n = b.size() - 1;
  vec x = XO;
  double sum;
  vec temp;

  while (k <= N)
    {
      for (int i = 0; i <= n; i++)
        {
          sum = 0;
          /* Using access with range checking here cf. Eigen documentation */
          for (int j = 0; j <=n; j++)
            if (j != i)
              sum = sum + A(i, j)*XO(j);
          x(i) = (b(i) - sum)/A(i, i);
        }
      if (norm(x - XO, "inf") < TOL)
        return x;
      XO = x;
    }
  /* I am like Oprah, but with exceptions */
  throw std::runtime_error("Maximum number of iterations exceeded");
}

#endif
