#ifndef JACOBI_H
#define JACOBI_H

#include <stdexcept>

#include <eigen3/Eigen/Core>
#include "typedefs.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

/*
 Linear equation solution by Jacobi iterative technique.
 Input: Number of equations; a matrix A; a vector b; initial guess;
        maximum number of iterations; error tolerance
 Output: Approximate solution of the system or an error message

 Ref: Burden, R., Faires, D., Burden, A. 2016. Numerical Analysis, 10rd ed.
      Cengage: Somewhere awful.
 */
VectorXd jacobi(const MatrixXd &A, const VectorXd &b, VectorXd &XO, const Int N, const Doub TOL)
{
  /* First draft: Too much copying of vectors */
  Int k = 0;

  /* Because who doesn't like explosions */
  if (A.rows() > b.size()) throw std::runtime_error("System overdetermined.");
  if (A.rows() < b.size()) throw std::runtime_error("System underdetermined.");

  Int n = b.size() - 1;
  VectorXd x = XO;
  Doub sum;
  VectorXd temp;

  while (k <= N)
    {
      for (Int i = 0; i <= n; i++)
        {
          sum = 0;
          /* Using access with range checking here cf. Eigen documentation */
          for (Int j = 0; j <=n; j++)
            if (j != i)
              sum = sum + A(i, j)*XO(j);
          x(i) = (b(i) - sum)/A(i, i);
        }
      temp = x - XO;
      if (temp.norm() < TOL)
        return x;
      XO = x;
    }
  /* I am like Oprah, but with exceptions */
  throw std::runtime_error("Maximum number of iterations exceeded");
}

#endif
