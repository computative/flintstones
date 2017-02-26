#include "gc.h"
#include <algorithm>
#include <string>
#include <omp.h>

using namespace std;

gc::gc(string polynomial, int points)
{
    for (int i=0; polynomial[i]; i++) polynomial[i] = tolower(polynomial[i]);
    type = polynomial;
    n = points;
}

double gc::quad4(double (*f)(double ,double ,double ,double ))
{
    if (type == "hermite")
        return hermite4(f);
}

double gc::hermite4(double (*f)(double ,double ,double ,double ))
{
    double *x = new double [n+1];
    double *w = new double [n];
    int i,its,j,m,k,l;
    double p1,p2,p3,pp,z,z1;
    double Epsilon = 3.0e-14, PIM4 = 0.7511255444649425;
    int MaxIterations = 10;
    m=(n+1)/2;

    for (i=1;i<=m;i++) {
      if (i == 1) {
        z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
      } else if (i == 2) {
        z -= 1.14*pow((double)n,0.426)/z;
      } else if (i == 3) {
        z=1.86*z-0.86*(x[0]);
      } else if (i == 4) {
        z=1.91*z-0.91*(x[1]);
      } else {
        z=2.0*z-x[i-3];
      }
      for (its=1;its<=MaxIterations;its++) {
        p1=PIM4;
        p2=0.0;
        for (j=1;j<=n;j++) {
          p3=p2;
          p2=p1;
          p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
        }
        pp=sqrt((double)2*n)*p2;
        z1=z;
        z=z1-p1/pp;
        if (fabs(z-z1) <= Epsilon) break;
      }
      if (its > MaxIterations) cout << "too many iterations in Hermite quadrature" << endl;
      x[i-1]=z;
      x[n-i] = -z;
      w[i-1]=2.0/(pp*pp);
      w[n-i]=w[i-1];
    }

    double Integral = 0.0;
    # pragma omp parallel default(shared) private (i, j, k, l) reduction(+:Integral)
    {
    # pragma omp for
      for ( i = 0;  i < n; i++){
        for ( j = 0;  j < n; j++){
          for ( k = 0;  k < n; k++){
            for ( l = 0;  l < n; l++){
              Integral += w[i]*w[j]*w[k]*w[l]*f(x[i],x[j],x[k],x[l]);
            }
          }
        }
      }
    }
    delete [] x;
    delete [] w;
    return Integral;
}
