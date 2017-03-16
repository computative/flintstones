#include <iostream>
#include <fstream>
#include <armadillo>
#include <map>
#include <omp.h>
#include "basis.h"
#include "gc.h"
#include "math.h"
#include "Coulomb_Functions.hpp"

using namespace std;
using namespace arma;


int m(int a);
int n(int a);
int E(int a);
int A(int m, int l);
int B(int N);
void invA(int A,int &k, int &l);
int sigma(int a);
int kronecker(int a, int b);
double onebody(int n, int m, double hw);
int merge(int alpha, int beta, int N);


bool exists(const std::map<int, double>& m, int key)
{
    return m.find(key) != m.end();
}

template <typename... Args, typename ...Krgs>
bool exists(const std::map<int, Args...>& m, int key, Krgs ...keys)
{
    auto it = m.find(key);
    return it != m.end() && exists(it->second, keys...);
}


int main(int, char *[])
{
    int F = 6; // fermi level
    double hw = 1;
    int R = 4; // num shells
    int N = R*(R+1); // num of states less than R
    double M, Ms;
    for (int alpha = 0; alpha < N; alpha++) {
        for (int beta = 0; beta < N; beta++ ) {
            for (int gamma = 0; gamma < N; gamma++) {
                for (int delta = 0; delta < N; delta ++) {
                    int ml1 = m(alpha+1);
                    int ml2 = m(beta+1);
                    int ml3 = m(gamma+1);
                    int ml4 = m(delta+1);
                    if (sigma(alpha+1) + sigma(beta+1) == sigma(gamma+1) + sigma(delta+1) and ml1 + ml2 == ml3 + ml4) {
                        int n1 = n(alpha+1);
                        int n2 = n(beta+1);
                        int n3 = n(gamma+1);
                        int n4 = n(delta+1);
                        M = ml1 + ml2;
                        Ms = sigma(alpha+1) + sigma(beta+1);

                        double A = kronecker(sigma(alpha+1), sigma(gamma+1))*kronecker(sigma(beta+1), sigma(delta+1))*Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4)
                                - kronecker(sigma(alpha+1), sigma(delta+1))*kronecker(sigma(gamma+1), sigma(beta+1))*Coulomb_HO(hw, n1, ml1, n2, ml2, n4, ml4, n3, ml3);

                        double B = kronecker(sigma(gamma+1), sigma(alpha+1))*kronecker(sigma(delta+1), sigma(beta+1))*Coulomb_HO(hw, n3, ml3, n4, ml4, n1, ml1, n2, ml2)
                                - kronecker(sigma(delta+1), sigma(alpha+1))*kronecker(sigma(beta+1), sigma(gamma+1))*Coulomb_HO(hw, n4, ml4, n3, ml3, n1, ml1, n2, ml2);
                        cout << A - B << endl;
                    }
                }
            }
        }
    }
    return 0;
}

int A(int m, int l) { // map from 2-touples to a section of the natural numbers
    double G = (m+l);
    return ceil(G/2.*(G/2 - 1)) + G - max(m,l);
}

int B(int N) { // hoyeste index fra funksjonen A. N er antallet tilstander
    int m = floor((N+1)/2.);
    int l = m + (N+1)%2;
    return A(m,l);
}

void invA(int A,int &k, int &l) {
    double m = floor(sqrt(A));
    double M = ceil(sqrt(A));
    double G;
    if (M == m) {
        k = m;
        l = m;
        return;
    }
    if (A < (M*M + m*m)/2.) {
         G = 2*m+1;
    }
    else {
         G = 2*M;
    }
    k = ceil(G*G/4. - G/2.) + G - A;
    l = G-k;
    return;
}

int sigma(int a) {
    if (a%2)
        return -1;
    return 1;
}
int kronecker(int a, int b) {
    if (a == b)
        return 1;
    return 0;
}

double onebody(int n, int m, double hw){
    return hw*(2*n + abs(m) + 1);
}

int m(int a) {
    return - abs(E(a)-1) + 2*floor((a-1 - E(a)*(E(a)-1))/2.);
}

int n(int a) {
    return 0.5*(E(a) - abs(m(a)) - 1 );
}

int E(int a) {
    return ceil(0.5*sqrt(1+4*a) - 0.5);
}

int merge(int alpha, int beta, int N) {
    return alpha*(2*N - alpha - 1)/2. + beta;
}
