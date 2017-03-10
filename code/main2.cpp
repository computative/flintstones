#include <iostream>
#include <fstream>
#include <armadillo>
#include "basis.h"
#include "gc.h"
#include "math.h"
#include "Coulomb_Functions.hpp"

using namespace std;
using namespace arma;


double hw = 1;
int m(int a);
int n(int a);
int E(int a);
int A(int m, int l);
int B(int N);
void invA(int A,int &k, int &l);
int sigma(int a);
int kronecker(int a, int b);
double onebody(int n, int m);

int main(int, char *[])
{

    int F = 6; // fermi level
    // truncation limit
    int R = 4; // num shells
    int N = R*(R+1); // num of states less than R

    vec singleparticleH = zeros<vec>(N);
    for (int i = 0; i<N; i++) {
        singleparticleH[i] = onebody(n(i+1), m(i+1));
    }
    // inital C- and rho-matrix
    mat C = eye<mat>(N,N);
    mat rho = zeros<mat>(N,N);
    // density matrix, konvensjon som i pythonkode
    for (int gamma = 0; gamma <N; gamma++) {
        for (int delta = 0; delta <N; delta++) {
            for (int a = 0; a <F; a++)
                rho(gamma,delta) += C(gamma,a)*C(delta,a);
        }
    }
    int maxHFiter = 100;
    double epsilon =  1.0e-5;
    double difference = 1.0;
    int hf_count = 0;
    vec oldenergies = zeros<vec>(N);
    vec newenergies = zeros<vec>(N);
    vec spenergies;

    mat nninteraction = zeros<mat>(N*N+1,N*N+1);
    for (int k = 1; k <= N*N; k++) {
        for (int l = 1; l <= N*N; l++) {
            int alpha, beta, gamma, delta;
            invA(k,alpha,beta);
            invA(l,gamma,delta);
            int n1 = n(alpha);
            int ml1 = m(alpha);
            int n2 = n(beta);
            int ml2 = m(beta);
            int n3 = n(gamma);
            int ml3 = m(gamma);
            int n4 = n(delta);
            int ml4 = m(delta);
            nninteraction(k,l) = Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4);
                               //- kronecker(sigma(alpha), sigma(delta))*kronecker(sigma(beta), sigma(gamma))*Coulomb_HO(hw, n1, ml1, n2, ml2, n4, ml4, n3, ml3);
        }
    }

    /*
    double **** nninteraction = new double***[N];

    for (int alpha = 0; alpha < N; alpha++) {
        nninteraction[alpha] = new double**[N];
        for (int beta = 0; beta < N; beta++ ) {
            nninteraction[alpha][beta] = new double*[N];
            for (int gamma = 0; gamma < N; gamma++) {
                nninteraction[alpha][beta][gamma] = new double[N];
                for (int delta = 0; delta < N; delta ++){
                    int n1 = n(alpha+1);
                    int ml1 = m(alpha+1);
                    int n2 = n(beta+1);
                    int ml2 = m(beta+1);
                    int n3 = n(gamma+1);
                    int ml3 = m(gamma+1);
                    int n4 = n(delta+1);
                    int ml4 = m(delta+1);
                    nninteraction[alpha][beta][gamma][delta] = 0;
                    if (kronecker(sigma(alpha+1), sigma(gamma+1))*kronecker(sigma(beta+1), sigma(delta+1)) and ml1 + ml2 == ml3 + ml4)
                        nninteraction[alpha][beta][gamma][delta] += Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4);
                    if (kronecker(sigma(alpha+1), sigma(delta+1))*kronecker(sigma(beta+1), sigma(gamma+1)) and ml1 + ml2 == ml3 + ml4)
                        nninteraction[alpha][beta][gamma][delta] += -Coulomb_HO(hw, n1, ml1, n2, ml2, n4, ml4, n3, ml3);
                }
            }
        }
    }
    */
    while (hf_count < maxHFiter and difference > epsilon) {
        mat HFmatrix = zeros<mat>(N,N);
        for (int alpha = 0; alpha < N; alpha++) {
            for (int beta = 0; beta < N; beta++ ) {
                // If tests for three-dimensional systems, including isospin conservation
                // if l[alpha] != l[beta] and j[alpha] != j[beta] and mj[alpha] != mj[beta] and tz[alpha] != tz[beta]: continue
                // Setting up the Fock matrix using the density matrix and antisymmetrized NN interaction in m-scheme """
                //HFmatrix(alpha,beta) = 0.0;
                for (int gamma = 0; gamma < N; gamma++) {
                    for (int delta = 0; delta < N; delta ++){
                        // continue hvis integral er null
                        int n1 = n(alpha+1);
                        int ml1 = m(alpha+1);
                        int n2 = n(gamma+1);
                        int ml2 = m(gamma+1);
                        int n3 = n(beta+1);
                        int ml3 = m(beta+1);
                        int n4 = n(delta+1);
                        int ml4 = m(delta+1);
                        /*
                        if (kronecker(sigma(alpha+1), sigma(beta+1))*kronecker(sigma(gamma+1), sigma(delta+1)) and ml1 + ml2 == ml3 + ml4) {
                            HFmatrix(alpha,beta) += rho(gamma,delta)*Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4);
                        }
                        if (kronecker(sigma(alpha+1), sigma(delta+1))*kronecker(sigma(gamma+1), sigma(beta+1)) and ml1 + ml2 == ml3 + ml4) {
                            HFmatrix(alpha,beta) += -rho(gamma,delta)*Coulomb_HO(hw, n1, ml1, n2, ml2, n4, ml4, n3, ml3);
                        }
                        */
                        //HFmatrix(alpha,beta) += rho(gamma,delta)*nninteraction[alpha][gamma][beta][delta];
                        HFmatrix(alpha,beta) += rho(gamma,delta)*(
                                    kronecker(sigma(alpha+1), sigma(beta+1))*kronecker(sigma(gamma+1), sigma(delta+1))*nninteraction(A(alpha+1,gamma+1) , A(beta+1,delta+1))
                                  - kronecker(sigma(alpha+1), sigma(delta+1))*kronecker(sigma(gamma+1), sigma(beta+1))*nninteraction(A(alpha+1,gamma+1) , A(delta+1,beta+1))
                                    );
                    }
                }
            }
            //  Adding the one-body term, here plain harmonic oscillator
            HFmatrix(alpha,alpha) += singleparticleH[alpha];
        }
        eig_sym( spenergies, C, HFmatrix);
        rho = zeros<mat>(N,N);
        for (int gamma = 0; gamma <N; gamma++) {
            // density matrix, konvensjon som i pythonkode
            for (int delta = 0; delta <N; delta++) {
                //rho(gamma,delta) = 0;
                for (int a=0; a <F; a++)
                    rho(gamma,delta) += C(gamma,a)*C(delta,a);
            }
        }

        newenergies = spenergies;
        // Brute force computation of difference between previous and new sp HF energies
        difference = 0.0;
        for (int i = 0; i<N; i++) {
            difference += abs(newenergies[i]-oldenergies[i])/N;
        }
        oldenergies = newenergies;
        cout << "Single-particle energies, ordering may have changed" << endl;
        for (int i = 0; i <N; i++ ) {
            printf("%i %g\n", i, oldenergies[i]);
        }
        cout << difference << endl;
        hf_count += 1;
    }

    // get energy
    double EHF = 0;
    for (int alpha = 0; alpha < F; alpha++)
        EHF += oldenergies[alpha];

    for (int alpha = 0; alpha < N; alpha++) {
        for (int beta = 0; beta < N; beta++ ) {
            for (int gamma = 0; gamma < N; gamma++) {
                for (int delta = 0; delta < N; delta ++){
                    int n1 = n(alpha+1);
                    int ml1 = m(alpha+1);
                    int n2 = n(beta+1);
                    int ml2 = m(beta+1);
                    int n3 = n(gamma+1);
                    int ml3 = m(gamma+1);
                    int n4 = n(delta+1);
                    int ml4 = m(delta+1);
                    if (kronecker(sigma(alpha+1), sigma(gamma+1))*kronecker(sigma(beta+1), sigma(delta+1)))
                    EHF -= 0.5*rho(alpha,gamma)*rho(beta,delta)*Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4);
                    if (kronecker(sigma(alpha+1), sigma(delta+1))*kronecker(sigma(gamma+1), sigma(beta+1)))
                    EHF -= -0.5*rho(alpha,gamma)*rho(beta,delta)*Coulomb_HO(hw, n1, ml1, n2, ml2, n4, ml4, n3, ml3);
                }
            }
        }
    }
    cout << EHF << endl;

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

double onebody(int n, int m){
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
