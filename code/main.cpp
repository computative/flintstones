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
int sigma(int a);
int kronecker(int a, int b);
double onebody(int n, int m);

int main(int, char *[])
{
    //cout << kronecker(1,2) << " " <<sigma(1) <<" " << sigma(7) <<" " << kronecker(1,1);
    //return 0;
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
                        HFmatrix(alpha,beta) += rho(gamma,delta)*(
                                    kronecker(sigma(alpha+1), sigma(beta+1))*kronecker(sigma(gamma+1), sigma(delta+1))*Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4)
                                    - kronecker(sigma(alpha+1), sigma(delta+1))*kronecker(sigma(gamma+1), sigma(beta+1))*Coulomb_HO(hw, n1, ml1, n2, ml2, n4, ml4, n3, ml3) );
                    }
                }
            }
                //  Adding the one-body term, here plain harmonic oscillator
            HFmatrix(alpha,alpha) += singleparticleH[alpha];
        }
        //HFmatrix.print();
        eig_sym( spenergies, C, HFmatrix);

        //mat rhoOld = rho;
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
    return 0;
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
