#include <iostream>
#include <fstream>
#include <armadillo>
#include "basis.h"
#include "gc.h"
#include "math.h"
#include "maths.h"
#include "Coulomb_Functions.hpp"

using namespace std;
using namespace arma;


int * ns; //global variabel fordi jeg ønsker det
double hw = 1;
basis B = basis(4);
gc integrator = gc("hermite", 20);
double energy(int p,int q,int r, int s);
double f(double x1,double y1,double x2,double y2);
double g(double x1,double y1,double x2,double y2);
int m(int a);
int n(int a);
int E(int a);
double onebody(int n, int m);

int main(int, char *[])
{
    int F = 2; // fermi level
    // truncation limit
    int R = 3; // num shells
    int N = R*(R+1); // num of states less than R

    // inital C- and rho-matrix
    mat C = eye<mat>(N,N);
    mat rho = zeros<mat>(N,N);
    // density matrix, konvensjon som i pythonkode
    for (int gamma = 0; gamma <N; gamma++) {
        for (int delta = 0; delta <N; delta++) {
            for (int a; a <F; a++)
                rho(gamma,delta) = C(gamma,a)*C(delta,a);
        }
    }
    // fock matrix
    mat HF = zeros<mat>(N,N);
    for (int gamma = 0; gamma <N; gamma++) {
        for (int delta = 0; delta <N; delta++) {
            for (int a = 0; a <F; a++)
                rho(gamma,delta) = C(gamma,a)*C(delta,a);
        }
    }
    //cout << energy(1,1,1,1) << endl;
    int maxHFiter = 100;
    double epsilon =  1.0e-5;
    double difference = 1.0;
    int hf_count = 0;
    int sumFockTerm;
    vec oldenergies = zeros<vec>(N);
    vec newenergies = zeros<vec>(N);
    vec singleparticleH = zeros<vec>(N);
    for (int i = 0; i<N; i++) {
            singleparticleH[i] = onebody(n(i+1), m(i+1));
    }
    vec spenergies;
    while (hf_count < maxHFiter and difference > epsilon) {
        mat HFmatrix = zeros<mat>(N,N);
        for (int alpha = 0; alpha < N; alpha++) {
            for (int beta = 0; beta < N; beta++ ) {
                // If tests for three-dimensional systems, including isospin conservation
                // if l[alpha] != l[beta] and j[alpha] != j[beta] and mj[alpha] != mj[beta] and tz[alpha] != tz[beta]: continue
                // Setting up the Fock matrix using the density matrix and antisymmetrized NN interaction in m-scheme """
                sumFockTerm = 0.0;
                for (int gamma = 0; gamma < N; gamma++) {
                    for (int delta = 0; delta < N; delta ++){
                        // continue hvis integral er null
                        int n1 = n(alpha+1);
                        int ml1 = m(alpha+1);
                        int n2 = n(beta+1);
                        int ml2 = m(beta+1);
                        int n3 = n(gamma+1);
                        int ml3 = m(gamma+1);
                        int n4 = n(delta+1);
                        int ml4 = m(delta+1);
                        cout << rho(gamma,delta)*(Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4) -Coulomb_HO(hw, n1, ml1, n2, ml2, n4, ml4, n3, ml3) ) << endl;
                        sumFockTerm += rho(gamma,delta)*(
                                    Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4)
                                  - Coulomb_HO(hw, n1, ml1, n2, ml2, n4, ml4, n3, ml3) );
                    HFmatrix(alpha,beta) = sumFockTerm;
                    //  Adding the one-body term, here plain harmonic oscillator
                    if (beta == alpha) HFmatrix(alpha,alpha) += singleparticleH[alpha];
                    }
                }
            }
        }
        eig_sym( spenergies, C, HFmatrix);
        //HFmatrix.print();
        C.print();
        // density matrix, konvensjon som i pythonkode
        for (int gamma = 0; gamma <N; gamma++) {
            for (int delta = 0; delta <N; delta++) {
                rho(gamma,delta) = 0;
                for (int a=0; a <F; a++)
                    rho(gamma,delta) += C(gamma,a)*C(delta,a);
            }
        }
        rho.print();
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

double f(double x1,double y1,double x2,double y2) {
    if (x1 != x2 && y1 != y2) {
        return B.hermite(1,x1)*B.hermite(2,y1)*B.hermite(1,x2)*B.hermite(1,y2)*B.hermite(1,x1)*B.hermite(2,y1)*B.hermite(1,x2)*B.hermite(1,y2)/sqrt(((x1)-(x2))*((x1)-(x2)) + ((y1)-(y2))*((y1)-(y2)));
    }
    return 0.0;
}

double g(double x1,double y1,double x2,double y2) {
    if (x1 != x2 && y1 != y2) {
        return B.hermite(ns[0],x1)*B.hermite(ns[1],y1)*B.hermite(ns[2],x2)*B.hermite(ns[3],y2)*B.hermite(ns[4],x1)*B.hermite(ns[5],y1)*B.hermite(ns[6],x2)*B.hermite(ns[7],y2)/sqrt(((x1)-(x2))*((x1)-(x2)) + ((y1)-(y2))*((y1)-(y2)));
    }
    return 0.0;
}

double energy(int p,int q,int r, int s) {
    ns = new int[8];
    double omega = 1;
    double e = 1;
    double epsilon0 = 1;
    ns[0] = B.get_state(p)[0]; ns[1] = B.get_state(p)[1];
    ns[2] = B.get_state(q)[0]; ns[3] = B.get_state(q)[1];
    ns[4] = B.get_state(r)[0]; ns[5] = B.get_state(r)[1];
    ns[6] = B.get_state(s)[0]; ns[7] = B.get_state(s)[1];
    double C = 1;
    for (int i = 0; i< 8; i++) {
        C *= sqrt(pow(2, ns[i])*factorial(ns[i]));
    }
    C = (sqrt(omega)*e*e/(4*M_PI*M_PI*M_PI*epsilon0))/C;
    C = -1;
    return -C*integrator.quad4(&g);
}
