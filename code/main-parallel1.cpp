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
    map < int , map < int, map < int, map< int, double > > > > nninteraction;
    int F = 6; // fermi level
    // truncation limit
    int R = 6; // num shells
    int N = R*(R+1); // num of states less than R
    //cout << R << endl;
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
    mat jobs = zeros<mat>(3,4*(R-1) + 1);
    /*
    jobs(0,0) = 1;
    jobs(0,1) = 1;
    jobs.print();
    int Ms;
    int M;
    cout << all(jobs) << endl;
    while ( accu(jobs) < 3*(4*(R-1)+1) ) {
        for (Ms = -2; Ms <= 2; Ms += 2) {
            for (M = -2*(R-1); M <= 2*(R-1); M++ ) {
                if (!jobs(Ms/2 + 1,M+2*(R-1)))
                    goto theEnd;
            }
        }
        theEnd:
        cout << Ms << " " << M << endl;
        jobs(Ms/2 + 1,M+2*(R-1)) = 1;
    }
    return 0;
    */

    #pragma omp parallel num_threads(8)
    {
        int Ms;
        int M;
        bool bit = false;
        int id = omp_get_thread_num();
        int threads = omp_get_num_threads();
        while (accu(jobs) < 3*(4*(R-1)+1)) {
            #pragma omp critical
            {
                if (accu(jobs) >= 3*(4*(R-1)+1))
                    bit = true;
                for (Ms = -2; Ms <= 2; Ms += 2) {
                    for (M = -2*(R-1); M <= 2*(R-1); M++ ) {
                        if (!jobs(Ms/2 + 1,M+2*(R-1)))
                            goto theEnd;
                    }
                }
                theEnd:;
                cout << Ms << " " << M << endl;
                jobs(Ms/2 + 1,M+2*(R-1)) = 1;
            }
            if(bit)
                break;
            for (int alpha = 0; alpha < N; alpha++) {
                //nninteraction[alpha] = new double**[N];
                for (int beta = 0; beta < N; beta++ ) {
                    //nninteraction[alpha][beta] = new double*[N];
                    for (int gamma = 0; gamma < N; gamma++) {
                        //nninteraction[alpha][beta][gamma] = new double[N];
                        for (int delta = 0; delta < N; delta ++){
                            int ml1 = m(alpha+1);
                            int ml2 = m(beta+1);
                            int ml3 = m(gamma+1);
                            int ml4 = m(delta+1);
                            if(not (M == ml1 + ml2 and Ms == sigma(alpha+1) + sigma(beta+1)) )
                                continue;
                            if (sigma(alpha+1) + sigma(beta+1) == sigma(gamma+1) + sigma(delta+1)  and ml1 + ml2 == ml3 + ml4) {
                                if(exists(nninteraction,beta,alpha,delta,gamma))
                                    continue;
                                if(exists(nninteraction,alpha,beta,delta,gamma))
                                    continue;
                                if(exists(nninteraction,beta,alpha,gamma,delta))
                                    continue;
                                if(exists(nninteraction,gamma,delta,alpha,beta))
                                    continue;
                                int n1 = n(alpha+1);
                                int n2 = n(beta+1);
                                int n3 = n(gamma+1);
                                int n4 = n(delta+1);
                                nninteraction[alpha][beta][gamma][delta] = kronecker(sigma(alpha+1), sigma(gamma+1))*kronecker(sigma(beta+1), sigma(delta+1))*Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4)
                                        -kronecker(sigma(alpha+1), sigma(delta+1))*kronecker(sigma(beta+1), sigma(gamma+1))*Coulomb_HO(hw, n1, ml1, n2, ml2, n4, ml4, n3, ml3);
                            }
                        }
                    }
                }
            }
        }
    }
    while (hf_count < maxHFiter and difference > epsilon) {
        mat HFmatrix = zeros<mat>(N,N);
        for (int alpha = 0; alpha < N; alpha++) {
            for (int beta = 0; beta < N; beta++ ) {
                for (int gamma = 0; gamma < N; gamma++) {
                    for (int delta = 0; delta < N; delta ++){
                        int ml1 = m(alpha+1);
                        int ml2 = m(gamma+1);
                        int ml3 = m(beta+1);
                        int ml4 = m(delta+1);
                        if (sigma(alpha+1) + sigma(gamma+1) == sigma(beta+1) + sigma(delta+1)  and ml1 + ml2 == ml3 + ml4) {
                            if(exists(nninteraction,gamma,alpha,delta,beta)) {
                                HFmatrix(alpha,beta) += rho(gamma,delta)*nninteraction[gamma][alpha][delta][beta];
                                continue;
                            } else if(exists(nninteraction,alpha,gamma,delta,beta)) {
                                HFmatrix(alpha,beta) += -rho(gamma,delta)*nninteraction[alpha][gamma][delta][beta];
                                continue;
                            } else if(exists(nninteraction,gamma,alpha,beta,delta)) {
                                HFmatrix(alpha,beta) += -rho(gamma,delta)*nninteraction[gamma][alpha][beta][delta];
                                continue;
                            } else if(exists(nninteraction,beta,delta,alpha,gamma)) {
                                HFmatrix(alpha,beta) += rho(gamma,delta)*nninteraction[beta][delta][alpha][gamma];
                                continue;
                            }
                            HFmatrix(alpha,beta) += rho(gamma,delta)*nninteraction[alpha][gamma][beta][delta];
                        }
                    }
                }
            }
            HFmatrix(alpha,alpha) += singleparticleH[alpha];
        }
        eig_sym( spenergies, C, HFmatrix);
        rho = zeros<mat>(N,N);
        for (int gamma = 0; gamma <N; gamma++) {
            for (int delta = 0; delta <N; delta++) {
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
                    int ml1 = m(alpha+1);
                    int ml2 = m(beta+1);
                    int ml3 = m(gamma+1);
                    int ml4 = m(delta+1);
                    if (sigma(alpha+1) + sigma(beta+1) == sigma(gamma+1) + sigma(delta+1)  and ml1 + ml2 == ml3 + ml4)
                        if(exists(nninteraction,beta,alpha,delta,gamma)) {
                            EHF -= 0.5*rho(alpha,gamma)*rho(beta,delta)*nninteraction[beta][alpha][delta][gamma];
                            continue;
                        } else if(exists(nninteraction,alpha,beta,delta,gamma)) {
                            EHF -= -0.5*rho(alpha,gamma)*rho(beta,delta)*nninteraction[alpha][beta][delta][gamma];
                            continue;
                        } else if(exists(nninteraction,beta,alpha,gamma,delta)) {
                            EHF -= -0.5*rho(alpha,gamma)*rho(beta,delta)*nninteraction[beta][alpha][gamma][delta];
                            continue;
                        } else if(exists(nninteraction,gamma,delta,alpha,beta)) {
                            EHF -= 0.5*rho(alpha,gamma)*rho(beta,delta)*nninteraction[gamma][delta][alpha][beta];
                            continue;
                        }
                        EHF -= 0.5*rho(alpha,gamma)*rho(beta,delta)*nninteraction[alpha][beta][gamma][delta];
                }
            }
        }
    }
    cout <<  setprecision(9) << EHF << endl;
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
