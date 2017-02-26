#ifndef BASIS_H
#define BASIS_H

#include <iostream>
using namespace std;

class basis
{
public:
    basis(int n);
    int * get_state(int n);
    int get_index(int * config);
    void print();
    void to_file(string filename);
    double hermite(int n, double x);
    double psi_n(int n, double x);
    double psi(int nx,int ny, double x, double y);

private:
    int ** state;
    int numStates;
};

#endif // BASIS_H
