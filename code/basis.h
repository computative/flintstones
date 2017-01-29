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
    int print();
    int to_file(string filename);
    double hermite(int n, double * x);

private:
    int ** state;
    int numStates;
};

#endif // BASIS_H
