#ifndef GC_H
#define GC_H

#include <iostream>
using namespace std;

class gc
{
public:
    gc(string type, int n);
    double quad4(double (*f)(double , double , double , double ));

private:
    double hermite4(double (*f)(double , double , double , double ));
    string type;
    int n;
};

#endif // GC_H
