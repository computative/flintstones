#include <iostream>
#include "basis.h"

using namespace std;

int main(int, char *[])
{
    basis B = basis(4);
    B.print();
    B.to_file("../resources/output.txt");
    // B.psi(nx,ny,x)
    return 0;
}
