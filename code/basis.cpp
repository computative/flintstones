#include <iostream>
#include <fstream>
#include "basis.h"
#include "maths.h"
#include "math.h"
#include "adapt.h"

using namespace std;

basis::basis(int n) // number of shells is n
{
    numStates = n*(n+1);                                // create table of numStates
    state = new int * [numStates];
    for (int i = 0; i < n; i++) {                       // for each shell
        for (int j = 0; j < 2*(i+1); j++) {             // for each state in shell
            state[i*(i+1) + j] = new int[4];
            state[i*(i+1) + j][0] = (j - j%2)/2;        // set nx
            state[i*(i+1) + j][1] = i - (j - j%2)/2;    // set ny
            state[i*(i+1) + j][2] = j % 2;              // set spinprojection
            state[i*(i+1) + j][3] = i + 1;              // set energy/omega
        }
    }
}

int basis::get_index(int * config)
{
    int i = config[3] - 1;              // get shell - 1
    int nx = config[0];                 // get nx
    int sigma = config[2];              // get spinprojection
    return i*(i+1) + 2*nx + sigma+1;    // return statenumber = 1,2,3,...
}

int * basis::get_state(int ni)
{
    return state[ni-1];
}

int basis::print(){
    int E = state[0][3];                // the lowest energy level
    for(int i = 0; i < numStates; i++){ // for each state
        if (state[i][3] != E) {         // if higher energy level
            cout << endl;               // then print linebreak to table
            E = state[i][3];
        }
        cout << '(' << state[i][0] <<  ',' << state[i][1] <<  ',' << state[i][2] <<  ',' << state[i][3] << "),";
    }
    cout << endl << "Legend: (nx,ny,sigma,E). sigma = 1 => spin = +1/2" << endl;
}

int basis::to_file(string filename){
    int E = state[0][3];
    ofstream myfile;
    myfile.open (filename);
    for(int i = 0; i < numStates; i++){
        if (state[i][3] != E) {         // if higher energy level
            myfile << endl;             //
            E = state[i][3];            // then print linebreak to table
        }
        myfile << '(' << state[i][0] <<  ',' << state[i][1] <<  ',' << state[i][2] <<  ',' << state[i][3] << "),";
    }
    myfile.close();
}

double basis::hermite(int n, double * x)
{
    Adapt s(1e-16);
    ans = s.integrate(func,a,2*M_PI);
}
