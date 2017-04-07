//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_main.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 9/19/16.
//
//

#include "MC_main.hpp"

#include "MC_pathData.hpp"
#include "MC_Cycle.hpp"

// Arguments of the main function:
// 1. t = time
// 2. nRuns = number of Monte Carlo runs we start
// 3. nFwdPaths = number of forward paths
// 4. calcSpinQ = whether to calculate spin correlations or not
// 5. numberOfInitialSpinStates = number of spin states we average over if we run 'MC_hole_infiniteTSpinStateAvg()'.
int main(int argc,char *argv[])
{
    // Define variables
    double t;
    int numberOfRuns;
    int nPaths;
    int numberOfInitialSpinStates;
    bool calcSpinQ;
    vector<double> transitionProbs;
    
    // Check if there is enough arguments
    if(!(argc==4 || argc==5 || argc==6)) {
        cerr << "Error in main(): Incorrect number of arguments." << endl;
        return 1;
    }
    
    // Get inputs
    t = atof(argv[1]);
    numberOfRuns = atol(argv[2]);
    nPaths = atol(argv[3]);
    if(argc>=5) {
        calcSpinQ = atol(argv[4]);
    } else {
        calcSpinQ = true;
    }
    if(argc==6) {
        numberOfInitialSpinStates = atol(argv[5]);
        
    }
    
    cout << "Parameters: " << endl;
    cout << "time = " << t << endl;
    cout << "number of runs = " << numberOfRuns << endl;
    cout << "number of paths generated = " << nPaths << endl;
    cout << "calculate spin correlations: " << (calcSpinQ ? "Yes" : "No");
    if(argc==6) {
        cout << "number of random initial spin states = " << numberOfInitialSpinStates << endl;
    }
    cout << endl << endl << endl;
    
    // Run simulation
    transitionProbs = MC_hole_infiniteT(t, numberOfRuns, nPaths, calcSpinQ);
    transitionProbs = MC_hole_infiniteT_TraceFwdBwd(t, numberOfRuns, nPaths, calcSpinQ);
    // transitionProbs = MC_hole_infiniteTSpinStateAvg(t, numberOfRuns, nPaths, numberOfInitialSpinStates, calcSpinQ);
    // transitionProbs = MC_hole_infiniteTSpinStateAvgFwdBwd(t, numberOfRuns, nPaths, numberOfInitialSpinStates, calcSpinQ);
    // transitionProbs = MC_hole_infiniteT_infiniteSpin(t, numberOfRuns, nPaths);
    // transitionProbs = MC_hole_infiniteT_infiniteSpinFwdBwd(t, numberOfRuns, nPaths);
    // transitionProbs = MC_hole_Neel(t, numberOfRuns, nPaths, calcSpinQ);
    // transitionProbs = MC_hole_Ferromagnet(t, numberOfRuns, nPaths, calcSpinQ);
    
    
    // Print results
    for(int i=0; i<4*LATTICE_SIZE*LATTICE_SIZE; i++) {
        cout << "P(" << indexToX(i) << ", " << indexToY(i) << ") = " << transitionProbs[i] << endl;
    }
    
    return 0;
}
