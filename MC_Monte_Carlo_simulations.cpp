//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_Monte_Carlo_simulations.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/18/16.
//
//

#include "MC_Monte_Carlo_simulations.hpp"



// Propagation of the hole in an infinite temperature spin environment
//
// Monte Carlo calculation, calling holePropagation::run() nThreads times, and then averaging
// the transition probabilities and normalizing them to 1. We also perform a symmetry averaging
// of the transition probabilities, so that they reflect the symmetry properties of the lattice.
vector<double> MC_hole_infiniteT(double t, int numberOfRuns, int nFwdPaths, bool calcSpinQ, double accuracy) {
    // Define the class that will run the simulation and store our data
    holePropagation simulation(t, nFwdPaths, calcSpinQ);
    
    // Simulation results
    vector<double> transProbs, transProbsNew;
    vector< vector<double> > labFrameSpinCorr, labFrameSpinCorrNew;
    vector< vector<double> > holeFrameSpinCorr, holeFrameSpinCorrNew;
    
    // Errors w.r.t the last run
    double transProbsErr;
    double labFrameSpinCorrErr;
    double holeFrameSpinCorrErr;
    
    // Load data of previous runs, if any
    if(!simulation.load()) {
        simulation.saveAllFilenames();
        simulation.save_returnProb_and_RMSdistance();
    }
    
    // Run simulation 'numberOfRuns' times, and save partial results after each run
    while(simulation.getNumberOfRuns() < numberOfRuns) {
        simulation.run();
        simulation.backup();
        simulation.save();
        
        // Check errors and stop if they are smaller than 'accuracy'
        transProbsNew = simulation.getTransitionProbabilities();
        labFrameSpinCorrNew = simulation.getLabFrameSpinCorrelationFn();
        holeFrameSpinCorrNew = simulation.getHoleFrameSpinCorrelationFn();
        
        transProbsErr = getTransProbsErr(transProbs, transProbsNew);
        labFrameSpinCorrErr = getSpinCorrErr(labFrameSpinCorr, labFrameSpinCorrNew);
        holeFrameSpinCorrErr = getSpinCorrErr(holeFrameSpinCorr, holeFrameSpinCorrNew);
        
        printErrors(transProbsErr, labFrameSpinCorrErr, holeFrameSpinCorrErr);
        if(transProbsErr<accuracy && labFrameSpinCorrErr<accuracy && holeFrameSpinCorrErr<accuracy) {
            break;
        }
        
        transProbs = transProbsNew;
        labFrameSpinCorr = labFrameSpinCorrNew;
        holeFrameSpinCorr = holeFrameSpinCorrNew;
    }
    
    return simulation.getTransitionProbabilities();
}



// Propagation of the hole in an infinite temperature spin environment
//
// Same as MC_hole_infiniteT(), except that it takes a different set for forward and backward propagating paths.
// This doubles the runtime and the memory requirement, but it also means that the code can run for times t>2.
// At these times, MC_hole_infiniteT gave unreliable results, even for 200 million paths.
vector<double> MC_hole_infiniteT_TraceFwdBwd(double t, int numberOfRuns, int nFwdPaths, bool calcSpinQ, double accuracy) {
    // Define the class that will run the simulation and store our data
    holePropagationTraceFwdBwd simulation(t, nFwdPaths, calcSpinQ);
    
    // Simulation results
    vector<double> transProbs, transProbsNew;
    vector< vector<double> > labFrameSpinCorr, labFrameSpinCorrNew;
    vector< vector<double> > holeFrameSpinCorr, holeFrameSpinCorrNew;
    
    // Errors w.r.t the last run
    double transProbsErr;
    double labFrameSpinCorrErr;
    double holeFrameSpinCorrErr;
    
    // Load data of previous runs, if any
    if(!simulation.load()) {
        simulation.saveAllFilenames();
        simulation.save_returnProb_and_RMSdistance();
    }
    
    // Run simulation 'numberOfRuns' times, and save partial results after each run
    while(simulation.getNumberOfRuns() < numberOfRuns) {
        simulation.run();
        simulation.backup();
        simulation.save();
        
        // Check errors and stop if they are smaller than 'accuracy'
        transProbsNew = simulation.getTransitionProbabilities();
        labFrameSpinCorrNew = simulation.getLabFrameSpinCorrelationFn();
        holeFrameSpinCorrNew = simulation.getHoleFrameSpinCorrelationFn();
        
        transProbsErr = getTransProbsErr(transProbs, transProbsNew);
        labFrameSpinCorrErr = getSpinCorrErr(labFrameSpinCorr, labFrameSpinCorrNew);
        holeFrameSpinCorrErr = getSpinCorrErr(holeFrameSpinCorr, holeFrameSpinCorrNew);
        
        printErrors(transProbsErr, labFrameSpinCorrErr, holeFrameSpinCorrErr);
        if(transProbsErr<accuracy && labFrameSpinCorrErr<accuracy && holeFrameSpinCorrErr<accuracy) {
            break;
        }
        
        transProbs = transProbsNew;
        labFrameSpinCorr = labFrameSpinCorrNew;
        holeFrameSpinCorr = holeFrameSpinCorrNew;
    }
    
    return simulation.getTransitionProbabilities();
}



// Propagation of the hole in an infinite temperature spin environment
//
// Monte Carlo calculation, calling holePropagationSpinStateAvg::run() nThreads times.
// Instead of performing a trace over spin configurations exactly, it samples random initial
// spin configurations 'numberOfInitialSpinStates' times, and averages the quantities we are
// interested in over this spin configuration space
vector<double> MC_hole_infiniteTSpinStateAvg(double t, int numberOfRuns, int nFwdPaths, int numberOfInitialSpinStates, bool calcSpinQ, double accuracy) {
    // Define the class that will run the simulation and store our data
    holePropagationSpinStateAvg simulation(t, nFwdPaths, numberOfInitialSpinStates, calcSpinQ);
    
    // Simulation results
    vector<double> transProbs, transProbsNew;
    vector< vector<double> > labFrameSpinCorr, labFrameSpinCorrNew;
    vector< vector<double> > holeFrameSpinCorr, holeFrameSpinCorrNew;
    
    // Errors w.r.t the last run
    double transProbsErr;
    double labFrameSpinCorrErr;
    double holeFrameSpinCorrErr;

    // Load data of previous runs, if any
    if(!simulation.load()) {
        simulation.saveAllFilenames();
        simulation.save_returnProb_and_RMSdistance();
    }
    
    // Run simulation 'numberOfRuns' times, and save partial results after each run
    while(simulation.getNumberOfRuns() < numberOfRuns) {
        simulation.run();
        simulation.backup();
        simulation.save();
        
        // Check errors and stop if they are smaller than 'accuracy'
        transProbsNew = simulation.getTransitionProbabilities();
        labFrameSpinCorrNew = simulation.getLabFrameSpinCorrelationFn();
        holeFrameSpinCorrNew = simulation.getHoleFrameSpinCorrelationFn();
        
        transProbsErr = getTransProbsErr(transProbs, transProbsNew);
        labFrameSpinCorrErr = getSpinCorrErr(labFrameSpinCorr, labFrameSpinCorrNew);
        holeFrameSpinCorrErr = getSpinCorrErr(holeFrameSpinCorr, holeFrameSpinCorrNew);
        
        printErrors(transProbsErr, labFrameSpinCorrErr, holeFrameSpinCorrErr);
        if(transProbsErr<accuracy && labFrameSpinCorrErr<accuracy && holeFrameSpinCorrErr<accuracy) {
            break;
        }
        
        transProbs = transProbsNew;
        labFrameSpinCorr = labFrameSpinCorrNew;
        holeFrameSpinCorr = holeFrameSpinCorrNew;
    }
    cout << endl;
    
    return simulation.getTransitionProbabilities();
}


// Monte Carlo calculation, calling 'holePropagationSpinStateAvgFwdBwd::run()' nThreads times.
// Almost identical to 'MC_hole_infiniteTSpinStateAvg', the only difference is that we create a separate
// set of forward and backward paths, to make sure that the interference effects are taken into account appropriately.
//
// IMPORTANT: we took something for granted during this implementation in holePropagationSpinStateAvgFwdBwd,
// to avoid having to write a large portion of new code:
// When we take the overlap of quantum states fwdState and bwdState (which are obtained from the initial state by
// transforming it using a forward and a backward path), we should update the transition probabilities and spin
// correlations by <bwdState|fwdState>. However, we update them instead only by the real part of this expression
// real(<bwdState|fwdState>). This is correct in the limit of infinitely many paths, since the imaginary
// parts should cancel out. However, this is not necessarly true in the case of finitely many paths.
// The cancelation of imaginary parts would thus be a useful sanity check.
// However, implementing this would be a huge hassle, and would require me to write a completely new base class
// instead of 'holePropagation', that could handle complex bins for the transition probabilities. Since probably
// it's not worth the time, I staid with taking the real part of the overlaps above.
vector<double> MC_hole_infiniteTSpinStateAvgFwdBwd(double t, int numberOfRuns, int nFwdPaths, int numberOfInitialSpinStates, bool calcSpinQ, double accuracy) {
    // Define the class that will run the simulation and store our data
    holePropagationSpinStateAvgFwdBwd simulation(t, nFwdPaths, numberOfInitialSpinStates, calcSpinQ);
    
    // Simulation results
    vector<double> transProbs, transProbsNew;
    vector< vector<double> > labFrameSpinCorr, labFrameSpinCorrNew;
    vector< vector<double> > holeFrameSpinCorr, holeFrameSpinCorrNew;
    
    // Errors w.r.t the last run
    double transProbsErr;
    double labFrameSpinCorrErr;
    double holeFrameSpinCorrErr;
    
    // Load data of previous runs, if any
    if(!simulation.load()) {
        simulation.saveAllFilenames();
        simulation.save_returnProb_and_RMSdistance();
    }
    
    // Run simulation 'numberOfRuns' times, and save partial results after each run
    while(simulation.getNumberOfRuns() < numberOfRuns) {
        simulation.run();
        simulation.backup();
        simulation.save();
        
        // Check errors and stop if they are smaller than 'accuracy'
        transProbsNew = simulation.getTransitionProbabilities();
        labFrameSpinCorrNew = simulation.getLabFrameSpinCorrelationFn();
        holeFrameSpinCorrNew = simulation.getHoleFrameSpinCorrelationFn();
        
        transProbsErr = getTransProbsErr(transProbs, transProbsNew);
        labFrameSpinCorrErr = getSpinCorrErr(labFrameSpinCorr, labFrameSpinCorrNew);
        holeFrameSpinCorrErr = getSpinCorrErr(holeFrameSpinCorr, holeFrameSpinCorrNew);
        
        printErrors(transProbsErr, labFrameSpinCorrErr, holeFrameSpinCorrErr);
        if(transProbsErr<accuracy && labFrameSpinCorrErr<accuracy && holeFrameSpinCorrErr<accuracy) {
            break;
        }
        
        transProbs = transProbsNew;
        labFrameSpinCorr = labFrameSpinCorrNew;
        holeFrameSpinCorr = holeFrameSpinCorrNew;
    }
    cout << endl;
    
    return simulation.getTransitionProbabilities();
}


// Propagation of the hole in an infinite temperature spin environment, with infinitely many spin degrees of freedom
// It only works at times t < 2.0, since it generates the same set of forward and backward paths. At around t=2.0, the
// phase space explodes, and we are undersampling it. Therefore, the correlations in the two sets spoil the results.
//
// Monte Carlo calculation, calling holePropagationInfiniteSpin::run() nThreads times, and then averaging
// the transition probabilities and normalizing them to 1. We also perform a symmetry averaging
// of the transition probabilities, so that they reflect the symmetry properties of the lattice.
vector<double> MC_hole_infiniteT_infiniteSpin(double t, int numberOfRuns, int nFwdPaths, double accuracy) {
    // Define the class that will run the simulation and store our data
    holePropagationInfiniteSpin simulation(t, nFwdPaths);
    
    // Simulation results
    vector<double> transProbs, transProbsNew;
    
    // Errors w.r.t the last run
    double transProbsErr;
    
    // Load data of previous runs, if any
    if(!simulation.load()) {
        simulation.saveAllFilenames();
        simulation.save_returnProb_and_RMSdistance();
    }
    
    // Run simulation 'numberOfRuns' times, and save partial results after each run
    while(simulation.getNumberOfRuns() < numberOfRuns) {
        simulation.run();
        simulation.backup();
        simulation.save();
        
        // Check errors and stop if they are smaller than 'accuracy'
        transProbsNew = simulation.getTransitionProbabilities();
        
        transProbsErr = getTransProbsErr(transProbs, transProbsNew);
        
        printErrors(transProbsErr, 0.0, 0.0);
        if(transProbsErr<accuracy) {
            break;
        }
        
        transProbs = transProbsNew;
    }
    cout << endl;
    
    return simulation.getTransitionProbabilities();
}


// Propagation of the hole in an infinite temperature spin environment, with infinitely many spin degrees of freedom
// As compared to MC_hole_infiniteT_infiniteSpin(), it allows one to perform simulations at longer times as well (t > 2.0).
// It achieves this by generating a different set of forward and backward paths.
//
// Monte Carlo calculation, calling holePropagationInfiniteSpin::run() nThreads times, and then averaging
// the transition probabilities and normalizing them to 1. We also perform a symmetry averaging
// of the transition probabilities, so that they reflect the symmetry properties of the lattice.
vector<double> MC_hole_infiniteT_infiniteSpinFwdBwd(double t, int numberOfRuns, int nFwdPaths, double accuracy) {
    // Define the class that will run the simulation and store our data
    holePropagationInfiniteSpinFwdBwd simulation(t, nFwdPaths);
    
    // Simulation results
    vector<double> transProbs, transProbsNew;
    
    // Errors w.r.t the last run
    double transProbsErr;
    
    // Load data of previous runs, if any
    if(!simulation.load()) {
        simulation.saveAllFilenames();
        simulation.save_returnProb_and_RMSdistance();
    }
    
    // Run simulation 'numberOfRuns' times, and save partial results after each run
    while(simulation.getNumberOfRuns() < numberOfRuns) {
        simulation.run();
        simulation.backup();
        simulation.save();
        
        // Check errors and stop if they are smaller than 'accuracy'
        transProbsNew = simulation.getTransitionProbabilities();
        
        transProbsErr = getTransProbsErr(transProbs, transProbsNew);
        
        printErrors(transProbsErr, 0.0, 0.0);
        if(transProbsErr<accuracy) {
            break;
        }
        
        transProbs = transProbsNew;
    }
    cout << endl;
    
    return simulation.getTransitionProbabilities();
}



// Propagation of the hole in the Neel state
vector<double> MC_hole_Neel(double t, int numberOfRuns, int nFwdPaths, bool calcSpinQ, double accuracy) {
    // Set the initial Neel state, with a hole at the origin
    singleHoleQState neelState;
    neelState.Neel();
    
    // Define the class that will run the simulation and store our data
    holePropagationQState simulation(t, nFwdPaths, neelState, "Neel", calcSpinQ);
    
    // Simulation results
    vector<double> transProbs, transProbsNew;
    vector< vector<double> > labFrameSpinCorr, labFrameSpinCorrNew;
    vector< vector<double> > holeFrameSpinCorr, holeFrameSpinCorrNew;
    
    // Errors w.r.t the last run
    double transProbsErr;
    double labFrameSpinCorrErr;
    double holeFrameSpinCorrErr;
    
    // Load data of previous runs, if any
    if(!simulation.load()) {
        simulation.saveAllFilenames();
        simulation.save_returnProb_and_RMSdistance();
    }
    
    // Run simulation 'numberOfRuns' times, and save partial results after each run
    while(simulation.getNumberOfRuns() < numberOfRuns) {
        simulation.run();
        simulation.backup();
        simulation.save();
        
        // Check errors and stop if they are smaller than 'accuracy'
        transProbsNew = simulation.getTransitionProbabilities();
        labFrameSpinCorrNew = simulation.getLabFrameSpinCorrelationFn();
        holeFrameSpinCorrNew = simulation.getHoleFrameSpinCorrelationFn();
        
        transProbsErr = getTransProbsErr(transProbs, transProbsNew);
        labFrameSpinCorrErr = getSpinCorrErr(labFrameSpinCorr, labFrameSpinCorrNew);
        holeFrameSpinCorrErr = getSpinCorrErr(holeFrameSpinCorr, holeFrameSpinCorrNew);
        
        printErrors(transProbsErr, labFrameSpinCorrErr, holeFrameSpinCorrErr);
        if(transProbsErr<accuracy && labFrameSpinCorrErr<accuracy && holeFrameSpinCorrErr<accuracy) {
            break;
        }
        
        transProbs = transProbsNew;
        labFrameSpinCorr = labFrameSpinCorrNew;
        holeFrameSpinCorr = holeFrameSpinCorrNew;
    }
    cout << endl;
    
    return simulation.getTransitionProbabilities();
}


// Propagation of the hole in the ferromagnetic state
vector<double> MC_hole_Ferromagnet(double t, int numberOfRuns, int nFwdPaths, bool calcSpinQ, double accuracy) {
    // Set the initial Neel state, with a hole at the origin
    singleHoleQState fmState;
    fmState.Ferromagnet();
    
    // Define the class that will run the simulation and store our data
    holePropagationQState simulation(t, nFwdPaths, fmState, "Ferromagnet", calcSpinQ);
    
    // Simulation results
    vector<double> transProbs, transProbsNew;
    vector< vector<double> > labFrameSpinCorr, labFrameSpinCorrNew;
    vector< vector<double> > holeFrameSpinCorr, holeFrameSpinCorrNew;
    
    // Errors w.r.t the last run
    double transProbsErr;
    double labFrameSpinCorrErr;
    double holeFrameSpinCorrErr;
    
    // Load data of previous runs, if any
    if(!simulation.load()) {
        simulation.saveAllFilenames();
        simulation.save_returnProb_and_RMSdistance();
    }
    
    // Run simulation 'numberOfRuns' times, and save partial results after each run
    while(simulation.getNumberOfRuns() < numberOfRuns) {
        simulation.run();
        simulation.backup();
        simulation.save();
        
        // Check errors and stop if they are smaller than 'accuracy'
        transProbsNew = simulation.getTransitionProbabilities();
        labFrameSpinCorrNew = simulation.getLabFrameSpinCorrelationFn();
        holeFrameSpinCorrNew = simulation.getHoleFrameSpinCorrelationFn();
        
        transProbsErr = getTransProbsErr(transProbs, transProbsNew);
        labFrameSpinCorrErr = getSpinCorrErr(labFrameSpinCorr, labFrameSpinCorrNew);
        holeFrameSpinCorrErr = getSpinCorrErr(holeFrameSpinCorr, holeFrameSpinCorrNew);
        
        printErrors(transProbsErr, labFrameSpinCorrErr, holeFrameSpinCorrErr);
        if(transProbsErr<accuracy && labFrameSpinCorrErr<accuracy && holeFrameSpinCorrErr<accuracy) {
            break;
        }
        
        transProbs = transProbsNew;
        labFrameSpinCorr = labFrameSpinCorrNew;
        holeFrameSpinCorr = holeFrameSpinCorrNew;
    }
    cout << endl;
    
    return simulation.getTransitionProbabilities();
}




// ----------------------------
// MONITORING AND SAVING ERRORS
// ----------------------------
// Returns the RMS error of transition probabilities
double getTransProbsErr(vector<double> transProbs1, vector<double> transProbs2) {
    double err, diff;
    int i;
    
    // if the arrays are not of the same size, then we just return infinity
    if(transProbs1.size() != transProbs2.size()) {
        return numeric_limits<double>::infinity();
    }
    
    // determine RMS error
    for(err=0., i=0; i<transProbs1.size(); i++) {
        diff = transProbs1[i] - transProbs2[i];
        err += diff * diff;
    }
    err /= (double) transProbs1.size();
    err = sqrt(err);
    
    return err;
}

// Returns the RMS error of spin correlations
double getSpinCorrErr(vector< vector<double> > spinCorr1, vector< vector<double> > spinCorr2) {
    double err, diff;
    int i1, i2;
    int n;
    
    // if the arrays are not of the same size, then we just return infinity
    if(spinCorr1.size() != spinCorr2.size()) {
        return numeric_limits<double>::infinity();
    }

    // determine RMS error
    for(err=0., n=0, i1=0; i1<spinCorr1.size(); i1++) {
        for(i2=0; i2<spinCorr1[i1].size(); i2++) {
            diff = spinCorr1[i1][i2] - spinCorr2[i1][i2];
            err += diff * diff;
        }
        n += spinCorr1[i1].size();
    }
    err /= (double) n;
    err = sqrt(err);
    
    return err;
}

// Print errors
void printErrors(double transProbsErr, double labFrameSpinCorrErr, double holeFrameSpinCorrErr) {
    cout << endl << endl << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "RMS errors with respect to the last run:" << endl;
    cout << "Err(transition probabilities):     \t" << setprecision(5) << fixed << transProbsErr << endl;
    cout << "Err(lab frame spin correlations):  \t" << setprecision(5) << fixed << labFrameSpinCorrErr << endl;
    cout << "Err(hole frame spin correlations): \t" << setprecision(5) << fixed << holeFrameSpinCorrErr << endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << endl << endl;
}

