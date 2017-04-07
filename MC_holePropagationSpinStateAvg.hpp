//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_holePropagationSpinStateAvg.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/19/16.
//
//
//  The class in this file simulates the propagation of the hole in an infinite temperature spin
//  environment. It is very similar to Carlström's simulation in the sense that it performs averaging
//  over initial quantum states, instead of performing the exact spin trace, as Izabella did.

#ifndef MC_holePropagationSpinStateAvg_hpp
#define MC_holePropagationSpinStateAvg_hpp

#include <stdio.h>
#include "MC_clusterRunQ.hpp"
#include "MC_holePropagation.hpp"
#include "MC_singleHoleQState.hpp"
#include "MC_singleHoleQStateContainer.hpp"

using namespace std;


// -------------------------------------------------------
// HOLE PROPAGATION IN AN INFINITE TEMPERATURE ENVIRONMENT
// -------------------------------------------------------
class holePropagationSpinStateAvg : public holePropagation {
    int nInitialStates;   // the number of initial quantum states of the spin system we will be averaging over
    
    
    // BASIC CONSTRUCTORS AND DESTRUCTORS (constructors also initialize the random number generator)
public:
    holePropagationSpinStateAvg(void) : holePropagation() {
        this->nInitialStates=0;
        this->generateSaveFilenames();
    };
    holePropagationSpinStateAvg(double t, int nFwdPaths, int numberOfInitialSpinStates, bool spinCorrQ=true) : holePropagation(t, nFwdPaths,spinCorrQ) {
        this->nInitialStates=numberOfInitialSpinStates;
        this->generateSaveFilenames();
    };
    ~holePropagationSpinStateAvg(void) {};
    
    // GET BASIC DATA
public:
    int getNumberOfInitialStates(void) const {return nInitialStates;};
    
    // MONTE CARLO CALCULATION
public:
    // First, generate new pathContainer storing the permutation cycles generated by the forward paths
    // Then, take 'nInitialStates' different random initial spin configurations, and transform that state using
    // the paths in the container. The benefit of this approach is that we do not need to generate the paths
    // again for each new initial state, and we spare a lot of time as compared to Carlström's original simulation.
    // Furthermore, it is expected to be faster than our previous simulation as well, when the number of paths is
    // large. That simulation suffered from the problem that we had to take pairs of paths to calculate the transition
    // probabilities, which was not feasible once the number of paths became of the order of 100K.
    // In contrast, in the current simulation, the number of quantum states that we generate is only of the
    // order of the paths (~50% of them), and we do not need to take pairs.
    void run(void);
    
    // same as run(void), just does not automatically generate the container, but uses 'pC' instead
    void run(const pathContainer &pC);
    
private:
    // Updates the bins with contributions from the quantumState 'qS'
    double update(const singleHoleQState &qS);
    
    
    // DETERMINE SPIN CORRELATIONS
public:
    // calculate spin correlations using LabFrameSpinCorrelationBins (Overwritten function!)
    void calculateLabFrameSpinCorrelationFn (vector< vector<double> > &) const;
    
    // calculate spin correlations using HoleFrameSpinCorrelationBins (Overwritten function!)
    void calculateHoleFrameSpinCorrelationFn(vector< vector<double> > &) const;
    
    // OVERLOADED OPERATORS
public:
    // add the values of the bins in hP to *this, and add their nRuns values if the basic parameters of hP agree with *this
    holePropagationSpinStateAvg & operator+= (const holePropagationSpinStateAvg & hP);
    
    // assignment operator
    holePropagationSpinStateAvg & operator=  (const holePropagationSpinStateAvg & hP);

    
    
    // LOAD EARLIER BACKUP DATA
public:
    // loads parameters and returns true only if they match those of '*this'
    // (Overloaded function!)
    bool loadParameters(void);

    // SAVE DATA
public:
    // We will need to indicate in the save filenames that we do simulation in the specified quantum state,
    // therefore generateSaveFilenames() is unique to the child class.
    // (Overwritten function!)
    void generateSaveFilenames(void);
    
    // save tHole, nPaths, LATTICE_SIZE, NUMBER_OF_SPINS, returnProbability, RMS distance, nInitialStates
    // (Overloaded function!)
    bool saveParameters(void) const;
};



#endif /* MC_holePropagationSpinStateAvg_hpp */