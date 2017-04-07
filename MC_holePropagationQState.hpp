//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com 
//
//  MC_holePropagationQState.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/17/16.
//
//
// The class defined here is a child of holePropagation, and simulates the propagation of the hole
// in an initial quantum state

#ifndef MC_holePropagationQState_hpp
#define MC_holePropagationQState_hpp

#include <stdio.h>
#include "MC_clusterRunQ.hpp"
#include "MC_holePropagation.hpp"
#include "MC_pathContainer.hpp"
#include "MC_singleHoleQState.hpp"
#include "MC_singleHoleQStateContainer.hpp"

using namespace std;



// --------------------------------------------
// HOLE PROPAGATION IN AN INITIAL QUANTUM STATE
// --------------------------------------------
class holePropagationQState : public holePropagation {
private:
    singleHoleQState QState; // initial state
public:
    // BASIC CONSTRUCTORS AND DESTRUCTORS
    // All of these set the coefficient of QState to 1.
    holePropagationQState(void);
    holePropagationQState(double t, int nFwdPaths, const singleHoleQState &inputQState, string nameOfState, bool spinCorrQ=true);
    ~holePropagationQState(void);
    void resetQState(double t, int nFwdPaths, const singleHoleQState &newQState, string nameOfState); // reset time, #runs, #paths, QState,
                                                                                                      // and clears all values from the bins

    
    // MONTE CARLO CALCULATION
public:
    // generate new singleHoleQStateContainer and updates the bins using holePropagationQState::update() (Overwritten function!)
    void run(void);
    
    // same as run(void), just does not automatically generate the container, but uses 'qSC' instead
    void run(const pathContainer & pathCont);
    
private:
    // Calculates the contribution of each path in pathContainer to the transition probabilities, and it updates transProbsBins[]
    void update(const pathContainer & pathCont);
    
    // Updates the bins with contributions from the quantumState 'qS'
    double update(const singleHoleQState & qState);
    
    
    // DETERMINE SPIN CORRELATIONS
public:
    // calculate spin correlations using LabFrameSpinCorrelationBins (Overwritten function!)
    void calculateLabFrameSpinCorrelationFn (vector< vector<double> > &) const;
    
    // calculate spin correlations using HoleFrameSpinCorrelationBins (Overwritten function!)
    void calculateHoleFrameSpinCorrelationFn(vector< vector<double> > &) const;


    // LOAD EARLIER BACKUP DATA
public:
    // load all three bins, using the saving functions of the parent class holePropagation (Overwritten function!)
    bool load(void);
    
    bool loadParameters(void) { return this->holePropagation::loadParameters(); };
    
private:
    // load the initial quantum state 'QState' we use during the simulation
    bool loadQState(void);
   
    
    // SAVE DATA
public:
    // We will need to indicate in the save filenames that we do simulation in the specified quantum state,
    // therefore generateSaveFilenames() is unique to the child class.
    void generateSaveFilenames(string fnameBase = "");
    //void generateSaveFilenames() { generateSaveFilenames(""); }; // (Overwritten function!)
    
    // save everything, including the initial state 'QState' of the system (Overwritten function!)
    bool save(void) const;
    
private:    
    // save the initial quantum state 'QState' we use during the simulation
    bool saveQState(void) const;
    
    
    // OVERLOADED OPERATORS
public:
    holePropagationQState & operator+= (const holePropagationQState & hP); // add the values of the bins in hP to *this, and add their nRuns values
    holePropagationQState & operator=  (const holePropagationQState & hP); // assignment operator
};


#endif /* MC_holePropagationQState_hpp */
