//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_holePropagationSpinStateAvgFwdBwd.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 11/19/16.
//
//
//  The class in this file simulates the propagation of the hole in an infinite temperature spin
//  environment. The only difference between this and the class 'holePropagationSpinStateAvg' defined
//  in 'MC_holePropagationSpinStateAvg.hpp' is that in this case, we generate a different set of forward
//  and backward paths.
//
//  IMPORTANT: we took something for granted during this implementation, to avoid having to rewrite most of the code.
//  When we take the overlap of quantum states fwdState and bwdState (which are obtained from the initial state by
//  transforming it using a forward and a backward path), we should update the transition probabilities and spin
//  correlations by <bwdState|fwdState>. However, we update them instead only by the real part of this expression
//  real(<bwdState|fwdState>). This would be correct in the limit of infinitely many paths, since the imaginary
//  parts should cancel out. However, this is not necessarly true in the case of finitely many paths, and the
//  cancelation of imaginary parts would be a useful sanity check.
//  However, implementing this would be a huge hassle, and would require me to write a completely new base class
//  instead of 'holePropagation', that can handle complex bins for the transition probabilities. Since probably
//  it's not worth the time, I staid with the above approach, that is compatible with 'holePropagation'.

#ifndef MC_holePropagationSpinStateAvgFwdBwd_hpp
#define MC_holePropagationSpinStateAvgFwdBwd_hpp

#include <stdio.h>
#include "MC_clusterRunQ.hpp"
#include "MC_holePropagation.hpp"
#include "MC_singleHoleQState.hpp"
#include "MC_singleHoleQStateBin.hpp"

using namespace std;


// -------------------------------------------------------
// HOLE PROPAGATION IN AN INFINITE TEMPERATURE ENVIRONMENT
// -------------------------------------------------------
class holePropagationSpinStateAvgFwdBwd : public holePropagation {
    int nInitialStates;   // the number of initial quantum states of the spin system we will be averaging over
    
    
    // BASIC CONSTRUCTORS AND DESTRUCTORS (constructors also initialize the random number generator)
public:
    holePropagationSpinStateAvgFwdBwd(void) : holePropagation() {
        this->nInitialStates=0;
        this->generateSaveFilenames();
    };
    holePropagationSpinStateAvgFwdBwd(double t, int nFwdPaths, int numberOfInitialSpinStates, bool spinCorrQ=true)
    : holePropagation(t, nFwdPaths,spinCorrQ) {
        this->nInitialStates=numberOfInitialSpinStates;
        this->generateSaveFilenames();
    };
    ~holePropagationSpinStateAvgFwdBwd(void) {};
    
    // GET BASIC DATA
public:
    inline int getNumberOfInitialStates(void) const {return this->nInitialStates;};
    
    // MONTE CARLO CALCULATION
public:
    // First, generate two new 'pathContainer'-s storing the permutation cycles generated by the forward and backward paths.
    // Then, take 'nInitialStates' different random initial spin configurations, and transform that state using
    // the paths in the containers. The benefit of this approach is that we do not need to generate the paths
    // again for each new initial state, and we spare a lot of time as compared to Carlström's original simulation.
    void run(void);
    
    // same as run(void), just does not automatically generate the pathContainers, but uses 'pCfwd' and 'pCbwd' instead
    void run(const pathContainer &pCfwd, const pathContainer &pCbwd);
    
private:
    // Updates the bins with contributions from the quantum states 'qSFwd' and qSBwd
    // The function using this routine is supposed to make sure that the spin states in
    // 'qSFwd' and 'qSBwd' are the same, and they might only differ in their coefficients.
    // This is not checked in 'update' not to waste resources.
    // In this class, 'update' is currently only used by 'run(const pathContainer &, const pathContainer &)'.
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
    double update(const singleHoleQState &qSfwd, const singleHoleQState &qSbwd);
    
    
    // DETERMINE SPIN CORRELATIONS
public:
    // calculate spin correlations using LabFrameSpinCorrelationBins (Overwritten function!)
    void calculateLabFrameSpinCorrelationFn (vector< vector<double> > &) const;
    
    // calculate spin correlations using HoleFrameSpinCorrelationBins (Overwritten function!)
    void calculateHoleFrameSpinCorrelationFn(vector< vector<double> > &) const;
    
    // OVERLOADED OPERATORS
public:
    // add the values of the bins in hP to *this, and add their nRuns values if the basic parameters of hP agree with *this
    holePropagationSpinStateAvgFwdBwd & operator+= (const holePropagationSpinStateAvgFwdBwd & hP);
    
    // assignment operator
    holePropagationSpinStateAvgFwdBwd & operator=  (const holePropagationSpinStateAvgFwdBwd & hP);
    
    
    
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



#include <stdio.h>

#endif /* MC_holePropagationSpinStateAvgFwdBwd_hpp */
