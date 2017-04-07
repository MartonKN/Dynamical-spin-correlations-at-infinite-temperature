//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_holePropagationTraceFwdBwd.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 1/10/17.
//
//

#ifndef MC_holePropagationTraceFwdBwd_hpp
#define MC_holePropagationTraceFwdBwd_hpp

#include <stdio.h>
#include <string>
#include "MC_clusterRunQ.hpp"
#include "MC_holePropagation.hpp"

using namespace std;


// -----------------------------------------------------------------------------------------------------
// HOLE PROPAGATION IN AN INFINITE TEMPERATURE ENVIRONMENT, WITH SPIN 1/2 DEGREES OF FREEDOM
// We calculate the spin correlations and the transition probabilities by evaluating the many-body
// trace exactly. Furthermore, we use a different set for forward and backward propagating paths.
// This is in contrast to the parent class holePropagation, which uses the same set for forward and
// backward propagation. [Although that is faster, it suffers from the correlations between the paths
// at long times (t>2.0), when the phase space becomes very large.]
// This version of the code allows us to go to longer times, (t>2.0) as well.
// -----------------------------------------------------------------------------------------------------
class holePropagationTraceFwdBwd : public holePropagation {
    // BASIC CONSTRUCTORS AND DESTRUCTORS (constructors also initialize the random number generator)
public:
    holePropagationTraceFwdBwd(void) : holePropagation() {
        this->generateSaveFilenames();
    };
    holePropagationTraceFwdBwd(double t, int nFwdPaths, bool spinCorrQ=true) : holePropagation(t, nFwdPaths, spinCorrQ) {
        this->generateSaveFilenames();
    };
    ~holePropagationTraceFwdBwd(void) {};
    
    // MONTE CARLO CALCULATION
public:
    virtual void run(void);                   // generates new pathContainer and updates the bins using holePropagation::update()
private:
    void   update(const pathContainer &contBwd, const pathContainer &contFwd);
                                              // updates the bins with contributions from all pairs of paths in 'contFwd' and 'contBwd'
    double update(const pathData &pDbwd, const pathData &pDfwd, const vector<double> &invNSpinsPowers);
                                              // Updates the bins with contributions from the pair of paths 'pDfwd' and 'pDbwd',and returns the contribution of this pair to the transition probabilities.
                                              // invNSpinsPowers is an auxiliary array: invNSpinsPowers[i] = 1/NUMBER_OF_SPINS^i, that needs to be determined in update(pathContainer &)
                                              // Note that we should add the transition probability contribution of the pair of paths to all same site spin correlation bins, except for the hole since this would lead to unnecessary computational overhead, we actually remove this contribution from the hole site. We will compensate for this at the end, in getLabFrameSpinCorrelationFn() and getHoleFrameSpinCorrelationFn().
    
    // SAVE DATA
    // We will need to indicate in the save filenames that we do simulation in the specified quantum state,
    // therefore generateSaveFilenames() is unique to the child class. (Overwritten function!)
public:
    void generateSaveFilenames(void);
};

#endif /* MC_holePropagationTraceFwdBwd_hpp */
