//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_Monte_Carlo_simulations.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/18/16.
//
//

#ifndef MC_Monte_Carlo_simulations_hpp
#define MC_Monte_Carlo_simulations_hpp

#include <stdio.h>
#include "MC_clusterRunQ.hpp"
#include "MC_holePropagation.hpp"
#include "MC_holePropagationTraceFwdBwd.hpp"
#include "MC_holePropagationQState.hpp"
#include "MC_holePropagationSpinStateAvg.hpp"
#include "MC_holePropagationSpinStateAvgFwdBwd.hpp"
#include "MC_holePropagationInfiniteSpin.hpp"
#include "MC_holePropagationInfiniteSpinFwdBwd.hpp"

using namespace std;




// -----------------------
// MONTE CARLO SIMULATIONS
// -----------------------
// Input parameters:
// t = hole propagation time
// numberOfRuns = number of independent runs, that we will average over
// nFwdPaths = number of randomly distributed forward paths
// accuracy = stop the simulation if the RMS error of the transition probs and spin correlations w.r.t. the last run become smaller than this


// Propagation of the hole in an infinite temperature spin environment
//
// Monte Carlo calculation, calling holePropagation::run() nThreads times, and then averaging
// the transition probabilities and normalizing them to 1. We also perform a symmetry averaging
// of the transition probabilities, so that they reflect the symmetry properties of the lattice.
vector<double> MC_hole_infiniteT(double t, int numberOfRuns, int nFwdPaths, bool calcSpinQ=true, double accuracy = 0.000001);


// Propagation of the hole in an infinite temperature spin environment
//
// Same as MC_hole_infiniteT(), except that it takes a different set for forward and backward propagating paths.
// This doubles the runtime and the memory requirement, but it also means that the code can run for times t>2.
// At these times, MC_hole_infiniteT gave unreliable results, even for 200 million paths.
vector<double> MC_hole_infiniteT_TraceFwdBwd(double t, int numberOfRuns, int nFwdPaths, bool calcSpinQ=true, double accuracy = 0.000001);



// Propagation of the hole in an infinite temperature spin environment
//
// Monte Carlo calculation, calling holePropagationSpinStateAvg::run() nThreads times.
// Instead of performing a trace over spin configurations exactly, it samples random initial
// spin configurations 'numberOfInitialSpinStates' times, and averages the quantities we are
// interested in over this spin configuration space
vector<double> MC_hole_infiniteTSpinStateAvg(double t, int numberOfRuns, int nFwdPaths, int numberOfInitialSpinStates, bool calcSpinQ=true, double accuracy = 0.000001);



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
vector<double> MC_hole_infiniteTSpinStateAvgFwdBwd(double t, int numberOfRuns, int nFwdPaths, int numberOfInitialSpinStates, bool calcSpinQ=true, double accuracy = 0.00001);



// Propagation of the hole in an infinite temperature spin environment, with infinitely many spin degrees of freedom
// It only works at times t < 2.0, since it generates the same set of forward and backward paths. At around t=2.0, the
// phase space explodes, and we are undersampling it. Therefore, the correlations in the two sets spoil the results.
//
// Monte Carlo calculation, calling holePropagationInfiniteSpin::run() nThreads times, and then averaging
// the transition probabilities and normalizing them to 1. We also perform a symmetry averaging
// of the transition probabilities, so that they reflect the symmetry properties of the lattice.
vector<double> MC_hole_infiniteT_infiniteSpin(double t, int numberOfRuns, int nFwdPaths, double accuracy = 0.00001);


// Propagation of the hole in an infinite temperature spin environment, with infinitely many spin degrees of freedom
// As compared to MC_hole_infiniteT_infiniteSpin(), it allows one to perform simulations at longer times as well (t > 2.0).
// It achieves this by generating a different set of forward and backward paths.
//
// Monte Carlo calculation, calling holePropagationInfiniteSpin::run() nThreads times, and then averaging
// the transition probabilities and normalizing them to 1. We also perform a symmetry averaging
// of the transition probabilities, so that they reflect the symmetry properties of the lattice.
vector<double> MC_hole_infiniteT_infiniteSpinFwdBwd(double t, int numberOfRuns, int nFwdPaths, double accuracy = 0.00001);



// Propagation of the hole in the Neel state
vector<double> MC_hole_Neel(double t, int numberOfRuns, int nFwdPaths, bool calcSpinQ=true, double accuracy = 0.00001);



// Propagation of the hole in the ferromagnetic state
vector<double> MC_hole_Ferromagnet(double t, int numberOfRuns, int nFwdPaths, bool calcSpinQ=true, double accuracy = 0.00001);


// ----------------------------
// MONITORING AND SAVING ERRORS
// ----------------------------
// Returns the RMS error of transition probabilities
double getTransProbsErr(vector<double> transProbs1, vector<double> transProbs2);

// Returns the RMS error of spin correlations
double getSpinCorrErr(vector< vector<double> > spinCorr1, vector< vector<double> > spinCorr2);

// Print errors
void printErrors(double transProbsErr, double labFrameSpinCorrErr, double holeFrameSpinCorrErr);


#endif /* MC_Monte_Carlo_simulations_hpp */
