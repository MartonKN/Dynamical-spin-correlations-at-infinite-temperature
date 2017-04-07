//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_holePropagation.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/6/16.
//
//

#ifndef MC_holePropagation_hpp
#define MC_holePropagation_hpp

#include <stdio.h>
#include <cmath>
#include <complex>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "MC_clusterRunQ.hpp"
#include "MC_random.hpp"
#include "MC_coord.hpp"
#include "MC_Permutation.hpp"
#include "MC_PermutationCycles.hpp"
#include "MC_pathData.hpp"
#include "MC_pathContainer.hpp"


using namespace std;


// ------------------------------------------------------------------
// CLASS RUNNING THE INFINITE TEMPERATURE SIMULATION AND STORING DATA
// ------------------------------------------------------------------
// Stores the transition probabilities as well as the spin correlation functions, both in the laboratory frame,
// as well as in the hole's co-moving frame
class holePropagation {
protected: 
    static const int NUMBER_OF_SPINS = 2;                        // how many spin components there are
    static const int nSites = 4 * LATTICE_SIZE * LATTICE_SIZE;   // number of lattice sites
    double tHole;                                                // propagation time of the hole
    int nPaths;                                                  // number of paths generated at each run of the simulation
    int nRuns;                                                   // number of times the simulation has been run already
    bool calcSpinCorrelationsQ;                                  // whether it should calculate the spin correlations
    vector<long double> transProbsBins;                          // bins storing the transition probability contributions of the paths, size nSites
    vector< vector<long double> > LabFrameSpinCorrelationBins;   // bins for spin correlations in the lab frame, size nSitesxnSites
    vector< vector<long double> > HoleFrameSpinCorrelationBins;  // bins for spin correlations in hole's co-moving frame, size nSitesxnSites
    vector<string> saveFilenames;                            // filenames are dictated by the system parameters, and they are set at construction
                                                             // saveFilenames[0] = file containing all the filenames
                                                             // saveFilenames[1] = system parameters
                                                             // saveFilenames[2] = return probability and rms distance data of each run
                                                             // saveFilenames[3] = transition probabilities
                                                             // saveFilenames[4] = lab frame spin correlation function
                                                             // saveFilenames[5] = hole frame spin correlation function
                                                             // saveFilenames[6] = transProbsBins
                                                             // saveFilenames[7] = LabFrameSpinCorrelationBins
                                                             // saveFilenames[8] = HoleFrameSpinCorrelationBins
    // CONSTRUCTORS, DESTRUCTOR AND BASIC FUNCTIONS
public:
    holePropagation(void);                                          // constructor: fills all arrays with zeros, sets time to zero
    holePropagation(double t, int nFwdPaths, bool spinCorrQ=true);  // constructor: fills all arrays with zeros, sets time to t
    ~holePropagation(void);                                         // destructor
    void reset(double t, int nFwdPaths);                            // resets time, #runs, #paths, and clears all values from the bins
    void clearBins(void);                                           // clears all values from the bins and resets nRuns to zero
    virtual void generateSaveFilenames(void);                       // creates the filenames based on system parameters
    
    // MONTE CARLO CALCULATION
public:
    virtual void run(void);                                  // generates new pathContainer and updates the bins using holePropagation::update()
private:
    void   update(const pathContainer &cont);                // updates the bins with contributions from all pairs of paths in 'cont'
    double update(const pathData &pDbwd, const pathData &pDfwd, const vector<double> &invNSpinsPowers);
                                                             // Updates the bins with contributions from the pair of paths 'pDfwd' and 'pDbwd',and returns the contribution of this pair to the transition probabilities.
                                                             // invNSpinsPowers is an auxiliary array: invNSpinsPowers[i] = 1/NUMBER_OF_SPINS^i, that needs to be determined in update(pathContainer &)
                                                             // Note that we should add the transition probability contributionof the pair of paths to all same site spin correlation bins, except for the hole since this would lead to unnecessary computational overhead, we actually remove this contribution from the hole site. We will compensate for this at the end, in getLabFrameSpinCorrelationFn() and getHoleFrameSpinCorrelationFn().
                                                             // Note furthermore, that on the square lattice, the contribution of (pDfwd, pDbwd) shall be the same as (pDbwd, pDfwd). As holePropagation::update() is only used when we loop over all pathData elements of a pathContainer, it is a factor of 2 waste of time to take all pairs. Instead, we set holePropagation::update() such that it only updates '*this', when pDfwd<=pDbwd, but in this case it takes into account the contributions of both (pDfwd, pDbwd) and (pDbwd, pDfwd), which shall be equal.

    // DETERMINE SPIN CORRELATIONS
public:
    virtual void calculateLabFrameSpinCorrelationFn (vector< vector<double> > &) const; // get spin correlations from LabFrameSpinCorrelationBins
    virtual void calculateHoleFrameSpinCorrelationFn(vector< vector<double> > &) const; // get spin correlations from HoleFrameSpinCorrelationBins
    
    // GET BASIC DATA
public:
    int getNumberOfRuns(void) const;                         // get the number of runs performed so far
    vector<string> getFilenames(void) const;                 // return the file names we will be saving our data in
    vector<double> getTransitionProbabilities(void) const;   // returns the transition probabilities (symmetry averaged)
    vector< vector<double> > getLabFrameSpinCorrelationFn(void) const;  // returns the lab frame spin correlations
    vector< vector<double> > getHoleFrameSpinCorrelationFn(void) const; // returns the hole frame spin correlations
    double getReturnProbability(void) const;                 // probability of the hole to return to the origin
    double getRMSdistance(void) const;                       // root mean squared distance of the hole from the origin
    
    // BACK UP THE BINS
public:
    bool backup(void) const;                            // saves all three bins, using the three saving functions below
    bool saveTransitionProbabilityBins(void) const;     // save transition probability bins
    bool saveLabFrameSpinCorrelationBins(void) const;   // save lab frame spin correlation bins
    bool saveHoleFrameSpinCorrelationBins(void) const;  // save hole frame spin correlation bins
    bool save_returnProb_and_RMSdistance(void) const;   // updates the file storing the return probability and the RMS distance of the hole
                                                        // after each run
    
    // LOAD EARLIER BACKUP DATA
public:
    virtual bool load(void);                            // loads all three bins, using the three saving functions below
    bool loadTransitionProbabilityBins(void);           // load transition probability bins
    bool loadLabFrameSpinCorrelationBins(void);         // load lab frame spin correlation bins
    bool loadHoleFrameSpinCorrelationBins(void);        // load hole frame spin correlation bins
    virtual bool loadParameters(void);                  // loads parameters and returns true only if they match those of '*this'
    
    // SAVE RESULTS
public:
    virtual bool save(void) const;                      // save everything
    virtual bool saveParameters(void) const;            // save tHole, nPaths, LATTICE_SIZE, NUMBER_OF_SPINS, returnProbability, RMS distance
    bool saveAllFilenames(void) const;                  // save the contents of saveFilenames into a file
    bool saveTransitionProbabilities(void) const;       // save transition probabilities of the hole from the origin to each site
    bool saveLabFrameSpinCorrelationFn(void) const;     // save spin correlation function in the lab frame
    bool saveHoleFrameSpinCorrelationFn(void) const;    // save spin correlation function in the hole's comoving frame
    
    // OVERLOADED OPERATORS
public:
    virtual holePropagation & operator+= (const holePropagation & hP);  // add the values of the bins in hP to *this, and add their nRuns values
    virtual holePropagation & operator=  (const holePropagation & hP);  // assignment operator
};

// SYMMETRY AVERAGING OF TRANSITION PROBABILITIES
vector<double> symmetryAverage(vector<double> transProb); // Symmetry averaging of the transition probabilities
int rotate90 (int siteIndex); // 90 degree rotation of lattice indices
int reflect10(int siteIndex); // reflection of site indices to the (1,0) axis
int reflect01(int siteIndex); // reflection of site indices to the (0,1) axis
int reflect11(int siteIndex); // reflection of site indices to the (1,1) axis
vector<double> transformProbs(vector<double> probs, int (*transformation) (int)); // transforms the indices of transProb

// AUXILIARY FUNCTIONS
bool checkIfFileExists(string fname);                   // check if a file already exists





#endif /* MC_holePropagation_hpp */
