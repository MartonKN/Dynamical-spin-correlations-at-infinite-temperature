//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_singleHoleQStateBin.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/19/16.
//
//

#ifndef MC_singleHoleQStateBin_hpp
#define MC_singleHoleQStateBin_hpp

#include <stdio.h>


#include <stdio.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <set>
#include "MC_clusterRunQ.hpp"
#include "MC_singleHoleQState.hpp"
#include "MC_pathData.hpp"

using namespace std;




// --------------------------------------------------------------------------
// CONTAINER OF SPIN STATES (PREFERABLY, BUT NOT NECESSARILY AT A GIVEN SITE)
// --------------------------------------------------------------------------
// The class 'singleHoleQStateContainer' has a similar functionality, but it turned out to consume too much memory.
// Since in many cases we just need to convert 'pathData' elements ending at a given site into 'singleHoleQState'-s,
// we need a container to temporarily store all those states.
//
// Stores quantum states in a red-black self balanced binary search tree provided by 'std::set', to ensure easy and fast access
class singleHoleQStateBin {
private:
    // Stores the sets of 'singleHoleQState'-s
    set<singleHoleQState> bin;
public:
    // Constructor and destructor
    singleHoleQStateBin(void) {InitializeRandomNumberGenerator(); this->clear();};
    ~singleHoleQStateBin(void) {this->clear();};
    
    // Get basic data
    set<singleHoleQState>::iterator begin(void) const {return this->bin.begin();};
    set<singleHoleQState>::iterator end(void) const {return this->bin.end();};
    
    // clear all states from the container
    inline void clear(void) {this->bin.clear(); };
    
    // get memory allocation in bytes
    long int getMemoryBytes(void) const;
    
    
    // Fill container with quantum states
    void fill(const set<pathData> &sPD, const singleHoleQState & qSInitial);     // fills '*this' with quantum states that were
                                                                                 // permuted according to the action of the hole in
                                                                                 // pathContainer pC.
                                                                                 // It deletes all previous entries from '*this'
    void add (const set<pathData> &sPD, const singleHoleQState & qSInitial);     // Does the same thing as fill, but without deleting
                                                                                 // all previous entries first
    
    // Searching in the bin
    inline set<singleHoleQState>::iterator begin(const singleHoleQState & qS) {return this->bin.begin();};
    inline set<singleHoleQState>::iterator end(const singleHoleQState & qS) {return this->bin.end();};
    inline set<singleHoleQState>::iterator find(const singleHoleQState & qS) {return this->bin.find(qS);};
    
    
    friend ostream & operator<<(ostream &out, const singleHoleQStateBin &qSBin);// print the contents of the bin
};

#endif /* MC_singleHoleQStateBin_hpp */
