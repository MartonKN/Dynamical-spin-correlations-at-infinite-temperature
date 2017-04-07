//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_singleHoleQStateContainer.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/19/16.
//
//
//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.


#ifndef MC_singleHoleQStateContainer_hpp
#define MC_singleHoleQStateContainer_hpp

#include <stdio.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <set>
#include "MC_clusterRunQ.hpp"
#include "MC_singleHoleQState.hpp"
#include "MC_pathContainer.hpp"

using namespace std;




// -------------------------------------
// CONTAINER OF SPIN STATES AT ALL SITES
// -------------------------------------
// Stores quantum states in a red-black self balanced binary search tree provided by 'std::set', to ensure easy and fast access
class singleHoleQStateContainer {
protected:
    // Stores the sets of 'singleHoleQState'-s with the same hole site
    vector< set<singleHoleQState> > container;
public:
    // Constructor and destructor
    singleHoleQStateContainer(void);  // constructor creating empty sets for each site
    ~singleHoleQStateContainer(void); // destructor that clears all sets
    void  clear(void);                // clear all states from the container
    
    // Get basic data
    int size(void) const {return container.size();};  // get the number of sites
    int getNumberOfBins(void) const;                  // get the number of bins (i.e. the number of non-identical paths)
    long int getMemoryBytes(void) const;              // get memory allocation in bytes
    const set<singleHoleQState>& operator[] (const int & siteIndex) const {return container[siteIndex];};  // element (read only)
    
    
    // Fill container with quantum states
    void fill(const pathContainer &pC, const singleHoleQState &qSInitial);     // fills '*this' with quantum states that were
    // permuted according to the action of the hole in
    // pathContainer pC.
    // It deletes all previous entries from '*this'
    void add (const pathContainer &pC, const singleHoleQState &qSInitial);     // does the same thing as fill, but without deleting
    // all previous entries first
    
    // Printing
    void printBin(const int siteIndex) const;                                        // prints the elements of the bin at site 'siteIndex'
    friend ostream & operator<<(ostream &out, const singleHoleQStateContainer &qSC); // print the entire container
};


#endif /* MC_singleHoleQStateContainer_hpp */
