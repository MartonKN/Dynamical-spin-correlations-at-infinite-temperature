//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_pathContainer.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/19/16.
//
//

#ifndef MC_pathContainer_hpp
#define MC_pathContainer_hpp

#include <stdio.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <set>
#include "MC_clusterRunQ.hpp"
#include "MC_coord.hpp"
#include "MC_pathData.hpp"

using namespace std;


// ----------------------
// CONTAINER OF PATH DATA
// ----------------------
// Stores path data in a red-black self balanced binary search tree provided by 'std::set', to ensure easy and fast access
class pathContainer {
protected:
    // Stores the sets of paths ending at the same site.
    // The vector container is indexed by the site indices of the endpoint of the paths stored in it
    vector< set<pathData> > container;
    double tHole;                       // hole propagation time, set when we fill the container with paths
    long int nPaths;                    // number of paths stored in the container
public:
    pathContainer(void);                // constructor (also initializes the random number generator)
    ~pathContainer(void);               // destructor that clears all sets
    const set<pathData>& operator[] (const int siteIndex) const; // provides a handle to the set of paths ending at siteIndex (read only)
    void fill(double t, int nFwdPaths); // empty, then fill the container with nFwdPaths elements, randomly generated using pathData::random
    void clear();                       // clear all paths from the container
    int getNumberOfBins() const;        // get the number of bins (i.e. the number of non-identical paths)
    long int getMemoryBytes() const;    // get memory allocation in bytes
    int getNumberOfPaths(void) const {return nPaths;};
    int getTime(void) const {return tHole;};
    void printBin(int siteIndex) const; // prints bin corresponding to site 'siteIndex' to 'cout'
    friend ostream & operator<< (ostream & out, const pathContainer & cont); // print all bins that are not empty
};



#endif /* MC_pathContainer_hpp */
