//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_pathContainer.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/19/16.
//
//

#include "MC_pathContainer.hpp"



// ----------------------
// CONTAINER OF PATH DATA
// ----------------------
// constructor
pathContainer::pathContainer(void) {
    InitializeRandomNumberGenerator();
    this->tHole = 0.;
    this->nPaths = 0;
    this->container.reserve(4 * LATTICE_SIZE * LATTICE_SIZE);
    this->container.resize(4 * LATTICE_SIZE * LATTICE_SIZE);
    this->container.shrink_to_fit();
}


// destructor that clears all sets
pathContainer::~pathContainer(void) {
    this->clear();
}


// provides the handle to the set containing paths ending at siteIndex (read only)
const set<pathData>& pathContainer::operator[] (const int siteIndex) const {
    return container[siteIndex];
}


// fills the container with nFwdPaths elements, randomly generated using pathData::random.
// previous paths are deleted
// the lengths of the paths are chosen according to a Poisson distribution of mean z*t, where z=4 is the coordination number
// and t is the propgation time of the hole
void pathContainer::fill(double t, int nFwdPaths) {
    this->clear();
    this->tHole = t;             // set the time of the container to 't', so that we remember later
    this->nPaths = nFwdPaths;    // set the number of paths in the container to 'nFwdPaths'
    
    int z=4;                // coordination number of the lattice
    int len;                // length of the newly generated path
    int i;                  // iterator
    pathData pD;            // pathData element to be added to the appropriate bin
    int endPointLatestPath; // endpoint of the latest path that we generated
    complex<double> coefficientLatestPath;  // coefficient of the latest path that we generated -- should be (-ii)^len
    int onePercent;         // one percent of the total number of paths (we display a '.' when a percent is ready)
    
    onePercent = nFwdPaths/100;
    if(onePercent==0) {
        onePercent = 1;
    }
    
    // We return a pair of values to be returned by 'std::set' when we insert a new element of pathData
    // pair::first is set to an iterator pointing to either the newly inserted element or to the equivalent element already in the set
    // pair::second is set to true if a new element was inserted, or false if an equivalent element already existed
    pair<set<pathData>::iterator, bool> ret;
    
    // generate random paths one by one, and add them to the appropriate container in our pathContainer
    for(i=0; i<nFwdPaths; i++) {
#ifndef RUN_ON_CLUSTER
        len = RandomPoisson(z*t);
#else
        len = RandomPoissonDummy(z*t);
#endif
        
        pD.random(len);
        endPointLatestPath = pD.getEndpoint();
        coefficientLatestPath = pD.getCoefficient();
        
        // insert the new pathData element, and if it was already there (ret.second==false), then just update its coefficient
        // in this case, the iterator ret.first points at the pathData element that was already there
        ret = (this->container[endPointLatestPath]).insert(pD);
        if(!ret.second) {
            (*(ret.first)).addToCoefficient(coefficientLatestPath);
        }
        
        // Display '.' whenever 1% of the paths have been finished
        if((i+1) % onePercent == 0) {
            cout << "." << flush;
        }
    }
    cout << endl;
} // Tested and works


// clear all paths from the container
void pathContainer::clear() {
    for(int i=0; i<this->container.size(); i++) {
        this->container[i].clear();
    }
    this->nPaths = 0;
    this->tHole = 0.;
} // Tested and works


// get the number of bins (i.e. the number of non-identical paths)
int pathContainer::getNumberOfBins() const {
    int siteIndex;
    set<pathData>::iterator it;
    int nBins;
    
    for(nBins=0, siteIndex=0; siteIndex<container.size(); siteIndex++) {
        for(it=(this->container[siteIndex]).begin(); it!=(this->container[siteIndex]).end(); ++it) {
            nBins++;
        }
    }
    
    return nBins;
} // Tested and works


// get memory allocation in bytes
long int pathContainer::getMemoryBytes() const {
    int siteIndex;
    set<pathData>::iterator it;
    long int mem;
    
    for(mem=0, siteIndex=0; siteIndex<container.size(); siteIndex++) {
        for(it=(this->container[siteIndex]).begin(); it!=(this->container[siteIndex]).end(); ++it) {
            mem += (*it).memoryBytes();
        }
    }
    
    return mem;
} // Tested and works


// prints bin corresponding to site 'siteIndex' to 'out'
void pathContainer::printBin(int siteIndex) const {
    set<pathData>::iterator it;
    for(it=(this->container[siteIndex]).begin(); it!=(this->container[siteIndex]).end(); ++it) {
        cout << (*it) << endl;
    }
} // Tested and works


// print all bins that are not empty
std::ostream& operator<< (std::ostream &out, const pathContainer &cont) {
    int nSites = 4 * LATTICE_SIZE * LATTICE_SIZE;
    set<pathData>::iterator it;
    for(int i=0; i<nSites; i++) {
        if(!(cont.container[i]).empty()) {
            out << "****** Bin  " << i << " = (" << indexToX(i) << ", " << indexToY(i) << ") ******" << endl;
            for(it=(cont.container[i]).begin(); it!=(cont.container[i]).end(); ++it) {
                out << (*it) << endl;
            }
            out << endl << endl << endl;
        }
    }
    return out;
} // Tested and works


