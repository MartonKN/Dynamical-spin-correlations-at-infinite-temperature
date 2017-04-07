//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_singleHoleQStateContainer.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/19/16.
//
//

#include "MC_singleHoleQStateContainer.hpp"



// ------------------------
// CONTAINER OF SPIN STATES
// ------------------------
// constructor
singleHoleQStateContainer::singleHoleQStateContainer(void) {
    InitializeRandomNumberGenerator();
    this->container.reserve(4 * LATTICE_SIZE * LATTICE_SIZE);
    this->container.resize (4 * LATTICE_SIZE * LATTICE_SIZE);
    this->container.shrink_to_fit();
    this->clear();
} // Tested and works

// destructor that clears all sets
singleHoleQStateContainer::~singleHoleQStateContainer(void) {
    this->clear();
} // Tested and works

// clear all states from the container
void singleHoleQStateContainer::clear() {
    for(int i=0; i<this->container.size(); i++) {
        this->container[i].clear();
    }
    this->container.resize(4*LATTICE_SIZE*LATTICE_SIZE);
    this->container.shrink_to_fit();
} // Tested and works

// get the number of bins (i.e. the number of non-identical paths)
int singleHoleQStateContainer::getNumberOfBins(void) const {
    int siteIndex;
    set<singleHoleQState>::iterator it;
    int nBins;
    
    for(nBins=0, siteIndex=0; siteIndex<container.size(); siteIndex++) {
        for(it=(this->container[siteIndex]).begin(); it!=(this->container[siteIndex]).end(); ++it) {
            nBins++;
        }
    }
    
    return nBins;
} // Tested and works

// get memory allocation in bytes
long int singleHoleQStateContainer::getMemoryBytes(void) const {
    int siteIndex;
    set<singleHoleQState>::iterator it;
    long int mem;
    
    for(mem=0, siteIndex=0; siteIndex<container.size(); siteIndex++) {
        for(it=(this->container[siteIndex]).begin(); it!=(this->container[siteIndex]).end(); ++it) {
            mem += (*it).getMemoryBytes();
        }
    }
    return mem;
} // Tested and works

// fills '*this' with quantum states that were permuted according to the action of the hole in pathContainer pC
// It deletes all previous entries from '*this'
void singleHoleQStateContainer::fill(const pathContainer &pC, const singleHoleQState &qSInitial) {
    this->clear();
    this->add(pC, qSInitial);
} // Looked the code through and seems to be good

// does the same thing as 'this->fill()', but without deleting all previous entries first
void singleHoleQStateContainer::add(const pathContainer &pC, const singleHoleQState &qSInitial) {
    int nSites = 4*LATTICE_SIZE*LATTICE_SIZE;
    int siteIndex;
    
    set<pathData>::iterator it;
    pathData pD;
    
    singleHoleQState qSTransformed;
    complex<double> coefficientTransformed;
    
    // We return a pair of values to be returned by 'std::set' when we insert a new element of singleHoleQState
    // pair::first is set to an iterator pointing to either the newly inserted element or to the equivalent element already in the set
    // pair::second is set to true if a new element was inserted, or false if an equivalent element already existed
    pair<set<singleHoleQState>::iterator, bool> ret;
    
    // We take each bin of paths from pC, transform the initial state with it, and then add the new quantum state to (*this).
    for(siteIndex=0; siteIndex<nSites; siteIndex++) {
        for(it=(pC[siteIndex]).begin(); it!=(pC[siteIndex]).end(); ++it) {
            qSTransformed = qSInitial;
            qSTransformed.permute((*it).getPathCycles());
            coefficientTransformed = (*it).getCoefficient();
            qSTransformed.setCoefficient(coefficientTransformed);
            
            
            // insert the new singleHoleQState element, and if it was already there (ret.second==false), then just update its coefficient
            // in this case, the iterator ret.first points at the pathData element that was already there
            ret = (this->container[siteIndex]).insert(qSTransformed);
            if(!ret.second) {
                (*(ret.first)).addToCoefficient(coefficientTransformed);
            }
        }
    }
} // Looked the code through and seems to be good

// prints the elements of the bin at site 'siteIndex'
void singleHoleQStateContainer::printBin(const int siteIndex) const {
    set<singleHoleQState>::iterator it;
    
    cout << "Elements of the bin " << siteIndex << " = (" << indexToX(siteIndex) << ", " << indexToY(siteIndex) << "):" << endl;
    for(it=this->container[siteIndex].begin(); it!=this->container[siteIndex].end(); ++it) {
        cout << (*it) << endl;
    }
} // Tested and works

// print the entire container
ostream & operator<<(ostream &out, const singleHoleQStateContainer &qSC) {
    set<singleHoleQState>::iterator it;
    for(int i=0; i<qSC.size(); i++) {
        if(!(qSC.container[i]).empty()) {
            out << "****** Bin  " << i << " = (" << indexToX(i) << ", " << indexToY(i) << ") ******" << endl;
            for(it=(qSC.container[i]).begin(); it!=(qSC.container[i]).end(); ++it) {
                out << (*it) << endl;
            }
            out << endl << endl << endl;
        }
    }
    return out;
} // Tested and works



