//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_singleHoleQStateBin.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/19/16.
//
//

#include "MC_singleHoleQStateBin.hpp"


// --------------------------------------------------------------------------
// CONTAINER OF SPIN STATES (PREFERABLY, BUT NOT NECESSARILY AT A GIVEN SITE)
// --------------------------------------------------------------------------
// The class 'singleHoleQStateContainer' has a similar functionality, but it turned out to consume too much memory.
// Since in many cases we just need to convert 'pathData' elements ending at a given site into 'singleHoleQState'-s,
// we need a container to temporarily store all those states.
//
// Stores quantum states in a red-black self balanced binary search tree provided by 'std::set', to ensure easy and fast access


// get memory allocation in bytes
long int singleHoleQStateBin::getMemoryBytes(void) const {
    set<singleHoleQState>::iterator it;
    long int mem;
    
    for(mem=0, it=this->bin.begin(); it!=this->bin.end(); ++it) {
        mem += (*it).getMemoryBytes();
    }
    return mem;
}


// Fill '*this' with quantum states that were permuted according to the action of the hole in pathContainer pC.
// It deletes all previous entries from '*this'
void singleHoleQStateBin::fill(const set<pathData> &sPD, const singleHoleQState &qSInitial) {
    this->clear();
    this->add(sPD, qSInitial);
}

// Does the same thing as fill, but without deleting all previous entries first
void singleHoleQStateBin::add (const set<pathData> &sPD, const singleHoleQState &qSInitial) {
    set<pathData>::iterator it;
    pathData pD;
    singleHoleQState qSTransformed;
    complex<double> coefficientTransformed;
    
    // We return a pair of values to be returned by 'std::set' when we insert a new element of singleHoleQState
    // pair::first is set to an iterator pointing to either the newly inserted element or to the equivalent element already in the set
    // pair::second is set to true if a new element was inserted, or false if an equivalent element already existed
    pair<set<singleHoleQState>::iterator, bool> ret;
    
    // We take each bin of paths from pC, transform the initial state with it, and then add the new quantum state to (*this).
    for(it=sPD.begin(); it!=sPD.end(); ++it) {
        qSTransformed = qSInitial;
        qSTransformed.permute((*it).getPathCycles());
        coefficientTransformed = (*it).getCoefficient();
        qSTransformed.setCoefficient(coefficientTransformed);
        
        // insert the new singleHoleQState element, and if it was already there (ret.second==false), then just update its coefficient
        // in this case, the iterator ret.first points at the pathData element that was already there
        ret = this->bin.insert(qSTransformed);
        if(!ret.second) {
            (*(ret.first)).addToCoefficient(coefficientTransformed);
        }
    }
}


// print the contents of the bin
ostream & operator<<(ostream &out, const singleHoleQStateBin &qSBin) {
    set<singleHoleQState>::iterator it;
    for(it=qSBin.bin.begin(); it!=qSBin.bin.end(); ++it) {
        cout << (*it) << endl;
    }
    return out;
}
