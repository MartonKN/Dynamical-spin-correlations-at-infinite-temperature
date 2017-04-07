//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_singleHoleQState.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/19/16.
//
//

#include "MC_singleHoleQState.hpp"




// ---------------------------------------------------
// QUANTUM STATE OF THE SPIN SYSTEM WITH A SINGLE HOLE
// ---------------------------------------------------
// constructor: fill QState with random numbers and put a hole in the origin
singleHoleQState::singleHoleQState(void) {
    this->QState.reserve(nSites);
    this->QState.resize(nSites);
    this->QState.shrink_to_fit();
    
    for(int i=0; i<nSites; i++) {
        this->QState[i] = (signed char) (RandomInteger(0, 1) == 1 ? 1 : -1);
    }
    this->holePosition = coordToIndex(0,0);
    this->QState[this->holePosition] = (signed char) 0;
    
    this->coefficient.real(1.);
    this->coefficient.imag(0.);
} // Tested and works

// destructor
singleHoleQState::~singleHoleQState(void) {
    
} // Tested and works

// generate Neél order and a hole at the origin
void singleHoleQState::Neel(void) {
    int x,y,iSite;
    for(iSite=0; iSite<nSites; iSite++) {
        x = indexToX(iSite);
        y = indexToY(iSite);
        this->QState[iSite] = (signed char) ((abs(x+y) % 2) == 0 ? 1 : -1);
    }
    this->holePosition = coordToIndex(0,0);
    this->QState[this->holePosition] = (signed char) 0;
    
    this->coefficient.real(1.);
    this->coefficient.imag(0.);
} // Tested and works

// generate ferromagnetic spin configuration and a hole at the origin
void singleHoleQState::Ferromagnet(void) {
    for(int iSite=0; iSite<nSites; iSite++)
        this->QState[iSite] = (signed char) 1;
    this->holePosition = coordToIndex(0,0);
    this->QState[this->holePosition] = (signed char) 0;
    
    this->coefficient.real(1.);
    this->coefficient.imag(0.);
}; // Tested and works


// generate random spin configuration and a hole at the origin, set coefficient to 1.
void singleHoleQState::random() {
    for(int iSite = 0; iSite<nSites; iSite++) {
        this->QState[iSite] = (signed char) (RandomInteger(0,1) == 1 ? 1 : -1);
    }
    this->holePosition = coordToIndex(0,0);
    this->QState[this->holePosition] = (signed char) 0;
    
    this->coefficient.real(1.);
    this->coefficient.imag(0.);
} // Tested and works

// reset quantum state to the input values
bool singleHoleQState::reset(const vector<int> &qS, complex<double> coeff) {
    int i;
    int nHoles;
    int holeSite;
    if(qS.size() != nSites) {
        cerr << "Error in singleHoleQState::reset. The input vector is of inappropriate size. The quantum state has not been updated." << endl;
        return false;
    }
    for(nHoles=0, i=0; i<nSites; i++) {
        if(qS[i]!=-1 && qS[i]!=0 && qS[i]!=1) {
            cerr << "Error in singleHoleQState::reset. The quantum state can only have 1, -1 and 0 elements. ";
            cerr << "The quantum state has not been updated." << endl;
            return false;
        } else if(qS[i]==0) {
            nHoles ++;
            holeSite = i;
        }
    }
    if(nHoles!=1) {
        cerr << "Error in singleHoleQState::reset. The input vector should have a single hole. The quantum state has not been updated." << endl;
        return false;
    }
    
    // If everything is OK, we reset the qantum state, the hole position and the coefficient
    this->coefficient = coeff;
    this->holePosition = holeSite;
    for(i=0; i<nSites; i++) {
        this->QState[i] = (signed char) qS[i];
    }
    return true;
} // Tested and works

// let the permutation of sites 'pC' act on the spins of the state
void singleHoleQState::permute(const PermutationCycles & pC) {
    int i1, i2;
    signed char tmp;
    int site;
    int newHolePosition = this->holePosition;
    
    // Loop over all cycles of pC, and move the spins around accordingly
    for(i1=0; i1<pC.size(); i1++) {
        site = pC[i1][pC[i1].size() - 1]; // last element of the current cycle
        tmp = this->QState[site];         // save its value
        if(this->holePosition == site) {
            newHolePosition = pC[i1][0];
        }
        
        for(i2=(pC[i1].size()-1); i2>0; i2--) {
            this->QState[ pC[i1][i2] ] = this->QState[ pC[i1][i2-1] ];
            if(this->holePosition == pC[i1][i2-1]) {
                newHolePosition = pC[i1][i2];
            }
        }
        this->QState[ pC[i1][0] ] = tmp;
    }
    this->holePosition = newHolePosition;
    
} // Tested and works


// check if the state indeed has a single hole in it and only spins 1 and -1
bool singleHoleQState::check(void) const {
    int i;          // iterator
    int nHoles;     // number of holes
    
    for(nHoles=0, i=0; i<nSites; i++) {
        if ((this->QState[i] != -1) && (this->QState[i] != 0) && (this->QState[i] != 1)) {
            cerr << "Error in singleHoleQState::check(): ";
            cerr << "the quantum state can only have elements up (1), down (-1) or hole(0)." << endl;
            return false;
        }
        if(this->QState[i] == 0) {
            nHoles ++;
        }
    }
    if(nHoles!=1) {
        cerr << "Error in singleHoleQState::check(): ";
        cerr << "there should only be a single hole site." << endl;
        return false;
    }
    
    return true;
} // Tested and works

// Assignment operators
singleHoleQState& singleHoleQState::operator= (const singleHoleQState & qS) {
    for(int i=0; i<nSites; i++) {
        this->QState[i] = qS.QState[i];
    }
    this->coefficient = qS.coefficient;
    this->holePosition = qS.holePosition;
    return *this;
} // Tested and works

// Printing
ostream & operator<<(ostream& out, const singleHoleQState & qS) {
    int x,y;
    out << "Amplitude: (" << qS.coefficient.real() << ", " << qS.coefficient.imag() << ")" << endl;
    out << "Hole position: " << qS.holePosition << " = (" << indexToX(qS.holePosition) << ", " << indexToY(qS.holePosition) << ")" << endl;
    out << "State: " << endl;
    
    for(x=-LATTICE_SIZE; x<LATTICE_SIZE; x++) {
        for(y=-LATTICE_SIZE; y<LATTICE_SIZE; y++) {
            out << qS.QState[coordToIndex(x,y)] << ",\t";
        }
        out << endl;
    }
    return out;
} // Tested and works


// Ordering simply based on the ordering of the elements, starting at the origin = coordToIndex(0,0)
// Returns 1 if qS1<qS2, returns -1 if qS1>qS2 and returns 0 if qS1==qS2
int smallerSingleHoleQState(const singleHoleQState& qS1, const singleHoleQState& qS2) {
    int i;
    int origin = coordToIndex(0,0);
    for(i=origin; i<singleHoleQState::nSites; i++) {
        if(qS1[i]<qS2[i]) {
            return 1;
        } else if (qS1[i]>qS2[i]) {
            return -1;
        }
    }
    for(i=origin-1; i>=0; i--) {
        if(qS1[i]<qS2[i]) {
            return 1;
        } else if (qS1[i]>qS2[i]) {
            return -1;
        }
    }
    return 0;
} // Tested and works


