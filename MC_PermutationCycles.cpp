//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_PermutationCycles.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/4/16.
//

#include "MC_PermutationCycles.hpp"


// PROTECTED MEMBER FUNCTIONS, FOR WITHIN THE CLASS USE
// orders the elements of allCycles according to the ordering defined in MC_Cycle.hpp
void PermutationCycles::sort(void) {
    std::sort(allCycles.begin(), allCycles.end());
} // Tested and works

// checks if the cycles in 'allCycles' are in increasing order
bool PermutationCycles::isOrdered(void) const {
    if(allCycles.size() < 2)
        return true;
    for(int i=1; i<allCycles.size(); i++) {
        if(allCycles[i-1] >= allCycles[i]) {
            cerr << "Error in PermutationCycles::isOrdered. The permutation cycles are not ordered." << endl;
            return false;
        }
    }
    return true;
} // Tested and works

// checks if there are repeated elements within the cycles. Returns true if there are none, false otherwise.
bool PermutationCycles::noDuplicates(void) const {
    int i1,i2;
    for(i1=1; i1<allCycles.size(); i1++) {
        for(i2=0; i2<i1; i2++) {
            if(hasCommonElement(allCycles[i1], allCycles[i2])) {
                cerr << "Error in PermutationCycles::noDuplicates. The following cycles have common elements:" << endl;
                cerr << allCycles[i1] << endl;
                cerr << allCycles[i2] << endl;
                return false;
            }
        }
    }
    return true;
} // Tested and works

// check if there are trivial cycles in allCycles. Returns true if there is none, false otherwise.
bool PermutationCycles::noTrivialCycles(void) const {
    for(int i=0; i<allCycles.size(); i++) {
        if((allCycles[i]).size() < 2) {
            cerr << "Error in PermutationCycles::noTrivialCycles. The " << i << "th cycle is trivial." << endl;
            return false;
        }
    }
    return true;
} // Tested and works



// CONSTRUCTORS AND DESTRUCTOR
// creates an empty set of cycles
PermutationCycles::PermutationCycles(void) {
    this->allCycles.resize(0);
} // Tested and works

// creates the set of permutation cycles out of those in 'cycs'
// Makes sure that these are ordered, and there are no duplicate elements between them
// Furthermore, it throws away trivial cycles
PermutationCycles::PermutationCycles(const vector<Cycle> &cycs) {
    int i1, i2;     // indices
    bool noDup;     // whether there are duplicate elements in cycs
    int nNonTriv;   // how many non-trivial cycles are there
    
    // make sure there are no duplicate elements of the cycles contained in cycs
    // otherwise, we return false, and *this will contain no cycles
    noDup = true;
    for(i1=1; i1<cycs.size(); i1++) {
        for(i2=0; i2<i1; i2++) {
            if(hasCommonElement(cycs[i1], cycs[i2])) {
                noDup = false;
                cerr << "Error in PermutationCycles::PermutationCycles(): the following cycles have duplicate elements: " << endl;
                cerr << cycs[i1] << endl;
                cerr << cycs[i2] << endl;
            }
        }
    }
    
    // if there were no duplicates, we can proceed further
    if(noDup) {
        // check how many non-trivial elements there are in cycs
        for(i1=0, nNonTriv=0; i1<cycs.size(); i1++) {
            if(cycs[i1].size()>1)
                nNonTriv++;
        }
        
        // then, we copy the non-trivial cycles into 'allCycles'
        allCycles.resize(nNonTriv);
        for(i1=0, nNonTriv=0; i1<cycs.size(); i1++) {
            if(cycs[i1].size()>1) {
                allCycles[nNonTriv] = cycs[i1];
                nNonTriv++;
            }
        }
        
        // sort the cycles in 'allCycles'
        this->sort();
    }
    
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    /*
     // Check if the constructor worked properly
     this->check();
    */
    
} // Tested and works

// copy constructor
// Does not check if the elements of permcycs are ordered, and if there are
// multiple copies. permcycs is supposed to fulfill these requirements.
PermutationCycles::PermutationCycles(const PermutationCycles &permcycs) {
    this->allCycles.resize(permcycs.allCycles.size());
    for(int i=0; i<permcycs.allCycles.size(); i++) {
        this->allCycles[i] = permcycs.allCycles[i];
    }
    
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    /*
     this->check();
    */
    
} // Tested and works

// destructor
PermutationCycles::~PermutationCycles(void) {
     
} // Tested and works



// ADD EXTRA CYCLES
// adds the cycle to the list of permutation cycles, making sure that cyc
// contains no elements already present in 'allCycles'. Furthermore, it ensures
// that allCycles remains ordered
// Returns true if there were no common elements between cyc and the elements of *this
bool PermutationCycles::push_back(const Cycle &cyc) {
    int i1; // index
    
    // if 'cyc' is a trivial cycle, we do nothing
    if(cyc.size()<2)
        return true;
    
    // if there are elements in cyc that also appear in the cycles within 'allCycles',
    // then we cannot add 'cyc' to it, and we return false
    for(i1=0; i1<allCycles.size(); i1++) {
        if(hasCommonElement(allCycles[i1], cyc)) {
            cerr << "Error in PermutationCycles::push_back. The new cycle has common elements with (*this)." << endl;
            return false;
        }
    }
    
    // if there are no common elements, that we add cyc by keeping allCycles ordered
    int pos=0; // the position where we will insert cyc in allCycles
    for(i1=0; i1<allCycles.size(); i1++) {
        if(allCycles[pos]<cyc) {
            pos++;
        } else {
            break;
        }
    }
    allCycles.insert(allCycles.begin() + pos, cyc);
    
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    /*
     this->check();
    */
    
    return true;
} // Tested and works


// CHECK DATA
// checks if allCycles indeed satisfies what it should
// (ordered, no duplicates, no trivial cycles)
bool PermutationCycles::check(void) const {
    if(!noTrivialCycles()) {
        return false;
    }
    if(!noDuplicates()) {
        return false;
    }
    if(!isOrdered()) {
        return false;
    }
    return true;
} // Tested and works


// PRINTING
// prints all the cycles to the screen
std::ostream& operator<< (std::ostream &out, const PermutationCycles &permcycs) {
    int i, len;
    len = permcycs.allCycles.size();
    
    out << "{";
    for(i=0; i<len; i++) {
        out << (permcycs.allCycles)[i];
        if(i<(len-1))
            out << ",  ";
    }
    out << "}";
    return out;
} // Tested and works


// ASSIGNMENT OPERATOR
// Does not check if the elements of permcycs are ordered, and if there are
// multiple copies. permcycs is supposed to fulfill these requirements.
PermutationCycles& PermutationCycles::operator= (const PermutationCycles& permcycs) {
    this->allCycles = permcycs.allCycles;
    return *this;
} // Tested and works


// CIRCUMVENT MEMORY WASTE
// Requests the compiler to use only as much memory space as the size of all cycles
// This way we can avoid excess memory usage
void PermutationCycles::shrink_to_fit(void) {
    allCycles.shrink_to_fit();
    for(int i=0; i<allCycles.size(); i++) {
        allCycles[i].shrink_to_fit();
    }
} // Tested and works

// Estimates memory used by our permutation cycle
// Neglects additional capacity of allCycles that takes on extra memory space
// One can get rid of this extra capacity by invoking shrink_to_fit
size_t PermutationCycles::memoryBytes(void) const {
    int i;
    size_t mem;
    for(i=0, mem=0; i<allCycles.size(); i++) {
        mem += allCycles[i].memoryBytes();
    }
    return mem;
} // Tested and works

// MULTIPLICATION
// (p1*p2) corresponds to the permutation x -> p1(p2(x)). I.e. p2 is evaluated first, then p1.
PermutationCycles operator* (const PermutationCycles & permcyc1, const PermutationCycles & permcyc2) {
    PermutationCycles result;
    
    Permutation p1, p2;
    p1 = ToPermutation(permcyc1);
    p2 = ToPermutation(permcyc2);
    
    result=ToCycles(p1*p2);
    return result;
} // Tested and works

// INVERSION
PermutationCycles& PermutationCycles::invert(void) {
    int i1;
    for(i1=0; i1<allCycles.size(); i1++) {
        allCycles[i1].invert();
    }
    return *this;
} // Tested and works

// EFFECT ON A SET OF ELEMENTS
// act by the permutation (*this) on the elements of the vector 'el' 
// during the calculation, we assume that the elements of 'el' are likely to follow each other within the cycle
// Appologies for the hard to read code. I had to write this routine fast, even at the cost of readibility
void PermutationCycles::permute(vector<int> & el) {
    int iEl;                // iterator of the elements in 'el'
    int currentElement;     // the element we are currently looking for
    int i1, i2, i1mod;      // indices
    int iLastCycle;         // index of the cycle where we found the last element
    int iLastElement;       // index of the last element within the last cycle where we found 'currentElement'
    int imax;               // temporary index variable
    int nCycles;            // number of cycles in allCycles
    bool foundIt;           // whether we found the current element
    
    nCycles = this->allCycles.size();
    
    // If the permutation cycle is empty, then we don't have to do anything to 'el'
    if(nCycles == 0)
        return;
    
    // Since the element in 'el' are supposed to be close to each other, we start searching in the same cycle where we
    // found the element that we were previously looking for
    iLastCycle=0;
    iLastElement=0;
    for(iEl=0; iEl<el.size(); iEl++) {
        currentElement = el[iEl]; // look for this element in the Permutation cycle
        foundIt = false;
        
        // Look for currentElement
        // if you find it, then set el[iEl] to its permuted version
        // if you do not find it, then leave el[iEl] as it is
        
        // First, let us start searching at exactly the same element where we stopped last time
        for(i2=iLastElement+1, imax=iLastElement+1+allCycles[iLastCycle].size(); i2<imax; i2++) {
            if(allCycles[iLastCycle][i2 % allCycles[iLastCycle].size()] == currentElement) {
                foundIt = true;
                iLastElement = i2 % allCycles[iLastCycle].size();
                el[iEl] = allCycles[iLastCycle][(i2+1) % allCycles[iLastCycle].size()];
                break;
            }
        }
        
        // Then, go over all the other cycles
        if(!foundIt) {
            for(i1=iLastCycle+1, imax=iLastCycle+nCycles; i1<imax; i1++) {
                i1mod = i1 % nCycles; // i1 can go out of range, so we introduce this other variable as index
                for(i2=0; i2<allCycles[i1mod].size(); i2++) {
                    if(allCycles[i1mod][i2] == currentElement) {
                        foundIt = true;
                        iLastElement=i2;
                        el[iEl] = allCycles[i1mod][(i2+1) % allCycles[i1mod].size()];
                        break;
                    }
                }
                if(foundIt) {
                    iLastCycle = i1mod;
                    break;
                }
            }
        }
    }
} // Tested and works



// ORDERING
// Routine introducing the ordering among PermutationCycles
// Returns 1 if permcyc1<permcyc2, returns -1 if permcyc1>permcyc2, and returns 0 if they are equal
// Condition of being smaller:
// 1. whoever has fewer number of cycles
// 2. if they have the same number of cycles, we go over their cycles, and compare
// them one-by-one until one of them will be smaller.
int smallerPermutationCycles(const PermutationCycles &permcyc1, const PermutationCycles &permcyc2) {
    int len1=permcyc1.size();
    int len2=permcyc2.size();
    int sm;
    
    // Check if any of them has a smaller number of elements
    if(len1<len2) {
        return 1;
    } else {
        if(len2<len1) {
            return -1;
        }
    }
    
    // If they have equal length, compare their elements
    for(int i=0; i<len1; i++) {
        sm = smallerCycle(permcyc1.allCycles[i],permcyc2.allCycles[i]);
        // if (permcyc1.allCycles[i] < permcyc2.allCycles[i]) --> sm = 1
        // if (permcyc1.allCycles[i] > permcyc2.allCycles[i]) --> sm =-1
        // if (permcyc1.allCycles[i] ==permcyc2.allCycles[i]) --> sm = 0
        if(sm != 0) {
            return sm;
        }
    }
    
    // If all their elements were equal, we return zero
    return 0;

} // Tested and works





