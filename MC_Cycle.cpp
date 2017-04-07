//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_Cycle.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/2/16.
//
//
// Class of a single permutation cycle
// It stores the elements of the permutation cycle starting from the smallest one
// It introduces an ordering between permutation cycles, that will be important for storing them in binary trees.


#include "MC_Cycle.hpp"

// CONSTRUCTORS AND DESTRUCTOR
// Creates a cycle without any elements
Cycle::Cycle(void) {
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    // this->check();
} // Tested and works


// Copies the cycle in 'el', rotating the vector such that the cycle starts with the smallest element of 'el'
Cycle::Cycle(const vector<int> & el) {
    // If 'el' only has a single element, then it encodes the trivial cycle
    // Otherwise we go ahead
    if(el.size() >= 2) {
    
        Elements.resize(el.size());
        int i, iSmallest, smallest;
        
        // Find the smallest element in el
        iSmallest = 0;
        smallest = el[0];
        for(i=1; i<el.size(); i++) {
            if(el[i]<smallest) {
                smallest = el[i];
                iSmallest = i;
            }
        }
        
        // Copy elements from el to Elements cyclically, such that the smallest is the first
        for(i=0; i<el.size(); i++) {
            Elements[i] =  el[(i+iSmallest) % Elements.size()];
        }
    }
    
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    /*
     // Check if the cycle was correclty constructed
     // If not, make it the trivial cycle
     if(!this->check()) {
         this->Elements.resize(0);
     }
    */
} // Tested and works

// Destructor
Cycle::~Cycle(void) {

} // Tested and works


// CHECK IF THE ELEMENTS FULFILL THE REQUIREMENTS WE NEED
// Check if 'Elements' indeed starts with its smallest element, and if it has any duplicates
bool Cycle::check(void) const {
    int i, j;
    // If Elements has a single element, then the cycle should be identically zero
    if(Elements.size()==1) {
        cerr << "Error in Cycle::check. Cycle of a single element should be zero." << endl;
        return false;
    }
    
    // Check if the cycle starts with the smallest element
    for(i=1; i<Elements.size(); i++) {
        if(Elements[i]<Elements[0]) {
            cerr << "Error in Cycle::check. Cycle does not start with the smallest element. " << endl;
            return false;
        }
    }
    
    // Check if Elements[] has any duplicates
    for(i=1; i<Elements.size(); i++) {
        for(j=0; j<i; j++) {
            if(Elements[i] == Elements[j]) {
                cerr << "Error in Cycle::check. Multiple copies of " << Elements[i] << "." << endl;
                return false;
            }
        }
    }
    
    return true;
} // Tested and works



// GET BASIC DATA
int Cycle::size(void) const {
    return Elements.size();
}

vector<int> Cycle::elements(void) const {
    return Elements;
} // Tested and works


// PRINTING
std::ostream& operator<< (std::ostream &out, const Cycle &cyc) {
    int i, len;
    len = cyc.Elements.size();
    
    //out << "(";
    out << "{";
    for(i=0; i<len; i++) {
        out << (cyc.Elements)[i];
        if(i<(len-1))
            out << ", ";
            // out << " -> ";
    }
    //out << ")";
    out << "}";
    return out;
} // Tested and works


// CIRCUMVENT MEMORY WASTE
// Requests the compiler to use only as much memory space as the size of the cycle
// This way we can avoid excess memory usage
void Cycle::shrink_to_fit(void) {
    Elements.shrink_to_fit();
}

// Estimates memory used by our cycle in the optimal case
size_t Cycle::memoryBytes(void) const {
    return Elements.capacity() * sizeof(int);
}



// OVERLOADED OPERATORS
Cycle& Cycle::operator= (const Cycle &cyc) {
    this->Elements.resize(cyc.Elements.size());
    for(int i=0; i<cyc.Elements.size(); i++) {
        (this->Elements)[i] = cyc.Elements[i];
    }
    
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    // this->check();
    return *this;
} // Tested and works


Cycle& Cycle::operator= (const vector<int> &el) {
    // If 'el' only has a single element, then it encodes the trivial cycle
    // Otherwise we go ahead
    if(el.size() < 2) {
        this->Elements.resize(0);
    } else {
        Elements.resize(el.size());
        int i, iSmallest, smallest;
        
        // Find the smallest element in el
        iSmallest = 0;
        smallest = el[0];
        for(i=1; i<el.size(); i++) {
            if(el[i]<smallest) {
                smallest = el[i];
                iSmallest = i;
            }
        }
        
        // Copy elements from el to Elements cyclically, such that the smallest is the first
        for(i=0; i<el.size(); i++) {
            Elements[i] =  el[(i+iSmallest) % Elements.size()];
        }
    }
    
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    /*
     // Check if the cycle was correclty constructed
     // If not, make it the trivial cycle
     if(!this->check()) {
         this->Elements.resize(0);
     }
    */
    
    return *this;
} // Tested and works

// Routine introducing the ordering between different cycles
// Returns 1 if cyc1<cyc2, returns -1 if cyc1>cyc2, and returns 0 if they are equal
// Condition of being smaller:
// 1. whoever has less elements
// 2. if they have the same number of elements, we go over their elements, and compare
// them one-by-one until one of them will be smaller.
// All other ordering relations are based on the output of Cycles::smaller
int smallerCycle(const Cycle &cyc1, const Cycle &cyc2) {
    int len1=cyc1.Elements.size();
    int len2=cyc2.Elements.size();
    
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
        if(cyc1.Elements[i] != cyc2.Elements[i]) {
            return (cyc1.Elements[i]<cyc2.Elements[i] ? 1 : -1);
        }
    }
    
    // If all their elements were equal, we return zero
    return 0;
} // Tested and works


// INVERSION
Cycle & Cycle::invert(void) {
    int tmp;
    int len = this->Elements.size();
    
    for(int i=0; i<(len-1)/2; i++) {
        tmp = Elements[i+1];
        Elements[i+1] = Elements[len-i-1];
        Elements[len-i-1] = tmp;
    }
    
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    /* if(!(this->check()))
        cerr << "Error in Cycle::invert. The resulting cycle does not fulfill basic properties required by Cycles." << endl;
    */
    
    return *this;
} // Tested and works


// Routine to check if two cycles have a common element
bool hasCommonElement(const Cycle &cyc1, const Cycle &cyc2) {
    int i1, i2;
    for(i1=0; i1<cyc1.Elements.size(); i1++) {
        for(i2=0; i2<cyc2.Elements.size(); i2++) {
            if(cyc1.Elements[i1] == cyc2.Elements[i2]) {
                return true;
            }
        }
    }
    return false;
} // Tested and works


