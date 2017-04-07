//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_Permutation.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/5/16.
//
//
//  Class corresponding to the elements of the permutation group
//  Elements of the group are stored as 2xn matrices, where n is the number of elements shuffled by the permutation
//  E.g. (0, 1, 2, 3, 4, 5, 6)
//       (3, 4, 1, 5, 2, 6, 1)
//  Rules:
//  - each element can appear only once in each row
//  - there has to be the same set of elements in both rows
//  - there are no trivial cycles: an element is not allowed to be in the same column in both rows

#ifndef MC_Permutation_hpp
#define MC_Permutation_hpp

#include <stdio.h>
#include <vector>
#include <iostream>
#include "MC_clusterRunQ.hpp"
#include "MC_PermutationCycles.hpp"

using namespace std;

// Forward declare class PermutationCycles, used in a friend function
class PermutationCycles;

class Permutation{
protected:
    vector< vector<int> > perm;            // 2xn matrix containing the permutation
    
    // Routines to make sure that perm satisfies the conditions we want from a permutation
    bool sizeOK(void) const;               // Checks if perm is indeed a 2xn matrix
    bool elementsOK(void) const;           // Checks if the first and second row have the same set of elements, and each element appears only once
    bool noTrivialCycles(void) const;      // Checks if there are single element (trivial) permutation cycles
    void removeTrivialCycles(void);        // Removes single element (trivial) permutation cycles
    
public:
    // CONSTRUCTORS AND DESTRUCTOR
    Permutation(void);                          // trivial constructor
    Permutation(vector< vector<int> >);         // copies permutation from the 2xn matrix to *this
    Permutation(const Permutation & p);         // copy constructor
    ~Permutation(void);                         // destructor
    
    // CHECK IF THE PERMUTATION AT HAND SATISFIES THE RULES
    bool check(void) const;
    
    // GET BASIC DATA
    vector<int> getElements(void) const;                // get the first row of permutations, containing all elements
    vector< vector<int> > getPermutation(void) const;   // get the full permutation representation perm, as a 2xn matrix
    int size(void) const;                               // number of elements in the permutation
    
    // PRINTING
    friend std::ostream& operator<< (std::ostream &out, const Permutation &p);
    
    // ASSIGNMENT OPERATOR
    Permutation& operator= (const Permutation & p);
    
    // CONVERT TO PERMUTATION CYCLES
    friend PermutationCycles ToCycles(Permutation p);   // get the decomposition of the permutation into cycles
                                                        // the function is deliberately not called with a const handle
                                                        // we need to create a copy of p within the function, since we perform
                                                        // calculations in that memory space
    
    // CONVERT FROM PERMUTATION CYCLES
    // Converts the permutation cycles permcyc into a permutation. It is defined in MC_permutation.cpp
    friend Permutation ToPermutation(const PermutationCycles &permcyc);
    
    // MULTIPLICATION
    friend Permutation operator* (const Permutation & p1, const Permutation & p2);
    
    // INVERSION
    Permutation& invert(void);
    
protected:    
    // FIND ELEMENT
    // Looks for 'element' in the first row of 'perm'
    // Returns true if it found it, and puts its index into i
    // Returns false if it didn't find it
    bool findelement(const int & element, int & i) const;
};


#endif /* MC_Permutation_hpp */
