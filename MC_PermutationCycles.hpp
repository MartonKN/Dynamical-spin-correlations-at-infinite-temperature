//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_PermutationCycles.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/4/16.
//
//
//  Class describing the cycles representation of a permutation
//  Its elements are single cycles, defined in MC_cycles.hpp
//  The cycles cannot have common elements

#ifndef MC_PermutationCycles_hpp
#define MC_PermutationCycles_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "MC_clusterRunQ.hpp"
#include "MC_Cycle.hpp"
#include "MC_Permutation.hpp"

using namespace std;


// Forward declare class Permutation, used in a friend function
class Permutation;

class PermutationCycles {
protected:
    vector<Cycle> allCycles;
    // Ordered vector of single permutation cycles, allCycles = {cyc0, cyc1, cyc2, ...}, such that cyc0 < cyc1 < cyc2 < ...,
    // using the ordering reliations introduced in the Cycles class in MC_Cycles.hpp. This ordering is important, because it will
    // allow us to introduce an ordering relation between PermutationCycles. The ordering within the class PermutationCycles will
    // make it possible for us to store these in a binary search tree during the simulation.
    // Furthermore, we require that the cycles have no common elements
    // If one of the cycles happens to be trivial (it has zero or one elements), we simply ommit it
    
    // Auxiliary routines used by the member functions:
    void sort(void);                   // orders the elements of allCycles according to the ordering defined in MC_Cycle.hpp
    bool isOrdered(void) const;        // checks if the cycles in 'allCycles' are in increasing order
    bool noDuplicates(void) const;     // checks if there are repeated elements within the cycles. Returns true if there are none, false otherwise.
    bool noTrivialCycles(void)const ;  // check if there are trivial cycles in allCycles. Returns true if there is none, false otherwise.
    
public:
    // CONSTRUCTORS AND DESTRUCTOR
    PermutationCycles(void);                        // creates an empty set of cycles
    PermutationCycles(const vector<Cycle> &cycs);   // creates the set of permutation cycles out of those in 'cycs'
                                                    // Makes sure that these are ordered, and there are no duplicate elements between them
                                                    // Furthermore, it throws away trivial cycles
    PermutationCycles(const PermutationCycles &permcycs); // copy constructor
                                                          // Does not check if the elements of permcycs are ordered, and if there are
                                                          // multiple copies. permcycs is supposed to fulfill these requirements.
    ~PermutationCycles(void);                       // destructor
    
    // ADD EXTRA CYCLES
    bool push_back(const Cycle &cyc);               // adds the cycle to the list of permutation cycles, making sure that cyc
                                                    // contains no elements already present in 'allCycles'. Furthermore, it ensures
                                                    // that allCycles remains ordered
                                                    // Returns true if there were no common elements between cyc and the elements of *this
    
    // GET BASIC DATA
    inline int size(void) const {return this->allCycles.size();} // get number of cycles
    inline vector<int> getCycleLengths(void) const  // get the number of elements in each cycle within the permutation cycle
    {
        vector<int> cycLens(allCycles.size());
        for(int i=0; i<allCycles.size(); i++) {
            cycLens[i] = (allCycles[i]).size();
        }
        return cycLens;
    }
    
    // CHECK DATA
    bool check(void) const;                         // checks if allCycles indeed satisfies what it should
                                                    // (ordered, no duplicates, no trivial cycles)
    
    // PRINTING
    friend std::ostream& operator<< (std::ostream &out, const PermutationCycles &permcycs); // prints all the cycles to the screen
    
    // ASSIGNMENT OPERATOR
    PermutationCycles& operator= (const PermutationCycles& permcycs);   // assignment operator.
                                                                        // Does not check permcycs, assumes it satisfies all properties
    inline Cycle operator[] (const int &i) const                        // get element (read only)
    {
        return (this->allCycles)[i];
    }
    
    // CIRCUMVENT MEMORY WASTE
    void shrink_to_fit(void);           // Requests the compiler to use only as much memory space as the size of all cycles
                                        // This way we can avoid excess memory usage
    size_t memoryBytes(void) const;     // Estimates memory used by our permutation cycle
                                        // Neglects additional capacity of allCycles that takes on extra memory space
                                        // One can get rid of this extra capacity by invoking shrink_to_fit

    // CONVERT FROM PERMUTATION CYCLES TO PERMUTATION
    // Converts the permutation cycles permcyc into a permutation. It is defined in MC_permutation.cpp
    friend Permutation ToPermutation(const PermutationCycles &permcyc);
    
    // MULTIPLICATION
    // (perm1*perm2)(x) is the permutation that maps x -> perm1(perm2(x))
    friend PermutationCycles operator*(const PermutationCycles & permcyc1, const PermutationCycles & permcyc2);
    
    // INVERSION
    PermutationCycles& invert(void);
    
    // EFFECT ON A SET OF ELEMENTS
    // act by the permutation (*this) on the elements of the vector 'el'
    // during the calculation, we assume that the elements of 'el' are likely to follow each other within the cycle
    void permute(vector<int> & el);
    
    // MULTIPLICATION
    // (p1*p2) corresponds to the permutation x -> p1(p2(x)). I.e. p2 is evaluated first, then p1.
    friend PermutationCycles operator* (const PermutationCycles & permcyc1, const PermutationCycles & permcyc2);
    
    // ORDERING
    // Routine introducing the ordering among PermutationCycles
    // Returns 1 if permcyc1<permcyc2, returns -1 if permcyc1>permcyc2, and returns 0 if they are equal
    // Condition of being smaller:
    // 1. whoever has fewer number of cycles
    // 2. if they have the same number of cycles, we go over their cycles, and compare
    // them one-by-one until one of them will be smaller.
    friend int smallerPermutationCycles(const PermutationCycles &permcyc1, const PermutationCycles &permcyc2);
    
    // All other ordering relations are based on the results of PermutationCycles::smaller
    // All of these were tested and they work
    inline friend bool operator== (const PermutationCycles &permcyc1, const PermutationCycles &permcyc2) {
        return (smallerPermutationCycles(permcyc1,permcyc2) == 0);
    }
    inline friend bool operator!= (const PermutationCycles &permcyc1, const PermutationCycles &permcyc2) {
        return !(smallerPermutationCycles(permcyc1,permcyc2) == 0);
    }
    inline friend bool operator<  (const PermutationCycles &permcyc1, const PermutationCycles &permcyc2) {
        return (smallerPermutationCycles(permcyc1,permcyc2) == 1);
    }
    inline friend bool operator<=  (const PermutationCycles &permcyc1, const PermutationCycles &permcyc2) {
        return (smallerPermutationCycles(permcyc1,permcyc2) != -1);
    }
    inline friend bool operator>  (const PermutationCycles &permcyc1, const PermutationCycles &permcyc2) {
        return (smallerPermutationCycles(permcyc1,permcyc2) == -1);
    }
    inline friend bool operator>=  (const PermutationCycles &permcyc1, const PermutationCycles &permcyc2) {
        return (smallerPermutationCycles(permcyc1,permcyc2) != 1);
    }
};

#endif /* PermutationCycles_hpp */




/* 
// TESTING WHETHER THE RELATIONS DEFINED IN THIS CLASS SATISFY THE ORDERING AXIOMS
 int i1, i2;
 int l1, l2, l3;
 int pl1, pl2, pl3;
 vector<int> d1, d2, d3;
 Cycle cyc1, cyc2, cyc3;
 PermutationCycles permcyc1, permcyc2, permcyc3, emptypermcyc;
 bool isOK;
 
 
 for(int i1=0; i1<10000; i1++) {
 permcyc1 = emptypermcyc;
 permcyc2 = emptypermcyc;
 permcyc3 = emptypermcyc;
 
 pl1 = RandomInteger(0,3);
 pl2 = RandomInteger(0,3);
 pl3 = RandomInteger(0,3);
 
 for(i2=0; i2<pl1; i2++) {
 l1 = RandomInteger(0,3);
 d1 = RandomInteger(-10,10,l1);
 cyc1 = d1;
 permcyc1.push_back(cyc1);
 }
 for(i2=0; i2<pl2; i2++) {
 l2 = RandomInteger(0,3);
 d2 = RandomInteger(-10,10,l2);
 cyc2 = d2;
 permcyc2.push_back(cyc2);
 }
 for(i2=0; i2<pl3; i2++) {
 l3 = RandomInteger(0,3);
 d3 = RandomInteger(-10,10,l3);
 cyc3 = d3;
 permcyc3.push_back(cyc3);
 }
 
 cout << endl;
 
 isOK = true;
 isOK = isOK && (permcyc1 == permcyc1);
 isOK = isOK && !(permcyc1 < permcyc1);
 isOK = isOK && !(permcyc1 > permcyc1);
 isOK = isOK && (permcyc1 >= permcyc1);
 isOK = isOK && (permcyc1 <= permcyc1);
 isOK = isOK && !(permcyc1<permcyc2 && permcyc2<permcyc1);
 if(permcyc1 < permcyc2 && permcyc2 < permcyc3) {
 isOK = isOK && (permcyc1 < permcyc3);
 }
 if(permcyc1 <= permcyc2 && permcyc2 <= permcyc3) {
 isOK = isOK && (permcyc1 <= permcyc3);
 }
 if(permcyc1 <= permcyc2 && permcyc2 <= permcyc1) {
 isOK = isOK && (permcyc1 == permcyc2);
 }
 
 
 if(!isOK) {
 cout << endl << endl;
 cout << "*************************************************************************************************************** " << endl;
 cout << "ERROR" << endl;
 cout << "permcyc1: " << permcyc1 << endl;
 cout << "permcyc2: " << permcyc2 << endl;
 cout << "permcyc3: " << permcyc3 << endl;
 cout << "*************************************************************************************************************** " << endl;
 cout << endl << endl;
 }
 }
*/


/*
 // TESTING THE MULTIPLICATION OF PERMUTATIONS
 
 // Initialize random number generator
 struct timeval tp;
 gettimeofday(&tp, NULL);
 long int time_in_microsec = tp.tv_sec * 1000000 + tp.tv_usec;
 srand (mod(time_in_microsec, RAND_MAX));
 
 
 vector<int> d;
 Cycle cyc;
 PermutationCycles permcyc1, permcyc2, permcyc0;
 
 
 int l=RandomInteger(0,40);
 int i1, i2, i;
 
 
 permcyc1=permcyc0;
 permcyc2=permcyc0;
 
 vector<int> elements(l);
 for(i1=0;i1<l; i1++)
 elements[i1] = i1+1;
 
 i1=0;
 i2 = RandomInteger(0,10);
 while(i2<l) {
 d.resize(i2-i1);
 for(i=0; i<(i2-i1); i++) {
 d[i] = elements[i+i1];
 }
 RandomShuffle(d);
 i1 = i2;
 i2 += RandomInteger(0,6);
 
 cyc = d;
 permcyc1.push_back(cyc);
 }
 
 i1=0;
 i2 = RandomInteger(0,10);
 while(i2<l) {
 d.resize(i2-i1);
 for(i=0; i<(i2-i1); i++) {
 d[i] = elements[i+i1];
 }
 RandomShuffle(d);
 i1 = i2;
 i2 += RandomInteger(0,6);
 
 cyc = d;
 permcyc2.push_back(cyc);
 }
 
 
 cout << endl << endl << endl;
 cout << "c1 = Cycles[" << permcyc1 << "];" << endl;
 cout << "c2 = Cycles[" << permcyc2 << "];" << endl;
 cout << "PermutationProduct[c2, c1] == Cycles[" << permcyc1*permcyc2 << "]" << endl;
 cout << endl << endl << endl;
 
 return 0;

 */


/*
 // Test PermutationCycles::permute
 void printVec(const vector<int> & el) {
 for(int i=0; i<el.size(); i++) {
 cout << el[i] << ", ";
 };
 cout << endl;
 }
 
 PermutationCycles randomPermutation(int len=10, int imin=1, int imax=100) {
 // len=number of elements in the permutation cycle
 // (imin, imax): range of elements
 
 vector<int> rng(imax-imin+1);
 Cycle cyc;
 int nPieces;
 vector<int> delimiters;
 vector<int> cycVec;
 PermutationCycles permCyc;
 int i,j;
 
 if(len==0) {
 return permCyc;
 }
 
 for(j=0, i=imin; i<=imax; i++,j++) {
 rng[j] = i;
 }
 RandomShuffle(rng);
 
 nPieces = RandomInteger(0,len-1);
 if(nPieces >= 2) {
 delimiters = RandomInteger(1, len-2, nPieces-1);
 }
 delimiters.push_back(0);
 delimiters.push_back(len-1);
 sort(delimiters.begin(), delimiters.end());
 cycVec.resize(delimiters[1]-delimiters[0]+1);
 for(j=0; j<cycVec.size(); j++) {
 cycVec[j] = rng[delimiters[0]+j];
 }
 cyc = cycVec;
 permCyc.push_back(cyc);
 
 for(i=2; i<delimiters.size(); i++) {
 cycVec.resize(delimiters[i]-delimiters[i-1]);
 for(j=0; j<cycVec.size(); j++) {
 cycVec[j] = rng[delimiters[i-1]+j+1];
 }
 cyc = cycVec;
 permCyc.push_back(cyc);
 }
 return permCyc;
 }
 
 void permuteTest(const PermutationCycles & permCyc, vector<int> & el) {
 int nVec;
 vector<int> tmpVec;
 vector<int>::iterator it;
 int i,j,k;
 
 for(nVec=0, i=0; i<permCyc.size(); i++) {
 nVec += (permCyc[i].size() + 1);
 }
 
 // Copy the elements of the permutation cycle into a vector
 tmpVec.resize(nVec);
 for(k=0, i=0; i<permCyc.size(); i++) {
 for(j=0; j<permCyc[i].size(); j++) {
 tmpVec[k] = permCyc[i][j];
 k++;
 }
 tmpVec[k] = permCyc[i][0];
 k++;
 }
 
 // Find elements of the permutation cycle in tmpVec, and figure out what they become after the permutation
 for(i=0; i<el.size(); i++) {
 it = find(tmpVec.begin(), tmpVec.end(), el[i]);
 if(it != tmpVec.end()) {
 ++it;
 el[i] = *it;
 }
 }
 }
 
 int main(int argc,char *argv[])
 {
 InitializeRandomNumberGenerator();
 int imin, imax;
 PermutationCycles permCyc;
 vector<int> el, elPerm1, elPerm2;
 
 for(int i=0; i<1000000; i++) {
 imin=RandomInteger(0,100);
 imax=imin+RandomInteger(0,100);
 
 permCyc = randomPermutation(RandomInteger(0,imax-imin+1),imin,imax);
 
 el = RandomInteger(imin,imax,RandomInteger(imin,imax));
 elPerm1 = el;
 elPerm2 = el;
 
 permCyc.permute(elPerm1);
 permuteTest(permCyc, elPerm2);
 
 if(elPerm1 != elPerm2) {
 printVec(el);
 printVec(elPerm1);
 printVec(elPerm2);
 cout << permCyc << endl;
 cout << endl << endl << endl;
 }
 }
 
 return 0;
 }

 
 */
