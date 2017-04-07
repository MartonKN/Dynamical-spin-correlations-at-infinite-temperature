//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_singleHoleQState.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/19/16.
//
//

#ifndef MC_singleHoleQState_hpp
#define MC_singleHoleQState_hpp

#include <stdio.h>
#include <cmath>
#include <complex>
#include <iostream>
#include "MC_clusterRunQ.hpp"
#include "MC_random.hpp"
#include "MC_coord.hpp"
#include "MC_PermutationCycles.hpp"

using namespace std;



// ---------------------------------------------------
// QUANTUM STATE OF THE SPIN SYSTEM WITH A SINGLE HOLE
// ---------------------------------------------------
class singleHoleQState {
private:
    static const int NUMBER_OF_SPINS = 2;                        // how many spin components there are
    static const int nSites = 4 * LATTICE_SIZE * LATTICE_SIZE;   // number of lattice sites
    vector<signed char> QState;                                // quantum state of spins up (1) and down(-1), as well as a single hole (0)
    int holePosition;                                            // position of the hole
    mutable complex<double> coefficient;                         // amplitude
    
public:
    // Constructors, destructors and methods to change the spin configuration
    // All of these set 'this->coefficient' to 1.
    singleHoleQState(void);                                 // constructor: generates random spin configuration and a spin at the origin
    ~singleHoleQState(void);                                // destructor
    void Neel(void);                                        // generate Neél order and a hole at the origin
    void Ferromagnet(void);                                 // generate ferromagnetic spin configuration and a hole at the origin
    void random();                                          // generate random spin configuration and a hole at the origin
    bool reset(const vector<int> &qS, complex<double> coeff = complex<double>(1., 0.)); // reset quantum state to the input values
    // returns true if successful, false if not
    void permute(const PermutationCycles & pC);             // let the permutation of sites 'pC' act on the spins of the state
    
    // check if the state indeed has a single hole in it and only spins 1 and -1
    bool check(void) const;
    
    // Get basic data
    inline int size(void) const {return QState.size();};                                                            // Tested and works
    inline vector<int> getQState(void) const {
        vector<int> tmpQState(QState.size());
        for(int i=0; i<QState.size(); i++) {
            tmpQState[i] = (int) (QState[i]);
        }
        return tmpQState;
    };                                                                                                              // Tested and works
    inline int getHolePosition(void) const {return holePosition;};                                                  // Tested and works
    inline complex<double> getCoefficient(void) const {return coefficient;};                                        // Tested and works
    size_t getMemoryBytes(void) const {return (QState.capacity()*sizeof(char) + 2*sizeof(double) + sizeof(int));};  // Tested and works
    
    
    // update the value of coefficient ('coefficient' is mutable, so these routines can be declared const)
    inline void addToCoefficient(const complex<double> & c) const {coefficient += c;}                               // Tested and works
    inline void setCoefficient  (const complex<double> & c) const {coefficient = c;}                                // Tested and works
    
    // Assignment operators
    singleHoleQState& operator= (const singleHoleQState & qS);                                                      // Tested and works
    
    // Element operato
    inline int operator[] (const int &i) const { return (this->QState)[i]; };                                       // Tested and works
    
    // Printing
    friend ostream & operator<<(ostream& out, const singleHoleQState & qS);                                         // Tested and works
    
    // Ordering
    // simply based on the ordering of the elements, starting at the origin = coordToIndex(0,0)
    friend int smallerSingleHoleQState(const singleHoleQState& qS1, const singleHoleQState& qS2);                   // Tested and works
    // smallerSingleHoleQState defines all the ordering relations for 'singleHoleQState'-s
    // Returns 1 if qS1<qS2, returns -1 if qS1>qS2, and returns 0 if they are equal
    
    // The ordering routines below have been tested and they work properly
    inline friend bool operator== (const singleHoleQState& qS1, const singleHoleQState& qS2) {
        return (smallerSingleHoleQState(qS1,qS2) == 0);
    };
    inline friend bool operator!= (const singleHoleQState& qS1, const singleHoleQState& qS2) {
        return !(smallerSingleHoleQState(qS1,qS2)== 0);
    };
    inline friend bool operator<  (const singleHoleQState& qS1, const singleHoleQState& qS2) {
        return (smallerSingleHoleQState(qS1,qS2) == 1);
    };
    inline friend bool operator<= (const singleHoleQState& qS1, const singleHoleQState& qS2) {
        return (smallerSingleHoleQState(qS1,qS2) !=-1);
    };
    inline friend bool operator>  (const singleHoleQState& qS1, const singleHoleQState& qS2) {
        return (smallerSingleHoleQState(qS1,qS2) ==-1);
    };
    inline friend bool operator>= (const singleHoleQState& qS1, const singleHoleQState& qS2) {
        return (smallerSingleHoleQState(qS1,qS2) != 1);
    };
    
};



#endif /* MC_singleHoleQState_hpp */



// ROUTINE TO TEST permute()
/*
 InitializeRandomNumberGenerator();
 
 singleHoleQState qS, qSNew;
 pathData pD;
 PermutationCycles cycs;
 Permutation perm;
 int i1,i2,j;
 bool isMember;
 vector< vector<int> > permMx;
 
 for(j=0; j<100000; j++) {
 qS.random();
 pD.random(RandomInteger(0,100));
 cycs = pD.getPathCycles();
 
 qSNew = qS;
 qSNew.permute(cycs);
 perm = ToPermutation(cycs);
 permMx = perm.getPermutation();
 
 for(i1=0; i1<perm.size(); i1++) {
 if(qS[ permMx[0][i1] ] != qSNew[ permMx[1][i1] ]) {
 cout << "Err: permutations wrong!" << endl;
 cout << qS << endl;
 cout << qSNew << endl;
 cout << perm << endl;
 }
 }
 
 for(i1=0; i1<qS.size(); i1++) {
 isMember = false;
 for(i2=0; i2<permMx[0].size(); i2++) {
 if(permMx[0][i2] == i1) {
 isMember = true;
 break;
 }
 }
 if(!isMember) {
 if(qSNew[i1] != qS[i1]) {
 cout << "Err: permutation moved elements it should not have!" << endl;
 cout << qS << endl;
 cout << qSNew << endl;
 cout << perm << endl;
 }
 }
 }
 }
 */

