//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_pathData.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/19/16.
//
//

#ifndef MC_pathData_hpp
#define MC_pathData_hpp

#include <stdio.h>
#include <cmath>
#include <complex>
#include <iostream>
#include "MC_clusterRunQ.hpp"
#include "MC_random.hpp"
#include "MC_coord.hpp"
#include "MC_Permutation.hpp"
#include "MC_PermutationCycles.hpp"

using namespace std;


// -----------------------------
// DATA STRUCTURE TO STORE PATHS
// -----------------------------
// Data structure used to store the important parameters of the path in a std::set<pathData> container (self-balancing binary tree)
// In order to perform our computations, we need to store the hole's effect on the spins of the lattice, the endpoint of the path
// and the coefficient of the path, (-ii)^length(path). If we draw the path multiple times, then we add up these coefficients
//
// It is essential to have an ordering between pathData. Our choice here is based solely on the ordering between pathCycles elements.
// As these are permutation cycles, we can use their natural ordering, defined in MC_PermutationCycles.hpp.
// This ordering allows us to store the paths in a red-black self-balancing binary tree, realized by the std::set class.
class pathData {
protected:
    PermutationCycles pathCycles;
    int endpoint;
    mutable complex<long double> coefficient;
    // 'coefficient' will need to be updated when the pathData is stored in a 'std::set' container
    // Therefore, 'coefficient' has to be defined mutable, since 'std::set' only allows const routines to act on its members
    
public:
    pathData(void);        // empty constructor
    ~pathData(void);       // destructor
    void random(int len);  // generate random path of length 'len'
    
    // get basic data
    inline int getEndpoint(void) const {return endpoint;}                                   // Tested and works
    inline complex<long double> getCoefficient(void) const {return coefficient;}            // Tested and works
    inline const PermutationCycles & getPathCycles(void) const {return pathCycles;}         // Tested and works
    
    // update the value of coefficient ('coefficient' is mutable, so these routines can be declared const)
    inline void addToCoefficient(const complex<long double> & c) const {coefficient += c;}  // Tested and works
    inline void setCoefficient  (const complex<long double> & c) const {coefficient = c;}   // Tested and works
    
    // estimate memory used in bytes
    size_t memoryBytes(void) const {return (pathCycles.memoryBytes() + sizeof(int) + 2*sizeof(long double));}
    
    // assignment operator
    pathData & operator=  (const pathData & pD);
    
    // ordering, solely based on the natural ordering between the permutation cycles that the paths generate
    inline friend bool operator<  (const pathData& p1, const pathData& p2) {return (p1.pathCycles <  p2.pathCycles);} // Tested and works
    inline friend bool operator<= (const pathData& p1, const pathData& p2) {return (p1.pathCycles <= p2.pathCycles);} // Tested and works
    inline friend bool operator>  (const pathData& p1, const pathData& p2) {return (p1.pathCycles >  p2.pathCycles);} // Tested and works
    inline friend bool operator>= (const pathData& p1, const pathData& p2) {return (p1.pathCycles >= p2.pathCycles);} // Tested and works
    inline friend bool operator== (const pathData& p1, const pathData& p2) {return (p1.pathCycles == p2.pathCycles);} // Tested and works
    inline friend bool operator!= (const pathData& p1, const pathData& p2) {return (p1.pathCycles != p2.pathCycles);} // Tested and works
    inline friend int smallerPathData(const pathData& p1, const pathData& p2) {return smallerPermutationCycles(p1.pathCycles, p2.pathCycles);}
    // smallerPathData returns 1 if p1<p2, returns -1 if p1>p2, and returns 0 if they are equal. Its use is advisable when multiple comparisons of p1 and p2 are necessary, e.g. if(p1<p2){...} else if(p1==p2) {...}. By writing this as switch(smallerPathData(p1,p2)) {...}, we need to compare the paths only once instead of twice.
    
    // printing
    friend std::ostream & operator<< (ostream & out, const pathData & p) {
        cout << "endpoint:    " << p.endpoint << " = (" << indexToX(p.endpoint) << ", " << indexToY(p.endpoint) << ")" << endl;
        cout << "coefficient: " << p.coefficient << endl;
        cout << "pathCycles:  " << p.pathCycles << endl;
        return out;
    } // Tested and works
}; // All methods have been tested, and they work appropriately




// -------------------------------
// ROUTINES USED BY THE SIMULATION
// -------------------------------
// This routine determines the permutation representation of the effect of the hole moving over sites over 'siteIndices'
Permutation pathToPermutation(vector<int> siteIndices);



#endif /* MC_pathData_hpp */
