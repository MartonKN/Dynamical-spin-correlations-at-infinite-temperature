//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_Cycle.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/2/16.
//
//

#ifndef MC_Cycle_hpp
#define MC_Cycle_hpp

#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include "MC_clusterRunQ.hpp"

using namespace std;


// Class implementing a single permutation cycle
class Cycle {
protected:
    vector<int> Elements;
    // The elements of the cycle.
    // They need to be stored in a way that element[0] is the lowest one of all of them
    // The cycle cannot contain multiple copies of any element
    
public:
    // CONSTRUCTORS AND DESTRUCTOR
    Cycle(void);                    // Creates a cycle without any elements
    Cycle(const vector<int> & el);  // Copies the cycle in 'el', rotating the vector such that the cycle starts with the smallest element of 'el'
    ~Cycle(void);                   // Destructor
    
    // GET BASIC DATA
    int size(void) const;
    vector<int> elements(void) const;
    
    // CHECK IF THE ELEMENTS FULFILL THE REQUIREMENTS WE NEED
    bool check(void) const; // Check if 'Elements' indeed starts with its smallest element, and if it has any duplicates
    
    // PRINTING
    friend std::ostream& operator<< (std::ostream &out, const Cycle & cyc);
    
    // CIRCUMVENT MEMORY WASTE
    void shrink_to_fit(void);           // Requests the compiler to use only as much memory space as the size of the cycle
                                        // This way we can avoid excess memory usage
    size_t memoryBytes(void) const;     // Estimates memory used by our permutation cycle
    
    // OVERLOADED OPERATORS
    Cycle& operator= (const Cycle &cyc);        // assignment operator
    Cycle& operator= (const vector<int> &el);   // assignment operator
    inline const int & operator[] (const int &i) const // elements of the cycle (read only)
    {
        // We deliberately do not check if i is within the lenght of the Cycle.
        // If it is out of bound, we let the error thrown by std::vector to propagate through the system.
        return this->Elements[i];
    }
    
    // ORDERING
    // Routine introducing the ordering between different cycles
    // Returns 1 if cyc1<cyc2, returns -1 if cyc1>cyc2, and returns 0 if they are equal
    // Condition of being smaller:
    // 1. whoever has less elements
    // 2. if they have the same number of elements, we go over their elements, and compare
    // them one-by-one until one of them will be smaller.
    friend int smallerCycle(const Cycle &cyc1, const Cycle &cyc2);
    
    // All other ordering relations are based on the results of Cycle::smaller
    // All of these have been tested and work properly
    inline friend bool operator== (const Cycle &cyc1, const Cycle &cyc2) {
        return (smallerCycle(cyc1,cyc2) == 0);
    }
    inline friend bool operator!= (const Cycle &cyc1, const Cycle &cyc2) {
        return !(smallerCycle(cyc1,cyc2) == 0);
    }
    inline friend bool operator<  (const Cycle &cyc1, const Cycle &cyc2) {
        return (smallerCycle(cyc1,cyc2) == 1);
    }
    inline friend bool operator<=  (const Cycle &cyc1, const Cycle &cyc2) {
        return (smallerCycle(cyc1,cyc2) != -1);
    }
    inline friend bool operator>  (const Cycle &cyc1, const Cycle &cyc2) {
        return (smallerCycle(cyc1,cyc2) == -1);
    }
    inline friend bool operator>=  (const Cycle &cyc1, const Cycle &cyc2) {
        return (smallerCycle(cyc1,cyc2) != 1);
    }
    
    // INVERSION
    Cycle & invert(void);
    
    // Routine to check if two cycles have a common element
    // Returns true if there is a common element and false if there is not
    friend bool hasCommonElement(const Cycle &cyc1, const Cycle &cyc2);
};



#endif /* MC_Cycle_hpp */




/*
 // TEST THAT THE RELATIONS BETWEEN THE CYCLES SATISFY THE ORDERING AXIOMS
 struct timeval tp;
 gettimeofday(&tp, NULL);
 long int time_in_microsec = tp.tv_sec * 1000000 + tp.tv_usec;
 srand (mod(time_in_microsec, RAND_MAX));
 
 
 int l1, l2, l3;
 vector<int> d1, d2, d3;
 Cycle cyc1, cyc2, cyc3;
 bool isOK;
 
 
 for(int i=0; i<1000; i++) {
 l1 = RandomInteger(3,6);
 l2 = RandomInteger(3,6);
 l3 = RandomInteger(3,6);
 
 d1 = RandomInteger(-70,70,l1);
 d2 = RandomInteger(-70,70,l2);
 d3 = RandomInteger(-70,70,l3);
 
 cyc1 = d1;
 cyc2 = d2;
 cyc3 = d3;
 cout << endl;
 
 isOK = true;
 isOK = isOK && (cyc1 == cyc1);
 isOK = isOK && !(cyc1 < cyc1);
 isOK = isOK && !(cyc1 > cyc1);
 isOK = isOK && (cyc1 >= cyc1);
 isOK = isOK && (cyc1 <= cyc1);
 isOK = isOK && !(cyc1<cyc2 && cyc2<cyc1);
 if(cyc1 < cyc2 && cyc2 < cyc3) {
 isOK = isOK && (cyc1 < cyc3);
 }
 if(cyc1 <= cyc2 && cyc2 <= cyc3) {
 isOK = isOK && (cyc1 <= cyc3);
 }
 if(cyc1 <= cyc2 && cyc2 <= cyc1) {
 isOK = isOK && (cyc1 == cyc2);
 }
 
 
 if(!isOK) {
 cout << endl << endl;
 cout << "*************************************************************************************************************** " << endl;
 cout << "ERROR" << endl;
 cout << "cyc1: " << cyc1 << endl;
 cout << "cyc2: " << cyc2 << endl;
 cout << "cyc3: " << cyc3 << endl;
 cout << "*************************************************************************************************************** " << endl;
 cout << endl << endl;
 }
 }

 */
