//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_Permutation.cpp
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

#include "MC_Permutation.hpp"

// PROTECTED ROUTINES
// Checks if perm is indeed a 2xn matrix
bool Permutation::sizeOK(void) const {
    return (this->perm[0].size() == this->perm[1].size());
} // Tested and works


// Checks if the first and second row have the same set of elements, and each element appears only once
bool Permutation::elementsOK(void) const {
    int i1,i2;
    int nUp, nDown;
    
    // Check first if the first and second row have the same number of elements
    if(!this->sizeOK())
        return false;
    
    // Go over these elements and count how many there are of each of them
    for(i1=0; i1<perm[0].size(); i1++) {
        for(i2=0, nUp=0, nDown=0; i2<perm[1].size(); i2++) {
            if(perm[0][i1] == perm[1][i2])
                nDown++;
            if(perm[0][i1] == perm[0][i2])
                nUp++;
        }
        // all elements have to be there exactly once in each row
        // if this is not true, return false
        if(nUp!=1 || nDown!=1)
            return false;
    }
    
    return true;
} // Tested and works


// Checks if there are single permutation cycles
bool Permutation::noTrivialCycles(void) const {
    for(int i=0; i<perm[0].size(); i++) {
        if(perm[0][i] == perm[1][i])
            return false;
    }
    return true;
} // Tested and works


// Removes single permutation cycles
void Permutation::removeTrivialCycles(void) {
    int i;
    vector<int> singleCycleIndices(0);
    
    // find the indices where we have trivial permutation cycles
    for(i=0; i<perm[0].size(); i++) {
        if(perm[0][i] != perm[1][i])
            singleCycleIndices.push_back(i);
    }
    
    // now, copy those elements to the beginning of perm, and erase everything else
    for(i=0; i<singleCycleIndices.size(); i++) {
        perm[0][i] = perm[0][singleCycleIndices[i]];
        perm[1][i] = perm[1][singleCycleIndices[i]];
    }
    perm[0].resize(singleCycleIndices.size());
    perm[1].resize(singleCycleIndices.size());
} // Tested and works


// CONSTRUCTORS AND DESTRUCTOR
// trivial constructor
Permutation::Permutation(void){
    this->perm.resize(2);
} // Tested and works


// copies permutation from the 2xn matrix to *this
Permutation::Permutation(vector< vector<int> > mx) {
    this->perm.resize(2);
    // Check if mx is of appropriate size, otherwise do nothing
    if(mx.size() == 2) {
        if(mx[0].size() == mx[1].size()) {
            this->perm[0].resize(mx[0].size());
            this->perm[1].resize(mx[1].size());
            for(int i=0; i<mx[0].size(); i++) {
                this->perm[0][i] = mx[0][i];
                this->perm[1][i] = mx[1][i];
            }
            
            this->removeTrivialCycles();
            
            // TESTS HAVE BEEN REMOVED FOR SPEEDUP
            /*
            // Check if everything went well
             if(!this->check()) {
                 cerr << "Error in Permutation::Permutaiton. The permutation could not be copied." << endl;
                 this->perm[0].resize(0);
                 this->perm[1].resize(0);
             }
            */
        }
    }
} // Tested and works

// copy constructor
Permutation::Permutation(const Permutation & p) {
    this->perm.resize(2);
    this->perm[0] = p.perm[0];
    this->perm[1] = p.perm[1];
} // Tested and works

// destructor
Permutation::~Permutation(void) {
    
} // Tested and works


// CHECK IF THE PERMUTATION AT HAND SATISFIES THE RULES
bool Permutation::check(void) const {
    if(!(this->sizeOK())) {
        cerr << "Error in Permutations::check(). The size of the first and second row are not equal." << endl;
        return false;
    }
    if(!(this->noTrivialCycles())) {
        cerr << "Error in Permutations::check(). There are trivial cycles in the permutation." << endl;
        return false;
    }
    if(!(this->elementsOK())) {
        cerr << "Error in Permutations::check(). The permutation is ill defined." << endl;
        return false;
    }
    return true;
} // Tested and works


// GET BASIC DATA
// get the first row of permutations, containing all elements
vector<int> Permutation::getElements(void) const {
    return this->perm[0];
} // Tested and works

// get the full permutation representation perm, as a 2xn matrix
vector< vector<int> > Permutation::getPermutation(void) const {
    return this->perm;
} // Tested and works

// number of elements in the permutation
int Permutation::size(void) const {
    return this->perm[0].size();
} // Tested and works


// PRINTING
std::ostream& operator<< (std::ostream &out, const Permutation &p) {
    int i;
    
    // Print first line
    cout << "(";
    for(i=0; i<p.perm[0].size(); i++) {
        cout << p.perm[0][i];
        if(i < (p.perm[0].size()-1))
            cout << ", ";
    }
    cout << ")" << endl;
    
    // Print second line
    cout << "(";
    for(i=0; i<p.perm[1].size(); i++) {
        cout << p.perm[1][i];
        if(i < (p.perm[1].size()-1))
            cout << ", ";
    }
    cout << ")" << endl;
    
    return out;
} // Tested and works


// ASSIGNMENT OPERATOR
Permutation& Permutation::operator= (const Permutation & p) {
    int len = p.size();
    this->perm[0].resize(len);
    this->perm[1].resize(len);
    
    for(int i1=0; i1<len; i1++) {
        this->perm[0][i1] = p.perm[0][i1];
        this->perm[1][i1] = p.perm[1][i1];
    }
    
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    /*
     if(!this->check()) {
         cerr << "Error in Permutation::operator=. Permutation does not fulfill basic requirements." << endl;
         this->perm[0].resize(0);
         this->perm[1].resize(0);
     }
    */
    
    return *this;
} // Tested and works


// CONVERT TO PERMUTATION CYCLES
// get the decomposition of the permutation into cycles
PermutationCycles ToCycles(Permutation p) {
    int i1, i2;                     // indices
    int tmp;                        // temporary variable used when swapping two elements
    int len = p.perm[0].size();     // number of elements in the permuted set
    int firstElementOfCycle;        // first element of each cycle ...
    int nextElement;                // when we go through the cycle, we store its next element
    
    Cycle cyc;                      // latest cycle
    PermutationCycles permcycles;   // all cycles generated by the permutation p
    vector<int> myCyc;              // vector to store the elements of the latest cycle
    myCyc.reserve(len);

    
    for(i1=0; i1<len; ) {
        firstElementOfCycle = p.perm[0][i1];
        nextElement = p.perm[1][i1];
        myCyc.resize(1);
        myCyc[0]=firstElementOfCycle;
        
        if(nextElement != firstElementOfCycle) {
            while (nextElement!=firstElementOfCycle) {
                // we look for nextElement at indices after i1
                i2 = i1+1;
                while(p.perm[0][i2] != nextElement) {
                    i2++;
                }
                
                // if we found it, we update nextElement ...
                nextElement = p.perm[1][i2];
                
                // ... and we swap the (i1+1)th column of p with the i2th one
                tmp = p.perm[0][i1+1]; p.perm[0][i1+1] = p.perm[0][i2]; p.perm[0][i2] = tmp;
                tmp = p.perm[1][i1+1]; p.perm[1][i1+1] = p.perm[1][i2]; p.perm[1][i2] = tmp;
                
                // store the next element of the cycle in myCyc
                myCyc.push_back(p.perm[0][i1+1]);
                
                // increase i1
                i1++;
            }
            // generate a Cycle out of the permutation cycle we found
            cyc = myCyc;
            permcycles.push_back(cyc);
        }
        
        // we increase i1 when we look for another cycle
        i1++;
    }
    return permcycles;
} // Tested and works



// Converts the permutation cycles permcyc into a permutation. It is defined in MC_permutation.cpp
Permutation ToPermutation(const PermutationCycles &permcyc) {
    Permutation p;
    int i1,i2,i;
    vector<int> cyclens = permcyc.getCycleLengths();
    int len;
    
    // Get number of elements in the cycle
    for(i=0, len=0; i<cyclens.size(); i++)
        len += cyclens[i];
    
    // Resize the permutation accordingle
    p.perm[0].resize(len);
    p.perm[1].resize(len);
    
    // Fill permutation with elements of the cycle
    i=0;
    for(i1=0; i1<cyclens.size(); i1++) {
        for(i2=0; i2<cyclens[i1]; i2++) {
            p.perm[0][i+i2] = (permcyc.allCycles)[i1][i2];
            p.perm[1][i+i2] = (permcyc.allCycles)[i1][(i2+1)%cyclens[i1]];
        }
        i+=cyclens[i1];
    }
    
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    /*
     // Check if p fulfills the requirements of a permutation
     if(!p.check())
         cerr << "Error in Permutation::ToPermutation()." << endl;
    */
    
    return p;
} // Tested and works


// MULTIPLICATION
// (p1*p2) corresponds to the permutation x -> p1(p2(x)). I.e. p2 is evaluated first, then p1.
Permutation operator* (const Permutation & p1, const Permutation & p2) {
    int nElements; // length of the new permutation (number of elements)
    int element;   // auxiliary variable
    int i1, i2;    // iterators
    
    // Declare the result of the multiplication, and reserve enough memory for it
    Permutation result;
    result.perm[0].resize(p1.size() + p2.size());
    result.perm[1].resize(p1.size() + p2.size());
    
    // array to store if an element of the first row of p1.perm has been visited already
    bool visited[p1.size()];
    for(i1=0; i1<p1.size(); i1++)
        visited[i1] = false;
    
    // loop over p2 first, and concatenate the permutation with that in p1
    for(i2=0; i2<p2.size(); i2++) {
        result.perm[0][i2] = p2.perm[0][i2];
        element = p2.perm[1][i2];
        if(p1.findelement(element, i1)) {
            result.perm[1][i2] = p1.perm[1][i1];
            visited[i1] = true;
        } else {
            result.perm[1][i2] = element;
        }
    }
    
    // then, we go over those elements in p1, that we have not visited yet
    nElements = p2.size();
    for(i1=0; i1<p1.size(); i1++) {
        if(!visited[i1]) {
            result.perm[0][nElements] = p1.perm[0][i1];
            result.perm[1][nElements] = p1.perm[1][i1];
            nElements++;
        }
    }
    
    // resize array in result according to the number of elements, and free extra memory
    result.perm[0].resize(nElements);
    result.perm[1].resize(nElements);
    result.removeTrivialCycles();
    result.perm[0].shrink_to_fit();
    result.perm[1].shrink_to_fit();
    
    
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    /*
     if(!result.check()) {
            cerr << "Error in Permutation::operator*. The resulting array does not correspond to a permutation." << endl;
     }
    */
    
    return result;
} // Tested and works


// INVERSION
Permutation& Permutation::invert(void) {
    // We just need to flip the two rows of perm to invert the permutation
    int tmp;
    for(int i=0; i<this->perm[0].size(); i++) {
        tmp = perm[0][i];
        this->perm[0][i] = this->perm[1][i];
        this->perm[1][i] = tmp;
    }
    
    // TESTS HAVE BEEN REMOVED FOR SPEEDUP
    /*
     // Checking
     if(!this->check()) {
         cerr << "Error in Permutation::invert(). The resulting permutation does not fulfill basic requirements." << endl;
     }
    */
    
    return *this;
} // Tested and works


// FIND ELEMENT
// Looks for 'element' in the first row of 'perm'
// Returns true if it found it, and puts its index into i
// Returns false if it didn't find it
bool Permutation::findelement(const int &element, int & i) const {
    for(int i1=0; i1<perm[0].size(); i1++) {
        if(element == perm[0][i1]) {
            i = i1;
            return true;
        }
    }
    return false;
} // Tested and works


