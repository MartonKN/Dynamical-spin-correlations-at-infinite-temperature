//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_pathData.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/19/16.
//
//

#include "MC_pathData.hpp"


// -----------------------------
// DATA STRUCTURE TO STORE PATHS
// -----------------------------
pathData::pathData(void) {
    this->coefficient.real(0);
    this->coefficient.imag(0);
    this->endpoint = coordToIndex(0,0);
} // Tested and works


pathData::~pathData(void) {
    
} // Tested and works


// generate random path of length 'len'
void pathData::random(int len) {
    if(len<0) {
        cerr << "Error in pathData::random. The requested path length was negative. We set the path length to be zero.";
        len = 0;
    }
    
    if(len==0) {
        // If the length is zero, we just generate the trivial cycle
        PermutationCycles trivialCyc;
        this->pathCycles = trivialCyc;
        this->endpoint = coordToIndex(0,0);
        this->coefficient.real(1.);
        this->coefficient.imag(0.);
    } else {
        // Randomly generate a set of steps ('steps'), and calculate the corresponding path ('siteIndices')
        vector<char> choices = {'u','d','l','r'};
        vector<char> steps = RandomChoice(choices, len);
        vector<int>  siteIndices(len+1);
        int i;
        int x,y;
        
        x=0; y=0;
        siteIndices[0] = coordToIndex(x,y);
        for(i=0; i<len; i++) {
            switch(steps[i]) {
                case 'u':
                    y++;
                    break;
                case 'd':
                    y--;
                    break;
                case 'l':
                    x--;
                    break;
                case 'r':
                    x++;
                    break;
            }
            siteIndices[i+1] = coordToIndex(x,y);
        }
        
        // Calculate the hole's effect on the lattice spins
        Permutation p;
        p = pathToPermutation(siteIndices);
        this->pathCycles = ToCycles(p);
        
        // Make sure we are not using more memory than we should
        this->pathCycles.shrink_to_fit();
        
        // Set endpoint
        this->endpoint = siteIndices.back();
        
        // Set coefficient = (-ii)^len;
        switch(len % 4) {
            case 0:
                this->coefficient.real( 1.);
                this->coefficient.imag( 0.);
                break;
            case 1:
                this->coefficient.real( 0.);
                this->coefficient.imag(-1.);
                break;
            case 2:
                this->coefficient.real(-1.);
                this->coefficient.imag( 0.);
                break;
            case 3:
                this->coefficient.real( 0.);
                this->coefficient.imag( 1.);
                break;
        }
    }
} // Tested and works


// assignment operator
pathData & pathData::operator=  (const pathData & pD) {
    this->pathCycles = pD.pathCycles;
    this->endpoint = pD.endpoint;
    this->coefficient = pD.coefficient;
    return *this;
} // Tested and works







// -------------------------------
// ROUTINES USED BY THE SIMULATION
// -------------------------------
// This routine determines the permutation representation of the effect of the hole moving over sites over 'siteIndices'
Permutation pathToPermutation(vector<int> siteIndices) {
    Permutation p;
    
    // If the path only has a single element, then we return a trivial permutation
    if(siteIndices.size()<2) {
        return p;
    }
    
    vector< vector<int> > workArray(2, vector<int>(siteIndices.size()));
    // The first row of workArray stores all the indices of the sites visited, without repetition
    // The second row shows which site ended up at that position
    // Thus, workArray is the inverse of the hole's effect on the spins
    
    int iSite;          // index of the position of the hole in siteIndices
    int nVisitedSites;  // number of sites visited so far
    int latestSite;     // the latest site visited by the hole
    int iLatestSite;    // the index in workArray[0] of the latest site visited
    int iPreviousSite;  // the index in workArray[0] of the previous site visited
    bool isNew;         // whether the latestSite has been visited already or not
    int i1;             // auxiliary index
    int tmp;            // temporary variable
    
    // The first step is trivial - the hole is in siteIndices[0]
    workArray[0][0] = siteIndices[0];
    workArray[1][0] = siteIndices[0];
    iLatestSite = 0;
    iPreviousSite = 0;
    
    // Loop over all sites in siteIndices
    for(iSite=1, nVisitedSites=1; iSite<siteIndices.size(); iSite++) {
        latestSite = siteIndices[iSite]; // the latest site visited by the hole
        
        // check if the latest site has been visited already
        // the visited sites are stored in workArray[0][0...(nVisitedSites-1)]
        isNew = true;
        for(i1=0; i1<nVisitedSites; i1++) {
            if(latestSite == workArray[0][i1]) {
                isNew = false;
                iLatestSite = i1;
                break;
            }
        }
        
        
        if(isNew) {
            // if the site was new, add it to the first row of workArray
            iLatestSite = nVisitedSites;
            workArray[0][iLatestSite] = latestSite;
            workArray[1][iLatestSite] = workArray[1][iPreviousSite];
            workArray[1][iPreviousSite] = latestSite;
            nVisitedSites++;
        } else {
            // if it is not new, then we simply swap the elements at iLatestSite and iPreviousSite
            tmp = workArray[1][iPreviousSite];
            workArray[1][iPreviousSite] = workArray[1][iLatestSite];
            workArray[1][iLatestSite] = tmp;
        }
        
        // Update the position of the hole
        iPreviousSite = iLatestSite;
    }
    
    workArray[0].resize(nVisitedSites);
    workArray[1].resize(nVisitedSites);
    
    p = workArray;
    p.invert();
    return p;
} // Tested and works

