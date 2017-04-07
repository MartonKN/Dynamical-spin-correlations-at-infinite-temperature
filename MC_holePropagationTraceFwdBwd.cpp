//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_holePropagationTraceFwdBwd.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 1/10/17.
//
//

#include "MC_holePropagationTraceFwdBwd.hpp"


// creates the filenames based on system parameters
// saveFilenames[0] = file containing all the filenames
// saveFilenames[1] = system parameters
// saveFilenames[2] = return probability and rms distance data of each run
// saveFilenames[3] = transition probabilities
// saveFilenames[4] = lab frame spin correlation function
// saveFilenames[5] = hole frame spin correlation function
// saveFilenames[6] = transProbsBins
// saveFilenames[7] = LabFrameSpinCorrelationBins
// saveFilenames[8] = HoleFrameSpinCorrelationBins
void holePropagationTraceFwdBwd::generateSaveFilenames(void) {
    stringstream fnameBaseStream;
    string fnameBase;
    
    // Create the base of all the file names
    fnameBaseStream << "InfiniteT_TraceFwdBwd_ver1_10_t=" << setprecision(2) << fixed << this->tHole;
    fnameBaseStream << "_lattice=" << 2*LATTICE_SIZE << "x" << 2*LATTICE_SIZE;
    fnameBaseStream << "_nPaths=" << this->nPaths;
    if(this->calcSpinCorrelationsQ) {
        fnameBaseStream << "_wSpinCorr";
    } else {
        fnameBaseStream << "_noSpinCorr";
    }
    fnameBase = fnameBaseStream.str();
    
    // create all the filenames
    this->saveFilenames.resize(9);
    this->saveFilenames[0] = fnameBase + "_files.txt";
    this->saveFilenames[1] = fnameBase + "_params.txt";
    this->saveFilenames[2] = fnameBase + "_returnProb_and_RMS.txt";
    this->saveFilenames[3] = fnameBase + "_transProbs.txt";
    this->saveFilenames[4] = fnameBase + "_spinCorrLabFrame.txt";
    this->saveFilenames[5] = fnameBase + "_spinCorrHoleFrame.txt";
    this->saveFilenames[6] = fnameBase + "_transProbsBins.txt";
    this->saveFilenames[7] = fnameBase + "_spinCorrLabFrameBins.txt";
    this->saveFilenames[8] = fnameBase + "_spinCorrHoleFrameBins.txt";
} // Tested and works





// generates new pathContainer and updates the bins using holePropagationTraceFwdBwd::update()
void holePropagationTraceFwdBwd::run(void) {
    // Generate path container
    pathContainer contFwd;
    pathContainer contBwd;
    
    // Fill path container with randomly generated paths
    contFwd.fill(this->tHole, this->nPaths);
    contBwd.fill(this->tHole, this->nPaths);
    
    cout << "Running simulation " << (this->nRuns + 1) << endl;
    cout << "Number of pathData bins: (" << contFwd.getNumberOfBins() << ", " << contBwd.getNumberOfBins() << ")" << endl;
    cout << "Estimated memory used by paths: " << setprecision(2) << fixed << ((contFwd.getMemoryBytes() + contBwd.getMemoryBytes())/(1024.*1024.)) << " MB" << endl;
    cout << endl << endl;
    
    // Calculate the contribution of each pair of path to the transition probability
    this->update(contBwd, contFwd);
    
    // Update number of runs
    (this->nRuns)++;
} // Checked the code, looks good


// updates the bins with contributions from all pairs of paths in 'contFwd' and 'contBwd'
void holePropagationTraceFwdBwd::update(const pathContainer &contBwd, const pathContainer &contFwd) {
    int iSite;                          // site index
    set<pathData>::iterator it1, it2;   // iterator of pathData elements within cont[iSite]
    pathData pDfwd, pDbwd;              // forward and backward paths
    
    // Precalculate the values 1/NUMBER_OF_SPINS^n, so that this does not rob time during the simulation
    vector<double >invNSpinsPowers(61); // 0.5^60 ~ 1.E-18, we do not take into account contributions that are smaller
    invNSpinsPowers.shrink_to_fit();
    long double tmp = 1.0L;
    for(int i=0; i<invNSpinsPowers.size(); i++) {
        invNSpinsPowers[i] = tmp;
        tmp /= ((long double) NUMBER_OF_SPINS);
    }
    
    // Loop over all sites as well as all pairs of paths within the container at that site, and update *this accordingly
    for(iSite=0; iSite<(this->nSites); iSite++) {
        for(it1=contFwd[iSite].begin(); it1!=contFwd[iSite].end(); ++it1) {
            pDfwd = (*it1);
            for(it2=contBwd[iSite].begin(); it2!=contBwd[iSite].end(); ++it2) {
                pDbwd = (*it2);
                this->update(pDbwd,pDfwd,invNSpinsPowers);
            }
        }
    }
}


// Updates the bins with contributions from the pair of paths 'pDfwd' and 'pDbwd',and returns the contribution of this pair to the transition probabilities.
// invNSpinsPowers is an auxiliary array: invNSpinsPowers[i] = 1/NUMBER_OF_SPINS^i, that needs to be determined in update(pathContainer &)
// Note that we should add the transition probability contribution of the pair of paths to all same site spin correlation bins, except for the hole since this would lead to unnecessary computational overhead, we actually remove this contribution from the hole site. We will compensate for this at the end, in getLabFrameSpinCorrelationFn() and getHoleFrameSpinCorrelationFn().
double holePropagationTraceFwdBwd::update(const pathData &pDbwd, const pathData &pDfwd, const vector<double> &invNSpinsPowers) {
    
    // Check if the two paths have the same endpoint (they should)
    if(pDfwd.getEndpoint() != pDbwd.getEndpoint()) {
        cerr << "Error in holePropagationTraceFwdBwd::update. The two paths have different endpoints. " << endl;
        return 0.;
    }
    
    // Check if one of the coefficients is zero. Then, we don't have to do anything
    if((pDbwd.getCoefficient() == complex<long double>(0.,0.)) || (pDfwd.getCoefficient() == complex<long double>(0.,0.))) {
        return 0.;
    }
    
    // Define variables
    long double contribution;             // contribution of this pair of paths to the transition probability from the origin to holePosition
    int exponent;                         // exponent = sum{(length(permutation cycle) - 1)} over all cycles generated by the combined path
    int holePosition=pDfwd.getEndpoint(); // position of the hole
    PermutationCycles combinedPath;       // combination of the permutation effect of the forward and backward path
    PermutationCycles fwdPath;            // permutations generated by the forward path
    vector<int> cycleLengths;             // length of cycles in 'combinedPath'
    vector<int> cycleElements;            // auxiliary variable to store the elements of each cycle in 'combinedPath'
    vector<int> cycleElementsComoving;    // same as cycleElements, just in the co-moving frame of the hole
    Cycle cyc;                            // auxiliary Cycle variable
    int x, y, xHole, yHole;               // auxiliary coordinate variables to determine 'cycleElementsComoving' from 'cycleElements'
    complex<long double> tmp;             // auxiliary variable used for multiplication of complex numbers
    int i,i1,i2;                          // indices
    
    
    // Get the permutation cycles representation of the combined path
    fwdPath = pDfwd.getPathCycles();
    combinedPath = pDbwd.getPathCycles();
    combinedPath.invert();
    combinedPath = combinedPath * fwdPath;
    
    // Get the length of cycles in the combined path
    cycleLengths = combinedPath.getCycleLengths();
    for(exponent=0, i=0; i<cycleLengths.size(); i++)
        exponent += (cycleLengths[i]-1);
    
    // Calculate the contribution of this pair of paths to the transition probabilities
    tmp = conj(pDbwd.getCoefficient()) * pDfwd.getCoefficient();
    contribution = tmp.real();
    if(tmp.imag() != 0.) {
        cerr << "Error in holePropagationTraceFwdBwd::update(). The contribution of the coefficient is not real. " << endl;
    }
    if(exponent <= invNSpinsPowers.size()) {
        contribution *= invNSpinsPowers[exponent];
    } else {
        return 0.;
    }
    
    // Update transition probabilities
    transProbsBins[holePosition] += contribution;
    
    // Update spin correlations (neglecting same-site correlations, which are trivially 1)
    // Let us denote the permutation effect of the forward path as 'fwd' and that of the backward path as 'bwd'
    // We get non-zero spin correlations between the sites that are in the same permutation cycle of inverse(bwd)*fwd == combinedPath.
    // Let us denote the elements of one of these cycles as {a, b, c, d, e, f}.
    // Then, the sites among which we will find correlations are {fwd(a), fwd(b), fwd(c), fwd(d), fwd(e), fwd(f)}
    if(this->calcSpinCorrelationsQ) {
        LabFrameSpinCorrelationBins [holePosition][holePosition]           -= contribution;
        HoleFrameSpinCorrelationBins[coordToIndex(0,0)][coordToIndex(0,0)] -= contribution;
        for(i=0; i<combinedPath.size(); i++) {
            cyc = combinedPath[i];          // choose an element of the permutation cycle
            cycleElements = cyc.elements(); // get elements of the cycle
            fwdPath.permute(cycleElements); // permute the sites according to the forward path
            
            // 'combinedPath' contributes to the spin correlations between all unequal elements in 'cycleElements'
            // As for spin correlations on the same site 'i', those are equal to 1-Probability(the hole is at site 'i')
            // First, update the lab frame spin correlations
            for(i1=0; i1<cycleElements.size(); i1++) {
                for(i2=i1+1; i2<cycleElements.size(); i2++) {
                    LabFrameSpinCorrelationBins[cycleElements[i1]][cycleElements[i2]] += contribution;
                    LabFrameSpinCorrelationBins[cycleElements[i2]][cycleElements[i1]] += contribution;
                }
            }
            // second, update the spin correlations in the hole's comoving frame
            xHole = indexToX(holePosition);
            yHole = indexToY(holePosition);
            cycleElementsComoving.resize(cycleElements.size());
            for(i1=0; i1<cycleElements.size(); i1++) {
                x = indexToX(cycleElements[i1]);
                y = indexToY(cycleElements[i1]);
                cycleElementsComoving[i1] = coordToIndex(x-xHole, y-yHole);
            }
            for(i1=0; i1<cycleElements.size(); i1++) {
                for(i2=i1+1; i2<cycleElements.size(); i2++) {
                    HoleFrameSpinCorrelationBins[cycleElementsComoving[i1]][cycleElementsComoving[i2]] += contribution;
                    HoleFrameSpinCorrelationBins[cycleElementsComoving[i2]][cycleElementsComoving[i1]] += contribution;
                }
            }
        }
    }
    return contribution;
}

