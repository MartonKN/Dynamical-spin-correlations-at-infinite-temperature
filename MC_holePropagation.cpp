//  Copyright 2017, Marton Kanasz-Nagy, All rights reserved.
//  kanasz.nagy.marton@gmail.com
//
//  MC_HolePropagation.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/6/16.
//
//

#include "MC_holePropagation.hpp"



// ---------------------------------------------
// CLASS RUNNING THE SIMULATION AND STORING DATA
// ---------------------------------------------
// Stores the transition probabilities as well as the spin correlation functions, both in the laboratory frame,
// as well as in the hole's co-moving frame

// constructor: fills all arrays with zeros
holePropagation::holePropagation(void) {
    int i1,i2;
    
    InitializeRandomNumberGenerator();
    
    // set the size of the bins
    this->transProbsBins.reserve(this->nSites);
    this->transProbsBins.resize(this->nSites);
    this->transProbsBins.shrink_to_fit();
    this->LabFrameSpinCorrelationBins.reserve(this->nSites);
    this->LabFrameSpinCorrelationBins.resize(this->nSites);
    this->LabFrameSpinCorrelationBins.shrink_to_fit();
    this->HoleFrameSpinCorrelationBins.reserve(this->nSites);
    this->HoleFrameSpinCorrelationBins.resize(this->nSites);
    this->HoleFrameSpinCorrelationBins.shrink_to_fit();
    for(i1=0; i1<nSites; i1++) {
        this->LabFrameSpinCorrelationBins[i1].reserve(this->nSites);
        this->LabFrameSpinCorrelationBins[i1].resize(this->nSites);
        this->LabFrameSpinCorrelationBins[i1].shrink_to_fit();
        this->HoleFrameSpinCorrelationBins[i1].reserve(this->nSites);
        this->HoleFrameSpinCorrelationBins[i1].resize(this->nSites);
        this->HoleFrameSpinCorrelationBins[i1].shrink_to_fit();
    }
    
    // fill the bins with zeros
    for(i1=0; i1<this->nSites; i1++) {
        this->transProbsBins[i1] = 0.;
        for(i2=0; i2<this->nSites; i2++) {
            this->LabFrameSpinCorrelationBins[i1][i2] = 0.;
            this->HoleFrameSpinCorrelationBins[i1][i2] = 0.;
        }
    }
    this->tHole=0.;
    this->nPaths=0;
    this->nRuns=0;
    this->calcSpinCorrelationsQ=true;
    this->generateSaveFilenames();
} // Tested and works

// constructor: fills all arrays with zeros
holePropagation::holePropagation(double t, int nFwdPaths, bool spinCorrQ) {
    int i1,i2;
    
    InitializeRandomNumberGenerator();
    
    // set the size of the bins
    this->transProbsBins.resize(this->nSites);
    this->transProbsBins.reserve(this->nSites);
    this->transProbsBins.shrink_to_fit();
    this->LabFrameSpinCorrelationBins.reserve(this->nSites);
    this->LabFrameSpinCorrelationBins.resize(this->nSites);
    this->LabFrameSpinCorrelationBins.shrink_to_fit();
    this->HoleFrameSpinCorrelationBins.reserve(this->nSites);
    this->HoleFrameSpinCorrelationBins.resize(this->nSites);
    this->HoleFrameSpinCorrelationBins.shrink_to_fit();
    for(i1=0; i1<nSites; i1++) {
        this->LabFrameSpinCorrelationBins[i1].reserve(this->nSites);
        this->LabFrameSpinCorrelationBins[i1].resize(this->nSites);
        this->LabFrameSpinCorrelationBins[i1].shrink_to_fit();
        this->HoleFrameSpinCorrelationBins[i1].reserve(this->nSites);
        this->HoleFrameSpinCorrelationBins[i1].resize(this->nSites);
        this->HoleFrameSpinCorrelationBins[i1].shrink_to_fit();
    }
    
    // Fill the bins with zeros
    for(i1=0; i1<this->nSites; i1++) {
        this->transProbsBins[i1] = 0.;
        for(i2=0; i2<this->nSites; i2++) {
            this->LabFrameSpinCorrelationBins[i1][i2] = 0.;
            this->HoleFrameSpinCorrelationBins[i1][i2] = 0.;
        }
    }
    this->tHole=t;
    this->nPaths=nFwdPaths;
    this->nRuns=0;
    this->calcSpinCorrelationsQ = spinCorrQ;
    this->generateSaveFilenames();
} // Tested and works


// destructor
holePropagation::~holePropagation(void) {

} // Tested and works


// resets time, #runs, #paths, and clears all values from the bins
void holePropagation::reset(double t, int nFwdPaths) {
    int i1,i2;
    
    InitializeRandomNumberGenerator();
    
    // set the size of the bins
    this->transProbsBins.resize(this->nSites);
    this->transProbsBins.shrink_to_fit();
    this->LabFrameSpinCorrelationBins.resize(this->nSites);
    this->LabFrameSpinCorrelationBins.shrink_to_fit();
    this->HoleFrameSpinCorrelationBins.resize(this->nSites);
    this->HoleFrameSpinCorrelationBins.shrink_to_fit();
    for(i1=0; i1<this->nSites; i1++) {
        this->LabFrameSpinCorrelationBins[i1].resize(this->nSites);
        this->LabFrameSpinCorrelationBins[i1].shrink_to_fit();
        this->HoleFrameSpinCorrelationBins[i1].resize(this->nSites);
        this->HoleFrameSpinCorrelationBins[i1].shrink_to_fit();
    }
    
    // Fill the bins with zeros
    for(i1=0; i1<this->nSites; i1++) {
        this->transProbsBins[i1] = 0.;
        for(i2=0; i2<this->nSites; i2++) {
            this->LabFrameSpinCorrelationBins[i1][i2] = 0.;
            this->HoleFrameSpinCorrelationBins[i1][i2] = 0.;
        }
    }
    this->tHole=t;
    this->nPaths=nFwdPaths;
    this->nRuns=0;
    this->generateSaveFilenames();
} // Tested and works


// clears all values from the bins and resets nRuns to zero
void holePropagation::clearBins(void) {
    int i1,i2;
    for(i1=0; i1<this->nSites; i1++) {
        this->transProbsBins[i1] = 0.;
        for(i2=0; i2<this->nSites; i2++) {
            this->LabFrameSpinCorrelationBins[i1][i2] = 0.;
            this->HoleFrameSpinCorrelationBins[i1][i2] = 0.;
        }
    }
    this->nRuns=0;
} // Tested and works


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
void holePropagation::generateSaveFilenames(void) {
    stringstream fnameBaseStream;
    string fnameBase;
    
    // Create the base of all the file names
    fnameBaseStream << "Carlstrom_ver1_10_t=" << setprecision(2) << fixed << this->tHole;
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


// generates new pathContainer and updates the bins using holePropagation::update()
void holePropagation::run(void) {
    // Generate path container
    pathContainer cont;
    
    // Fill path container with randomly generated paths
    cont.fill(this->tHole, this->nPaths);
    
    cout << "Running simulation " << (this->nRuns + 1) << endl;
    cout << "Number of pathData bins: " << cont.getNumberOfBins() << endl;
    cout << "Estimated memory used by paths: " << setprecision(2) << fixed << (cont.getMemoryBytes()/(1024.*1024.)) << " MB" << endl;
    cout << endl << endl;
    
    // Calculate the contribution of each pair of path to the transition probability
    this->update(cont);
    
    // Update number of runs
    (this->nRuns)++;
} // Checked the code, looks good


// updates the bins with contributions from all pairs of paths in 'cont'
void holePropagation::update(const pathContainer &cont) {
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
        for(it1=cont[iSite].begin(); it1!=cont[iSite].end(); ++it1) {
            pDfwd = (*it1);
            for(it2=cont[iSite].begin(); it2!=cont[iSite].end(); ++it2) {
                pDbwd = (*it2);
                this->update(pDbwd,pDfwd,invNSpinsPowers);
            }
        }
    }
} // Checked the code, looks good


// Updates the bins with contributions from the pair of paths 'pDfwd' and 'pDbwd',
// and returns the contribution of this pair to the transition probabilities.
// Note that we should add the transition probability contributionof the pair of paths to all same site spin correlation bins,
// except for the hole since this would lead to unnecessary computational overhead, we actually remove this contribution
// from the hole site. We will compensate for this at the end, in getLabFrameSpinCorrelationFn() and getHoleFrameSpinCorrelationFn().
// Note furthermore, that on the square lattice, the contribution of (pDfwd, pDbwd) shall be the same as (pDbwd, pDfwd).
// As holePropagation::update() is only used when we loop over all pathData elements of a pathContainer, it is a factor of 2 waste
// of time to take all pairs. Instead, we set holePropagation::update() such that it only updates '*this', when pDfwd<=pDbwd,
// but in this case it takes into account the contributions of both (pDfwd, pDbwd) and (pDbwd, pDfwd), which shall be equal.
double holePropagation::update(const pathData &pDbwd, const pathData &pDfwd, const vector<double> &invNSpinsPowers) {
    // First, we check if pDfwd<pDbwd (sm==1), pDfwd==pDbwd (sm==0) or  pDfwd>pDbwd (sm==-1).
    // If pDfwd>pDbwd, we do not do anything
    int sm = smallerPathData(pDfwd, pDbwd);
    if (sm == -1)
        return 0.;
    
    // Check if the two paths have the same endpoint (they should)
    if(pDfwd.getEndpoint() != pDbwd.getEndpoint()) {
        cerr << "Error in holePropagation::update. The two paths have different endpoints. " << endl;
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
        cerr << "Error in holePropagation::update(). The contribution of the coefficient is not real. " << endl;
    }
    if(exponent <= invNSpinsPowers.size()) {
        contribution *= invNSpinsPowers[exponent];
    } else {
        return 0.;
    }
    
    // In case pDfwd<pDbwd, we need to take into account the contribution of the (pDfwd,pDbwd) pair, as well as the reversed path (pDbwd,pDfwd)
    // In this case, the contributions have to be doubled, to account for this
    if (sm == 1)
        contribution *= 2.;
    
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

// calculate spin correlations using data in LabFrameSpinCorrelationBins
void holePropagation::calculateLabFrameSpinCorrelationFn (vector< vector<double> > & labFrameSpinCorrelationFn) const {
    long double total = 0.L;
    long double invTotal;
    int i,j;
    
    // Resize array and make sure it does not waste memory
    labFrameSpinCorrelationFn.reserve(this->nSites);
    labFrameSpinCorrelationFn.resize(this->nSites);
    labFrameSpinCorrelationFn.shrink_to_fit();
    for(i=0; i<nSites; i++) {
        labFrameSpinCorrelationFn[i].reserve(this->nSites);
        labFrameSpinCorrelationFn[i].resize(this->nSites);
        labFrameSpinCorrelationFn[i].shrink_to_fit();
    }
    
    // Calculate spin correlations only if you have to
    if(!this->calcSpinCorrelationsQ) {
        for(i=0; i<this->nSites; i++) {
            for(j=0; j<this->nSites; j++) {
                labFrameSpinCorrelationFn[i][j] = 0.;
            }
        }
        return;
    }
    
    // Calculate the total of the transition probability bins
    // This will provide normalization for the spin correlations as well
    for(i=0; i<this->nSites; i++) {
        total += transProbsBins[i];
    }
    
    // Normalize spin correlation function
    if(total != 0.) {
        invTotal = 1.0L / total;
        for(i=0; i<this->nSites; i++) {
            for(j=0; j<this->nSites; j++) {
                labFrameSpinCorrelationFn[i][j] = invTotal * LabFrameSpinCorrelationBins[i][j];
            }
        }
    } else {
        for(i=0; i<this->nSites; i++) {
            for(j=0; j<this->nSites; j++) {
                labFrameSpinCorrelationFn[i][j] = 0.;
            }
        }
        cerr << "Error in holePropagation::calculateLabFrameSpinCorrelationFn. Spin correlation function ill defined." << endl;
        return;
    }
    
    // Add the diagonal elements of the spin correlations, that we neglected during the simulation
    for(i=0; i<this->nSites; i++) {
        labFrameSpinCorrelationFn[i][i] += 1.;
    }
}


// calculate spin correlations using data in HoleFrameSpinCorrelationBins
void holePropagation::calculateHoleFrameSpinCorrelationFn(vector< vector<double> > & holeFrameSpinCorrelationFn) const {
    long double total = 0.L;
    long double invTotal;
    int i,j;
    
    // Resize array and make sure it does not waste memory
    holeFrameSpinCorrelationFn.resize(this->nSites);
    holeFrameSpinCorrelationFn.reserve(this->nSites);
    holeFrameSpinCorrelationFn.shrink_to_fit();
    for(i=0; i<nSites; i++) {
        holeFrameSpinCorrelationFn[i].resize(this->nSites);
        holeFrameSpinCorrelationFn[i].reserve(this->nSites);
        holeFrameSpinCorrelationFn[i].shrink_to_fit();
    }
    
    // Calculate spin correlations only if you have to
    if(!this->calcSpinCorrelationsQ) {
        for(i=0; i<this->nSites; i++) {
            for(j=0; j<this->nSites; j++) {
                holeFrameSpinCorrelationFn[i][j] = 0.;
            }
        }
        return;
    }

    // Calculate the total of the transition probability bins
    // This will provide normalization for the spin correlations as well
    for(i=0; i<this->nSites; i++) {
        total += transProbsBins[i];
    }
    
    // Normalize spin correlation function
    if(total != 0.) {
        invTotal = 1.0L / total;
        for(i=0; i<this->nSites; i++) {
            for(j=0; j<this->nSites; j++) {
                holeFrameSpinCorrelationFn[i][j] = invTotal * HoleFrameSpinCorrelationBins[i][j];
            }
        }
    } else {
        for(i=0; i<this->nSites; i++) {
            for(j=0; j<this->nSites; j++) {
                holeFrameSpinCorrelationFn[i][j] = 0.;
            }
        }
        cerr << "Error in holePropagation::calculateholeFrameSpinCorrelationFn. Spin correlation function ill defined." << endl;
        return;
    }
    
    // Add the diagonal elements of the spin correlations, that we neglected during the simulation
    for(i=0; i<this->nSites; i++) {
        holeFrameSpinCorrelationFn[i][i] += 1.;
    }
}


// get the number of runs performed so far
int holePropagation::getNumberOfRuns(void) const {
    return this->nRuns;
}


// returns the file names we will be saving our data in
vector<string> holePropagation::getFilenames(void) const {
    return this->saveFilenames;
}; // Tested and works


// returns the transition probabilities
vector<double> holePropagation::getTransitionProbabilities(void) const {
    vector<double> transProbs(this->nSites);
    transProbs.shrink_to_fit();
    
    long double total = 0.L;
    long double invTotal;
    int i;
    
    for(i=0; i<this->nSites; i++) {
        total += transProbsBins[i];
    }
    
    if(total!=0.) {
        invTotal = 1.0L / total;
        for(i=0; i<this->nSites; i++) {
            transProbs[i] = invTotal * transProbsBins[i];
        }
    } else {
        for(i=0; i<this->nSites; i++) {
            transProbs[i] = 0.;
        }
        cerr << "Error in holePropagation::getTransitionProbabilities. Transition probabilities are ill defined." << endl;
    }
    
    // Symmetry averaging
    transProbs = symmetryAverage(transProbs);
    
    return transProbs;
} // Checked the code, looks good


// returns the lab frame spin correlations
vector< vector<double> > holePropagation::getLabFrameSpinCorrelationFn(void) const {
    vector< vector<double> > labFrameSpinCorrelationFn;
    this->calculateLabFrameSpinCorrelationFn(labFrameSpinCorrelationFn);
    return labFrameSpinCorrelationFn;
    // Note: since the compiler probably performs return value optimization, it is fine to return this large array
} // Checked the code, looks good


// returns the hole frame spin correlations
vector< vector<double> > holePropagation::getHoleFrameSpinCorrelationFn(void) const {
    vector< vector<double> > holeFrameSpinCorrelationFn;
    this->calculateHoleFrameSpinCorrelationFn(holeFrameSpinCorrelationFn);
    return holeFrameSpinCorrelationFn;
    // Note: since the compiler probably performs return value optimization, it is fine to return this large array
} // Checked the code, looks good


// probability of the hole to return to the origin
double holePropagation::getReturnProbability(void) const {
    long double total = 0.L;
    int origin = coordToIndex(0,0);
    double returnProb;
    
    for(int i=0; i<this->nSites; i++) {
        total += transProbsBins[i];
    }
    if(total!=0) {
        returnProb = ((double) transProbsBins[origin]) / total;
    } else {
        cerr << "Error in holePropagation::getReturnProbability. Return probability ill defined." << endl;
        returnProb = 0.;
    }
    
    return returnProb;
} // Checked the code, looks good


// root mean squared distance of the hole from the origin
double holePropagation::getRMSdistance(void) const {
    vector<double> tP; // transition probabilities
    int siteInd;       // site index
    int x,y;           // site coordinates
    double RMS;        // root mean squared distance of the hole from the origin
    
    tP = this->getTransitionProbabilities();
    for(RMS=0., siteInd=0; siteInd<this->nSites; siteInd++) {
        x = indexToX(siteInd);
        y = indexToY(siteInd);
        RMS += (x*x + y*y) * tP[siteInd];
    }
    RMS = sqrt(RMS);
    
    return RMS;
} // Checked the code, looks good


// saves all three bins, using the three saving functions below
bool holePropagation::backup(void) const {
    bool savingSuccessful = true;
    
    savingSuccessful = savingSuccessful && this->saveTransitionProbabilityBins();
    savingSuccessful = savingSuccessful && this->save_returnProb_and_RMSdistance();
    if(this->calcSpinCorrelationsQ) {
        savingSuccessful = savingSuccessful && this->saveLabFrameSpinCorrelationBins();
        savingSuccessful = savingSuccessful && this->saveHoleFrameSpinCorrelationBins();
    }
    
    if(!savingSuccessful) {
        cerr << "Error in holePropagation::backup. Could not save all data." << endl;
    }
    
    return savingSuccessful;
} // Checked the code, looks good


// save transition probability bins
// first element of the file is the number of runs performed so far (this->nRuns).
// all later elements are from transProbsBins.
bool holePropagation::saveTransitionProbabilityBins(void) const {
    string fname = this->saveFilenames[6];
    
    ofstream myfile (fname);
    if(!myfile.good()) {
        cerr << "Error in holePropagation::saveTransitionProbabilityBins writing to " << fname << endl << endl;
        return false;
    }
    
    myfile << (this->nRuns) << endl;
    for(int i1=0; i1<this->nSites; i1++) {
        myfile << setprecision(30) << this->transProbsBins[i1] << endl;
    }
    myfile.close();
    
    return true;
}   // Tested and works


// save lab frame spin correlation bins
bool holePropagation::saveLabFrameSpinCorrelationBins(void) const {
    string fname = this->saveFilenames[7];
    int i1, i2;
    
    ofstream myfile (fname);
    if(!myfile.good()) {
        cerr << "Error in holePropagation::saveLabFrameSpinCorrelationBins writing to " << fname << endl << endl;
        return false;
    }
    
    for(i1=0; i1<this->nSites; i1++) {
        for(i2=0; i2<this->nSites; i2++) {
            myfile << setprecision(30) << this->LabFrameSpinCorrelationBins[i1][i2] << "\t";
        }
        myfile << endl;
    }
    myfile.close();
    
    return true;
}


// save hole frame spin correlation bins
bool holePropagation::saveHoleFrameSpinCorrelationBins(void) const {
    string fname = this->saveFilenames[8];
    int i1, i2;
    
    ofstream myfile (fname);
    if(!myfile.good()) {
        cerr << "Error in holePropagation::saveHoleFrameSpinCorrelationBins writing to " << fname << endl << endl;
        return false;
    }
    
    for(i1=0; i1<this->nSites; i1++) {
        for(i2=0; i2<this->nSites; i2++) {
            myfile << setprecision(30) << this->HoleFrameSpinCorrelationBins[i1][i2] << "\t";
        }
        myfile << endl;
    }
    myfile.close();
    
    return true;
}


// updates the file storing the return probability and the RMS distance of the hole at each backup
bool holePropagation::save_returnProb_and_RMSdistance(void) const {
    string fname = this->saveFilenames[2];
    ofstream myfile;
    
    // If the file did not exist, we create a new one.
    if(!checkIfFileExists(fname)) {
        myfile.open(fname);
        if(!myfile.good()) {
            cerr << "Error writing to " << fname << endl;
            return false;
        }
        myfile << "# \t #runs \t P_return \t RMS_distance" << endl;
        return true;
    }
    
    // If it did exist, we append a new line to it
    myfile.open(fname, ios_base::app);
    double returnProb = this->getReturnProbability();
    double RMS = this->getRMSdistance();
    if(!myfile.good()) {
        cerr << "Error writing to " << fname << endl;
        cout << "Error writing to " << fname << endl;
        cout << "Return probability: \t" << returnProb << endl;
        cout << "RMS distance: \t" << RMS << endl;
        cout << endl << endl << endl << endl;
        return false;
    }
    
    
    myfile << this->nRuns << "\t" << returnProb << "\t" << RMS << endl;
    myfile.close();
    
    return true;
}


// loads all three bins, using the three saving functions below
bool holePropagation::load(void) {
    // If the parameters in the file do not equal to those in '*this', or if the file does not exist, we do nothing
    // This is not necessarily an error, could be that the simulation has not been run yet.
    if(!this->loadParameters()) {
        cout << "Could not load parameter file " << this->saveFilenames[1] << endl;
        return false;
    }
    
    // Check if files that we want to load data from indeed exist
    // Otherwise do not start loading data
    if(!checkIfFileExists(saveFilenames[6])) {
        cerr << "Error in holePropagation::load. Could not load file " << saveFilenames[6] << endl;
        return false;
    }
    if(this->calcSpinCorrelationsQ) {
        if(!checkIfFileExists(saveFilenames[7])) {
            cerr << "Error in holePropagation::load. Could not load file " << saveFilenames[7] << endl;
            return false;
        }
        if(!checkIfFileExists(saveFilenames[8])) {
            cerr << "Error in holePropagation::load. Could not load file " << saveFilenames[8] << endl;
            return false;
        }
    }
    
    // If everything looks good, we update our bins
    if(!this->loadTransitionProbabilityBins()) {
        cerr << "Error in holePropagation::load. Could not load file " << saveFilenames[6] << endl;
        return false;
    }
    
    if(this->calcSpinCorrelationsQ) {
        // Load data from files
        // If something goes wrong beyong this point, the bins may be corrupted, since some of them were already updated,
        // whereas the outhers have not. Therefore, we clear all bins and restart the simulation.
        if(!this->loadLabFrameSpinCorrelationBins()) {
            cerr << "**********************************************************************************************" << endl;
            cerr << "Catastrophic error in holePropagation::load. Could not load file " << saveFilenames[7] << endl;
            cerr << "Bins could be corrupted. We delete all bins and restart simulation." << endl;
            cerr << "**********************************************************************************************" << endl;
            this->clearBins();
            return false;
        }
        if(!this->loadHoleFrameSpinCorrelationBins()) {
            cerr << "**********************************************************************************************" << endl;
            cerr << "Catastrophic error in holePropagation::load. Could not load file " << saveFilenames[8] << endl;
            cerr << "Bins could be corrupted. We delete all bins and restart simulation." << endl;
            cerr << "**********************************************************************************************" << endl;
            this->clearBins();
            return false;
        }
    }
    
    // If everything went well, we return true
    return true;
}   // Tested and works


// load transition probability bins
bool holePropagation::loadTransitionProbabilityBins(void) {
    bool fileExists;
    long double binValue;
    
    ifstream fileBins(this->saveFilenames[6]);
    if(!fileBins.good())
        return false;
    
    // Load data: first element of fileBins contains the number of threads performed so far.
    fileBins >> binValue;
    this->nRuns = (int) binValue;
    
    // All later elements correspond to binsTotal.
    int i = 0;
    while((fileBins >> binValue) && (i<this->nSites)) {
        transProbsBins[i] = binValue;
        i++;
    }
    
    return true;
}   // Tested and works


// load lab frame spin correlation bins
bool holePropagation::loadLabFrameSpinCorrelationBins(void) {
    // open file
    ifstream myfile(this->saveFilenames[7]);
    if(!myfile.good()) {
        return false;
    }
    
    // start reading
    int i1, i2;  // indices
    string line; // one line of the text file
    
    // read lines one-by-one from the file
    for(i1=0; i1<this->nSites; i1++) {
        getline(myfile, line);
        istringstream stream(line);
        
        // read elements of the line one-by-one
        for(i2=0; i2<this->nSites; i2++) {
            if(!(stream >> (this->LabFrameSpinCorrelationBins[i1][i2]))) {
                return false;
            }
        }
    }
    
    return true;
}   // Tested and works


// load hole frame spin correlation bins
bool holePropagation::loadHoleFrameSpinCorrelationBins(void) {
    // open file
    ifstream myfile(this->saveFilenames[8]);
    if(!myfile.good()) {
        return false;
    }
    
    // start reading
    int i1, i2;  // indices
    string line; // one line of the text file
    
    // read lines one-by-one from the file
    for(i1=0; i1<this->nSites; i1++) {
        if(!getline(myfile, line)) {
            return false;
        }
        istringstream stream(line);
        
        // read elements of the line one-by-one
        for(i2=0; i2<this->nSites; i2++) {
            if(!(stream >> this->HoleFrameSpinCorrelationBins[i1][i2])) {
                return false;
            }
        }
    }
    
    return true;
}   // Tested and works


// loads parameters and returns true only if they match those of '*this'
bool holePropagation::loadParameters(void) {
    string fname = saveFilenames[1];
    string line;
    ifstream myfile(fname);
    
    // Data we read from the file
    double tHole_FromFile;
    int nPaths_FromFile;
    int LATTICE_SIZE_FromFile;
    int NUMBER_OF_SPINS_FromFile;
    
    // Check if file exists and could be opened
    if(!myfile.good()) {
        return false;
    }
    
    // First line is just comments
    if(!getline(myfile, line))
        return false;
    
    // Second line contains the data in the following order: time, nPaths, LATTICE_SIZE, NUMBER_OF_SPINS
    if(!getline(myfile, line))
        return false;
    
    // Read data
    istringstream stream(line);
    stream >> tHole_FromFile;
    stream >> nPaths_FromFile;
    stream >> LATTICE_SIZE_FromFile;
    stream >> NUMBER_OF_SPINS_FromFile;
    
    // Check if it is the same as the one we have in *this
    // If any of the data members are incorrect, return false
    if(tHole_FromFile != this->tHole)
        return false;
    if(nPaths_FromFile != this->nPaths)
        return false;
    if(LATTICE_SIZE_FromFile != LATTICE_SIZE)
        return false;
    if(NUMBER_OF_SPINS_FromFile != this->NUMBER_OF_SPINS)
        return false;
    
    myfile.close();
    
    // If everything went well, return true
    return true;
}   // Tested and works


// Check if the file already exists
bool checkIfFileExists(string fname) {
    ifstream myfile(fname);
    return myfile.good();
} // Tested and works


// save everything
bool holePropagation::save(void) const {
    bool saveOK = true;
    
    saveOK = saveOK && this->saveAllFilenames();
    saveOK = saveOK && this->saveParameters();
    saveOK = saveOK && this->saveTransitionProbabilities();
    if(this->calcSpinCorrelationsQ) {
        saveOK = saveOK && this->saveLabFrameSpinCorrelationFn();
        saveOK = saveOK && this->saveHoleFrameSpinCorrelationFn();
    }
    if(!saveOK) {
        cerr << "Error in holePropagation::save. Could not save all files." << endl;
    }
    
    return saveOK;
}

// save the contents of saveFilenames into a file
bool holePropagation::saveAllFilenames(void) const {
    string fname = this->saveFilenames[0];
    int i1;
    
    ofstream myfile (fname);
    if(!myfile.good()) {
        cerr << "Error in holePropagation::saveAllFilenames writing to " << fname << endl << endl;
        return false;
    }
    
    for(i1=0; i1<saveFilenames.size(); i1++) {
        myfile << this->saveFilenames[i1] << endl;
    }
    myfile.close();
    
    return true;
}   // Tested and works


// save tHole, nPaths, LATTICE_SIZE, NUMBER_OF_SPINS, returnProbability, RMS distance
bool holePropagation::saveParameters(void) const {
    string fname = this->saveFilenames[1];
    ofstream myfile (fname);
    if(!myfile) {
        cerr << "Error in holePropagation::saveParameters writing to " << fname << endl << endl;
        return false;
    }
    
    myfile << "# time \t #paths \t LATTICE_SIZE \t NUMBER_OF_SPINS \t #runs_finished \t return_probability \t RMS_distance" << endl;
    myfile << this->tHole << "\t";
    myfile << this->nPaths << "\t";
    myfile << LATTICE_SIZE << "\t";
    myfile << this->NUMBER_OF_SPINS << "\t";
    myfile << this->nRuns << "\t";
    myfile << this->getReturnProbability() << "\t";
    myfile << this->getRMSdistance() << "\t";
    myfile << endl;
    
    myfile.close();
    
    return true;    
}   // Tested and works


// save transition probabilities of the hole from the origin to each site
bool holePropagation::saveTransitionProbabilities(void) const {
    string fname = this->saveFilenames[3];
    ofstream myfile (fname);
    if(!myfile) {
        cerr << "Error in holePropagation::saveTransitionProbabilities writing to " << fname << endl << endl;
        return false;
    }
    
    vector<double> tP;
    tP = this->getTransitionProbabilities();
    for(int i=0; i<this->nSites; i++) {
        myfile << tP[i] << endl;
    }
    
    myfile.close();
    return true;
}


// save spin correlation function in the lab frame
bool holePropagation::saveLabFrameSpinCorrelationFn(void) const {
    string fname = this->saveFilenames[4];
    ofstream myfile (fname);
    if(!myfile) {
        cerr << "Error in holePropagation::saveLabFrameSpinCorrelationFn writing to " << fname << endl << endl;
        return false;
    }
    
    vector< vector<double> > spinCorr;
    spinCorr = this->getLabFrameSpinCorrelationFn();
    int i1, i2;
    for(i1=0; i1<this->nSites; i1++) {
        for(i2=0; i2<this->nSites; i2++) {
            myfile << spinCorr[i1][i2] << "\t";
        }
        myfile << endl;
    }
    
    myfile.close();
    return true;
}


// save spin correlation function in the hole's comoving frame
bool holePropagation::saveHoleFrameSpinCorrelationFn(void) const {
    string fname = this->saveFilenames[5];
    ofstream myfile (fname);
    if(!myfile) {
        cerr << "Error in holePropagation::saveHoleFrameSpinCorrelationFn writing to " << fname << endl << endl;
        return false;
    }
    
    vector< vector<double> > spinCorr;
    spinCorr = this->getHoleFrameSpinCorrelationFn();
    int i1, i2;
    for(i1=0; i1<this->nSites; i1++) {
        for(i2=0; i2<this->nSites; i2++) {
            myfile << spinCorr[i1][i2] << "\t";
        }
        myfile << endl;
    }
    
    myfile.close();
    return true;
}


// add the values of the bins in hP to *this, and add their nRuns values
holePropagation & holePropagation::operator+= (const holePropagation & hP) {
    bool sameParametersQ;
    
    // we only add 'hP' to '*this' if they have the same physical parameters, and they only differ by their bins
    // let us first test this
    sameParametersQ = true;
    sameParametersQ = sameParametersQ && (this->NUMBER_OF_SPINS == hP.NUMBER_OF_SPINS);
    sameParametersQ = sameParametersQ && (this->nSites          == hP.nSites);
    sameParametersQ = sameParametersQ && (this->tHole           == hP.tHole);
    sameParametersQ = sameParametersQ && (this->nPaths          == hP.nPaths);
    if(sameParametersQ == false)
        return *this;
    
    // if they really have equal parameters, then we can add the values of the bins in 'hP' to those in '*this'
    int i1, i2;
    this->nRuns += hP.nRuns;
    for(i1=0; i1<this->nSites; i1++) {
        this->transProbsBins[i1] += hP.transProbsBins[i1];
        for(i2=0; i2<this->nSites; i2++) {
            this->LabFrameSpinCorrelationBins[i1][i2]  += hP.LabFrameSpinCorrelationBins[i1][i2];
            this->HoleFrameSpinCorrelationBins[i1][i2] += hP.HoleFrameSpinCorrelationBins[i1][i2];
        }
    }
    return *this;
} // Tested and works


// assignment operator
holePropagation & holePropagation::operator=  (const holePropagation & hP) {
    int i1,i2;
    this->tHole = hP.tHole;
    this->nPaths = hP.nPaths;
    this->nRuns = hP.nRuns;
    
    for(i1=0; i1<this->nSites; i1++) {
        this->transProbsBins[i1] = hP.transProbsBins[i1];
        for(i2=0; i2<this->nSites; i2++) {
            this->LabFrameSpinCorrelationBins[i1][i2]  = hP.LabFrameSpinCorrelationBins[i1][i2];
            this->HoleFrameSpinCorrelationBins[i1][i2] = hP.HoleFrameSpinCorrelationBins[i1][i2];
        }
    }
    return *this;
} // Tested and works


// Symmetry averaging of the transition probabilities
vector<double> symmetryAverage(vector<double> probs) {
    vector<double> probsTrf(probs.size());
    vector<double> probsAvg = probs;
    int i1;
    
    // Average over reflections
    probsTrf = transformProbs(probsAvg, reflect10);
    for(i1=0; i1<probs.size(); i1++) {
        probsAvg[i1] = 0.5 * (probsAvg[i1] + probsTrf[i1]);
    }
    probsTrf = transformProbs(probsAvg, reflect01);
    for(i1=0; i1<probs.size(); i1++) {
        probsAvg[i1] = 0.5 * (probsAvg[i1] + probsTrf[i1]);
    }
    probsTrf = transformProbs(probsAvg, reflect11);
    for(i1=0; i1<probs.size(); i1++) {
        probsAvg[i1] = 0.5 * (probsAvg[i1] + probsTrf[i1]);
    }
    
    // Average over rotations
    probsTrf = transformProbs(probsAvg, rotate90);
    for(i1=0; i1<probs.size(); i1++) {
        probsAvg[i1] += probsTrf[i1];
    }
    probsTrf = transformProbs(probsTrf, rotate90);
    for(i1=0; i1<probs.size(); i1++) {
        probsAvg[i1] += probsTrf[i1];
    }
    probsTrf = transformProbs(probsTrf, rotate90);
    for(i1=0; i1<probs.size(); i1++) {
        probsAvg[i1] += probsTrf[i1];
    }
    for(i1=0; i1<probs.size(); i1++) {
        probsAvg[i1] *= 0.25;
    }
    
    return probsAvg;
} // Tested and works


// 90 degree rotation of lattice indices
int rotate90(int siteIndex) {
    int x = indexToX(siteIndex);
    int y = indexToY(siteIndex);
    
    int xNew = -y;
    int yNew = x;
    
    return coordToIndex(xNew, yNew);
} // Tested and works


// reflection of site indices to the (1,0) axis
int reflect10(int siteIndex) {
    int x = indexToX(siteIndex);
    int y = indexToY(siteIndex);
    
    int xNew = x;
    int yNew = -y;
    
    return coordToIndex(xNew, yNew);
} // Tested and works


// reflection of site indices to the (0,1) axis
int reflect01(int siteIndex) {
    int x = indexToX(siteIndex);
    int y = indexToY(siteIndex);
    
    int xNew = -x;
    int yNew = y;
    
    return coordToIndex(xNew, yNew);
} // Tested and works


// reflection of site indices to the (1,1) axis
int reflect11(int siteIndex) {
    int x = indexToX(siteIndex);
    int y = indexToY(siteIndex);
    
    int xNew = y;
    int yNew = x;
    
    return coordToIndex(xNew, yNew);
} // Tested and works


// transforms the indices of transProb
vector<double> transformProbs(vector<double> probs, int (*transformation) (int)) {
    vector<double> probsTrf(probs.size());
    
    for(int iSite=0; iSite<probs.size(); iSite++) {
        probsTrf[(*transformation)(iSite)] = probs[iSite];
    }
    
    return probsTrf;
} // Tested and works




