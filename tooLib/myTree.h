//
//  myTree.h
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/6/7.
//

#ifndef myTree_h
#define myTree_h

#include "TChain.h"
#define MAXIMUM 999

class myTree {

private:
public:

    TChain * pTree;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    RichHitR
    int nHits;

    float     pointX[MAXIMUM];
    float     pointY[MAXIMUM];
    float     pointZ[MAXIMUM];

    bool     isCrossed[MAXIMUM];
    float     Npe[MAXIMUM];
    float     Cpe[MAXIMUM];
    int     Channel[MAXIMUM];

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    RichRingR
    int     nRings;
    
    float     richBeta[MAXIMUM];
    bool     isNaF[MAXIMUM];
    int     tileIndex[MAXIMUM];

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    Tracker
    int     nTrack;

    double    trPointX;
    double    trPointY;
    double    trPointZ;
    double     trRadX;
    double    trRadY;
    double     trRadZ;

    double    trRigidity;
    float    trCharge;
    float    trInnerCharge;
    float    trLayerCharge[9];
    double    trChi2X;
    double    trChi2Y;
    double    trTheta;
    double    trPhi;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    TOF
    float    tofQL[4];
    float    tofX[4];
    float    tofY[4];
    float    tofZ[4];
    double    tofRadX;
    double    tofRadY;
    double    tofRadZ;
    double    tofTheta;
    double    tofPhi;
    float    tofBeta;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    Rec

    double recR;
    double recBeta;
    double recMass;

    myTree (const char * path) {

        pTree = new TChain("dailyData");
        pTree -> Add(path);
        SetRichAddress(pTree);
        SetTrackerAddress(pTree);
        SetTofAddress(pTree);
    }

    void SetRichAddress(TChain * inTree) {

        inTree -> SetBranchAddress("nHits",            &nHits);
        inTree -> SetBranchAddress("pointX",        pointX);
        inTree -> SetBranchAddress("pointY",        pointY);
        inTree -> SetBranchAddress("pointZ",        pointZ);
        inTree -> SetBranchAddress("isCrossed",        isCrossed);
        inTree -> SetBranchAddress("Npe",             Npe);
        inTree -> SetBranchAddress("Cpe",             Cpe);
        inTree -> SetBranchAddress("Channel",         Channel);

        inTree -> SetBranchAddress("nRings",         &nRings);
        inTree -> SetBranchAddress("richBeta",         richBeta);
        inTree -> SetBranchAddress("isNaF",         isNaF);
        inTree -> SetBranchAddress("tileIndex",     tileIndex);
    }
    void SetTrackerAddress(TChain * inTree) {

        inTree -> SetBranchAddress("nTrack",         &nTrack);
        inTree -> SetBranchAddress("trPointX",         &trPointX);
        inTree -> SetBranchAddress("trPointY",         &trPointY);
        inTree -> SetBranchAddress("trPointZ",         &trPointZ);
        inTree -> SetBranchAddress("trRadX",         &trRadX);
        inTree -> SetBranchAddress("trRadY",         &trRadY);
        inTree -> SetBranchAddress("trRadZ",         &trRadZ);
        inTree -> SetBranchAddress("trRigidity",     &trRigidity);
        inTree -> SetBranchAddress("trCharge",         &trCharge);
        inTree -> SetBranchAddress("trInnerCharge", &trInnerCharge);
        inTree -> SetBranchAddress("trLayerCharge", trLayerCharge);
        inTree -> SetBranchAddress("trChi2X",         &trChi2X);
        inTree -> SetBranchAddress("trChi2Y",         &trChi2Y);
        inTree -> SetBranchAddress("trTheta",         &trTheta);
        inTree -> SetBranchAddress("trPhi",         &trPhi);
    }
    void SetTofAddress(TChain * inTree) {

        inTree -> SetBranchAddress("tofQL",         tofQL);
        inTree -> SetBranchAddress("tofX",            tofX);
        inTree -> SetBranchAddress("tofY",            tofY);
        inTree -> SetBranchAddress("tofZ",            tofZ);
        inTree -> SetBranchAddress("tofRadX",        &tofRadX);
        inTree -> SetBranchAddress("tofRadY",        &tofRadY);
        inTree -> SetBranchAddress("tofRadZ",        &tofRadZ);
        inTree -> SetBranchAddress("tofTheta",        &tofTheta);
        inTree -> SetBranchAddress("tofPhi",        &tofPhi);
        inTree -> SetBranchAddress("tofBeta",        &tofBeta);
    }
    void SetRecAddress(TChain * inTree) {

        inTree -> SetBranchAddress("recR",             &recR);
        inTree -> SetBranchAddress("recBeta",         &recBeta);
        inTree -> SetBranchAddress("recMass",         &recMass);
    }

    void MakeRichAddress(TTree * inTree) {

        inTree -> Branch("nHits",            &nHits,         "nHits/I");
        inTree -> Branch("pointX",            pointX,         "pointX[nHits]/F");
        inTree -> Branch("pointY",            pointY,         "pointY[nHits]/F");
        inTree -> Branch("pointZ",            pointZ,         "pointZ[nHits]/F");
        inTree -> Branch("isCrossed",        isCrossed,         "isCrossed[nHits]/O");
        inTree -> Branch("Npe",             Npe,             "Npe[nHits]/F");
        inTree -> Branch("Cpe",             Cpe,             "Cpe[nHits]/F");
        inTree -> Branch("Channel",         Channel,         "Channel[nHits]/I");

        inTree -> Branch("nRings",             &nRings,         "nRings/I");
        inTree -> Branch("richBeta",         richBeta,         "richBeta[nRings]/F");
        inTree -> Branch("isNaF",             isNaF,             "isNaF[nRings]/O");
        inTree -> Branch("tileIndex",         tileIndex,         "tileIndex[nRings]/I");
    }
    void MakeTrackerAddress(TTree * inTree) {

        inTree -> Branch("nTrack",             &nTrack,         "nTrack/I");
        inTree -> Branch("trPointX",         &trPointX,         "trPointX/D");
        inTree -> Branch("trPointY",         &trPointY,         "trPointY/D");
        inTree -> Branch("trPointZ",         &trPointZ,         "trPointZ/D");
        inTree -> Branch("trRadX",             &trRadX,         "trRadX/D");
        inTree -> Branch("trRadY",             &trRadY,         "trRadY/D");
        inTree -> Branch("trRadZ",             &trRadZ,         "trRadZ/D");
        inTree -> Branch("trRigidity",         &trRigidity,     "trRigidity/D");
        inTree -> Branch("trCharge",         &trCharge,         "trCharge/F");
        inTree -> Branch("trInnerCharge",     &trInnerCharge, "trInnerCharge/F");
        inTree -> Branch("trLayerCharge",     trLayerCharge,     "trLayerCharge[9]/F");
        inTree -> Branch("trChi2X",         &trChi2X,         "trChi2X/D");
        inTree -> Branch("trChi2Y",         &trChi2Y,         "trChi2Y/D");
        inTree -> Branch("trTheta",         &trTheta,         "trTheta/D");
        inTree -> Branch("trPhi",             &trPhi,         "trPhi/D");
    }
    void MakeTofAddress(TTree * inTree) {

        inTree -> Branch("tofQL",             tofQL,             "tofQL[4]/F");
        inTree -> Branch("tofX",            tofX,             "tofX[4]F");
        inTree -> Branch("tofY",            tofY,             "tofY[4]F");
        inTree -> Branch("tofZ",            tofZ,             "tofZ[4]F");
        inTree -> Branch("tofRadX",            &tofRadX,         "tofRadX/D");
        inTree -> Branch("tofRadY",            &tofRadY,         "tofRadY/D");
        inTree -> Branch("tofRadZ",            &tofRadZ,         "tofRadZ/D");
        inTree -> Branch("tofTheta",        &tofTheta,         "tofTheta/D");
        inTree -> Branch("tofPhi",            &tofPhi,         "tofPhi/D");
        inTree -> Branch("tofBeta",            &tofBeta,         "tofBeta/F");
    }
    void MakeRecAddress(TTree * inTree) {

        inTree -> Branch("recR",             &recR,             "recR/D");
        inTree -> Branch("recBeta",         &recBeta,         "recBeta/D");
        inTree -> Branch("recMass",         &recMass,         "recMass/D");
    }

    long long  GetEntries() { return pTree -> GetEntries(); }
    void GetEntry(int i) { pTree -> GetEntry(i); }
    bool Select() {

        bool cut[15];
        for (int i = 0; i < 15; i ++) {

            cut[i] = 1;
        }

        cut[0] = (trTheta < 0.12);
        cut[1] = (trChi2X < 10);
        cut[2] = (trChi2Y < 10);
        cut[3] = (trRigidity > 0);
        cut[4] = (fabs(trRadX) > 17.0 || fabs(trRadY) > 17.0);
        cut[5] = (tofBeta > 0);
        cut[6] = (trInnerCharge > 0);
        cut[7] = (nHits > 3);

        for (int i = 0; i < 15; i ++) {
            
            if (!cut[i]) return false;
        }

        return true;
    }

};

#endif /* myTree_h */
