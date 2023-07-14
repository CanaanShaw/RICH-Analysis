//
//  myTree.h
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/6/7.
//

#ifndef myTree_h
#define myTree_h

#include "sharedHeader.h"

class myTree {

private:
public:

    TChain * pTree;
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    AMS
    unsigned int AMSTime;

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
    double recTheta;

    myTree (const char * path, const char * treeName) {

        pTree = new TChain(treeName);
        pTree -> Add(path);
    }
    
    void SetAMSAddress() {
        
        pTree -> SetBranchAddress("AMSTime", &AMSTime);
    }

    void SetRichAddress() {

        pTree -> SetBranchAddress("nHits",            &nHits);
        pTree -> SetBranchAddress("pointX",        pointX);
        pTree -> SetBranchAddress("pointY",        pointY);
        pTree -> SetBranchAddress("pointZ",        pointZ);
        pTree -> SetBranchAddress("isCrossed",        isCrossed);
        pTree -> SetBranchAddress("Npe",             Npe);
        pTree -> SetBranchAddress("Cpe",             Cpe);
        pTree -> SetBranchAddress("Channel",         Channel);

        pTree -> SetBranchAddress("nRings",         &nRings);
        pTree -> SetBranchAddress("richBeta",         richBeta);
        pTree -> SetBranchAddress("isNaF",         isNaF);
        pTree -> SetBranchAddress("tileIndex",     tileIndex);
    }
    void SetTrackerAddress() {

        pTree -> SetBranchAddress("nTrack",         &nTrack);
        pTree -> SetBranchAddress("trPointX",         &trPointX);
        pTree -> SetBranchAddress("trPointY",         &trPointY);
        pTree -> SetBranchAddress("trPointZ",         &trPointZ);
        pTree -> SetBranchAddress("trRadX",         &trRadX);
        pTree -> SetBranchAddress("trRadY",         &trRadY);
        pTree -> SetBranchAddress("trRadZ",         &trRadZ);
        pTree -> SetBranchAddress("trRigidity",     &trRigidity);
        pTree -> SetBranchAddress("trCharge",         &trCharge);
        pTree -> SetBranchAddress("trInnerCharge", &trInnerCharge);
        pTree -> SetBranchAddress("trLayerCharge", trLayerCharge);
        pTree -> SetBranchAddress("trChi2X",         &trChi2X);
        pTree -> SetBranchAddress("trChi2Y",         &trChi2Y);
        pTree -> SetBranchAddress("trTheta",         &trTheta);
        pTree -> SetBranchAddress("trPhi",         &trPhi);
    }
    void SetTofAddress() {

        pTree -> SetBranchAddress("tofQL",         tofQL);
        pTree -> SetBranchAddress("tofX",            tofX);
        pTree -> SetBranchAddress("tofY",            tofY);
        pTree -> SetBranchAddress("tofZ",            tofZ);
        pTree -> SetBranchAddress("tofRadX",        &tofRadX);
        pTree -> SetBranchAddress("tofRadY",        &tofRadY);
        pTree -> SetBranchAddress("tofRadZ",        &tofRadZ);
        pTree -> SetBranchAddress("tofTheta",        &tofTheta);
        pTree -> SetBranchAddress("tofPhi",        &tofPhi);
        pTree -> SetBranchAddress("tofBeta",        &tofBeta);
    }
    void SetRecAddress() {

        pTree -> SetBranchAddress("recR",             &recR);
        pTree -> SetBranchAddress("recTheta",        &recTheta);
        pTree -> SetBranchAddress("recBeta",         &recBeta);
        pTree -> SetBranchAddress("recMass",         &recMass);
    }
    
    void MakeAMSAddress(TTree * inTree) {
        inTree -> Branch("AMSTime", &AMSTime, "AMSTime/i");
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

        inTree -> Branch("recR",            &recR,              "recR/D");
        inTree -> Branch("recTheta",        &recTheta,          "recTheta/D");
        inTree -> Branch("recBeta",         &recBeta,           "recBeta/D");
        inTree -> Branch("recMass",         &recMass,           "recMass/D");
    }
    
    void MakeRichAddress(TChain * inTree) {

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
    void MakeTrackerAddress(TChain * inTree) {

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
    void MakeTofAddress(TChain * inTree) {

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
    void MakeRecAddress(TChain * inTree) {

        inTree -> Branch("recR",             &recR,             "recR/D");
        inTree -> Branch("recTheta",        &recTheta,          "recTheta/D");
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
        cut[5] = (trInnerCharge > 0);
        cut[6] = (nHits > 3);

        for (int i = 0; i < 15; i ++) {
            
            if (!cut[i]) return false;
        }

        return true;
    }

};

#endif /* myTree_h */
