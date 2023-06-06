//
//  main.cpp
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/4/14.
//

#include <iostream>
#include "weightModel.h"
#include "TTree.h"
#include "tooLib.h"

#include <cstring>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TH3.h"
#include "TMath.h"

#include "TString.h"
#include <vector>
#include "TRandom2.h"
#include "TH2.h"

int main(int argc, const char * argv[]) {
    // insert code here...
	string ss = "/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/data/*.root";
	myTree tree(ss.c_str());

	TString outFilename= "GaussFitEllipseFromXCode";
	TString path 	= "/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/rec.";
			path	+= outFilename;
			path	+= ".root";
	TFile * outFile 	= new TFile(path, "RECREATE");

	outFile -> cd();
	TTree * outTree = new TTree("treeRec", "Rec Result");
	tree.MakeRecAddress(outTree);
	tree.MakeRichAddress(outTree);

	TFile * map2 = new TFile("/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/RefractiveMap.root", "READ");
	TH2D * nMap2 = (TH2D * )map2 -> Get("a");

	const double dx = 0.09;
	const double dy = -0.075;
	for (int iEvent = 0; iEvent < tree.GetEntries(); iEvent ++) {
		
		tree.GetEntry(iEvent);
		if (iEvent % 10000 == 0)
			cout << "Current: " << iEvent << "	out of " << tree.GetEntries() << endl;

		if (tree.Select()) {

			double centerA = tree.trPointX + dx;
			double centerB = tree.trPointY + dy;
			vector < double > center {centerA, centerB};
			
			double radA = tree.trRadX + dx;
			double radB = tree.trRadY + dy;
			vector < double > radCenter {radA, radB};

			vector < vector < double > > hits;
			for (int iHit = 0; iHit < tree.nHits; iHit ++) {
				if (!tree.isCrossed[iHit]) {

					vector < double > tempHit {tree.pointX[iHit] - centerA, tree.pointY[iHit] - centerB, 0};
					hits.push_back(spinY(spinZ(tempHit, tree.trPhi), tree.trTheta));
				}
			}

			unsigned long hitsAvailable = hits.size();
			if (hitsAvailable < 5) {
				continue;
			}

			vector < double > rHits;
			for (int i = 0; i < hits.size(); i ++) {
				
				double r_i = distance3D(vector < double > {0, 0, 0}, hits[i]) * cos(tree.trTheta);
				if (r_i > 23 || r_i < 5) continue;
				rHits.push_back(r_i);
			}

			if (rHits.size() < 3 || rHits.size() * 1.0 / hits.size() < 0.5) continue;

			double maxWeight = 0;
			double rInitialGuess = -1.0;

			for (double r_i = 5.0; r_i < 19.0; r_i += 0.2) {
				
				double weight_i = computeWeight(r_i, rHits);
				if (maxWeight < weight_i) {
					maxWeight = weight_i;
					rInitialGuess = r_i;
				}
			}

			double maxWeightRefine = 0;
			double rRefine = -1.0;

			vector < double > rHitsClean;
			for (int rIndex = 0; rIndex < rHits.size(); rIndex ++) {

				if (fabs(rHits[rIndex] - rInitialGuess) < 1) rHitsClean.push_back(rHits[rIndex]);
			}

			if (rHitsClean.size() < 6 || rHitsClean.size() * 1.0 / rHits.size() < 0.5) continue;

			// This Range affects the number of final succesful reconstruction events.
			for (double r_i = max < double > (rInitialGuess - 0.2, 5.0); r_i < min < double > (rInitialGuess + 0.2, 19.0); r_i += 0.01) {

				double weight_i = computeWeight(r_i, rHitsClean);
				if (maxWeightRefine <= weight_i) {
					maxWeightRefine = weight_i;
					rRefine = r_i;
				}
			}

			tree.recR = rRefine;

			using namespace TMath;
			double n2 = nMap2 -> GetBinContent(nMap2 -> FindBin(radA, radB));
			double Theta = ATan(tree.recR / RichConst::NaFTransmissionHeight());
			tree.recBeta = 1.0 / n2 / cos(Theta);
			tree.recMass = tree.trRigidity * tree.trInnerCharge / tree.recBeta * sqrt(1.0 - square(tree.recBeta));

			outTree -> Fill();
		}
	}

	outFile -> cd();
	outTree -> Write();
	outFile -> Save();
	outFile -> Close();
	
	return 0;
}
