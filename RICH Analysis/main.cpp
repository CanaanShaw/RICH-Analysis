//
//  main.cpp
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/4/14.
//

#include "sharedHeader.h"
#include "weightModel.h"
#include "tooLib.h"
#include "myTree.h"

vector < double > FindCrossPoint(vector < double > hits, vector < double > ringCenter,
								 vector < double > radCenter, double theta, double phi, double z = 47.17) {
	// Intersection point of each line between photon hit and emission point, with spin plane going through ringCenter.
	
	double t = -1.0 * (hits[0] * sin(theta) * cos(phi) + hits[1] * sin(theta) * sin(phi)) / ((radCenter[0] - hits[0]) * sin(theta) * cos(phi) + (radCenter[1] - hits[1]) * sin(theta) * sin(phi) + z * cos(theta));
	
	double xOut = hits[0] + (radCenter[0] - hits[0]) * t;
	double yOut = hits[1] + (radCenter[1] - hits[1]) * t;
	double zOut = z * t;
	
	return vector < double > {xOut, yOut, zOut};
}

int main(int argc, const char * argv[]) {
	
	// File initialization
	string dataPath = "/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/data/*.root";
	myTree tree(dataPath.c_str(), "dailyData");
	tree.SetRichAddress();
	tree.SetTrackerAddress();
	tree.SetTofAddress();

	string outPath	= "/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/rec.GaussFitEllipseFromXCode.root";
	TFile * outFile = new TFile(outPath.c_str(), "RECREATE");
	outFile -> cd();
	TTree * outTree = new TTree("treeRec", "Rec Result");
	tree.MakeRecAddress(outTree);
	tree.MakeRichAddress(outTree);
	tree.MakeTrackerAddress(outTree);
	
	TFile * mapFile = new TFile("/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/RefractiveMap.root", "READ");
	TH2D * refractiveMap = (TH2D *) mapFile -> Get("a");

	// Static correction of RICH spacial position
	const double dx = 0.09;
	const double dy = -0.075;
	
	for (int iEvent = 0; iEvent < tree.GetEntries(); iEvent ++) {
		
		tree.GetEntry(iEvent);
		if (iEvent % 10000 == 0)
			cout << "Current: " << iEvent << "	out of " << tree.GetEntries() << endl;

		if (tree.Select()) {

			// Tracker extrapolated photon hit position with correction
			double centerA = tree.trPointX + dx;
			double centerB = tree.trPointY + dy;
			vector < double > center {centerA, centerB};
			
			// Tracker extrapolated radiation emission point with correction
			double radA = tree.trRadX + dx;
			double radB = tree.trRadY + dy;
			vector < double > radCenter {radA, radB};

			// Photon hits position without ones crossed by particles
			vector < vector < double > > hits;
			for (int iHit = 0; iHit < tree.nHits; iHit ++) {
				if (!tree.isCrossed[iHit]) {
					vector < double > tempHit {tree.pointX[iHit], tree.pointY[iHit], 0};
					hits.push_back(tempHit);
					
//					hits.push_back(FindCrossPoint(tempHit, center, radCenter, tree.trTheta, tree.trPhi));
				}
			}

			unsigned long hitsAvailable = hits.size();
			if (hitsAvailable < 5) {
				continue;
			}

			// Distances of each photon hit to ring center
			vector < double > rHits;
			for (int i = 0; i < hits.size(); i ++) {
				
				double r_i = distance3D(vector < double > {centerA, centerB, 0}, hits[i]) * cos(tree.trTheta);
				if (r_i > 30 || r_i < 5) continue;
				rHits.push_back(r_i);
			}

			if (rHits.size() < 3 || rHits.size() * 1.0 / hits.size() < 0.5) continue;

			double maxWeight = 0;
			double rInitialGuess = -1.0;

			for (double r_i = 5.0; r_i < 23.0; r_i += 0.2) {
				
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

				if (fabs(rHits[rIndex] - rInitialGuess) < 3) rHitsClean.push_back(rHits[rIndex]);
			}

			if (rHitsClean.size() < 6 || rHitsClean.size() * 1.0 / rHits.size() < 0.5) continue;

			// This Range affects the number of final succesful reconstruction events.
			for (double r_i = max < double > (rInitialGuess - 0.2, 5.0); r_i < min < double > (rInitialGuess + 0.2, 23.0); r_i += 0.01) {

				double weight_i = computeWeight(r_i, rHitsClean);
				if (maxWeightRefine <= weight_i) {
					maxWeightRefine = weight_i;
					rRefine = r_i;
				}
			}

			tree.recR = rRefine;

			using namespace TMath;
			double n2 = refractiveMap -> GetBinContent(refractiveMap -> FindBin(radA, radB));
			double Theta = ATan(tree.recR / RichConst::NaFTransmissionHeight());
			tree.recBeta = 1.0 / n2 / cos(Theta);
			tree.recMass = tree.trRigidity * tree.trInnerCharge / tree.recBeta * sqrt(1.0 - square(tree.recBeta));

			// Correction for position changes from light guide to PMT pixel?
			
			outTree -> Fill();
		}
	}

	outFile -> cd();
	outTree -> Write();
	outFile -> Save();
	outFile -> Close();
	
	return 0;
}
