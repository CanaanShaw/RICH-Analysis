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
		if (iEvent > 1000000) break;

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
				}
			}

			unsigned long hitsAvailable = hits.size();
			if (hitsAvailable < 5) {
				continue;
			}
			
			double maxWeight = 0.0;
			double theta_c = 0.0;
			
			// Find the most possible Cherenkov angle theta_i
			for (double theta_i = 0.2; theta_i < 0.4; theta_i += 0.01) {
				
				double semiMajorAxis = 0.5 * RichConst::aglTransmissionHeight() * (tan(tree.trTheta + theta_i) - tan(tree.trTheta - theta_i));
				double semiMinorAxis = RichConst::aglTransmissionHeight() / cos(tree.trTheta) * tan(theta_i);
				
				// Compute the center of the ellipse
//				vector < double > ellipseCenter {centerA + semiMajorAxis * cos(tree.trPhi), centerB + semiMajorAxis * sin(tree.trPhi)};
				double rad2Center = RichConst::aglTransmissionHeight() * tan(tree.trTheta - theta_i) + semiMajorAxis;
				double phi = atan((centerB - radB) / (centerA - radA));
				if (centerA - radA < 0) phi += TMath::Pi();
				vector < double > ellipseCenter {radA + rad2Center * cos(phi) , radB + rad2Center * sin(phi)};
//				cout << "(" << radA << ", " << radB << "), " << "(" << centerA << ", " << centerB << "), " << "(" << ellipseCenter[0] << ", " << ellipseCenter[1] << ")\n";
//				cout << "(" << radA << ", " << radB << "), " << "(" << centerA << ", " << centerB << "), " << tree.trPhi << endl;
//				double rad2Center = distance(radCenter, ellipseCenter);
//				cout << (ellipseCenter[0] - centerA) / (centerA - radA) << endl;
//				cout << rad2Center - distance(radCenter, center) << endl;
				
				// Compute the likelihood of theta_i
				double weight = 1.0;
				for (int i = 0; i < hits.size(); i ++) {
					
					double hit2Rad = distance(radCenter, hits[i]);
					double hit2Center = distance(ellipseCenter, hits[i]);
//					cout << fabs(distance(ellipseCenter, center)) << endl;
					
					if (hit2Center > 25.0) continue;
					// Angle phi being radCenter-ellipseCenter-hit
					double cosPhi = (square(rad2Center) + square(hit2Center) - square(hit2Rad)) / 2.0 / rad2Center / hit2Center;
					
					// rExpected being the expected distance of the point on the ellipse to the center, which is on the direction of photon hit to center.
					double rExpected = sqrt(square(semiMajorAxis * cosPhi) + square(semiMinorAxis * sqrt(1.0 - cosPhi * cosPhi)));
					weight *= computeWeight2(hit2Center, rExpected);
				}
				
				if (maxWeight < weight) {
					maxWeight = weight;
					theta_c = theta_i;
				}
			}

			using namespace TMath;
			double n2 = refractiveMap -> GetBinContent(refractiveMap -> FindBin(radA, radB));
			n2 = 1.055;
//			double Theta = ATan(tree.recR / RichConst::aglTransmissionHeight());
//			tree.recBeta = 1.0 / n2 / cos(Theta);
			tree.recBeta = 1.0 / n2 / cos(theta_c);
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
