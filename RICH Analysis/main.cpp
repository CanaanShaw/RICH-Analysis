//
//  main.cpp
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/4/14.
//

#include "../tooLib/sharedHeader.h"
#include "weightModel.h"
#include "../tooLib/tooLib.h"
#include "../tooLib/myTree.h"

int main(int argc, const char * argv[]) {
	
	// File initialization
	
	string dataPath = "/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/data/*.root";
	myTree tree(dataPath.c_str(), "dailyData");
	
//	cout << "Reading from file: " << argv[1] << endl;
//	myTree tree(argv[1], "dailyData");
	
	tree.SetRichAddress();
	tree.SetTrackerAddress();
	tree.SetTofAddress();

	string outPath	= "/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/rec.GaussFitEllipseFromXCode.root";
	TFile * outFile = new TFile(outPath.c_str(), "RECREATE");
	
//	cout << "Writing to file: " << argv[2] << endl;
//	TFile * outFile = new TFile(argv[2], "RECREATE");
	
	outFile -> cd();
	TTree * outTree = new TTree("treeRec", "Rec Result");
	tree.MakeRecAddress(outTree);
	tree.MakeRichAddress(outTree);
	tree.MakeTrackerAddress(outTree);
	double recTheta;		// Reconstructed cherenkov angle
	outTree -> Branch("recTheta", &recTheta, "recTheta/D");
	double meanInlierDistance;
	outTree -> Branch("meanInlierDistance", &meanInlierDistance, "meanInlierDistance/D");
	double meanErr;
	outTree -> Branch("meanErr", &meanErr, "meanErr/D");
	double absMeanInlierDistance;
	outTree -> Branch("absMeanInlierDistance", &absMeanInlierDistance, "absMeanInlierDistance/D");
	double absMeanErr;
	outTree -> Branch("absMeanErr", &absMeanErr, "absMeanErr/D");
	
	TFile * mapFile = new TFile("/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/RefractiveMap.root", "READ");
	TH2D * refractiveMap = (TH2D *) mapFile -> Get("a");

	// Static correction of RICH spacial position
	const double dx = 0.09;
	const double dy = -0.075;
	
	// Start processing data
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
			int crossedCenterSize = 0;
			double crossedCenterX = 0;
			double crossedCenterY = 0;
			for (int iHit = 0; iHit < tree.nHits; iHit ++) {
				
				vector < double > tempHit {tree.pointX[iHit], tree.pointY[iHit], 0};
				if (!tree.isCrossed[iHit]) {
					hits.push_back(tempHit);
				} else if(distance(center, tempHit) < 10) {
					crossedCenterSize ++;
					crossedCenterX += tree.pointX[iHit];
					crossedCenterY += tree.pointY[iHit];
				}
			}
			
//			if (crossedCenterSize > 3) {
//				crossedCenterX /= crossedCenterSize;
//				crossedCenterY /= crossedCenterSize;
//				cout << "(" << centerA - crossedCenterX << ",	" << centerB - crossedCenterY << ")\n";
//			}
			
			unsigned long hitsAvailable = hits.size();
			if (hitsAvailable < 5) {
				continue;
			}
			
			double maxWeight = 0.0;
			recTheta = 0.0;
			double theta = atan(distance(center, radCenter) / RichConst::aglTransmissionHeight());
			double phi = atan((centerB - radB) / (centerA - radA));
			if (centerA - radA < 0) phi += TMath::Pi();
			
			// Find the most possible Cherenkov angle theta_i
			for (double theta_i = 0.1; theta_i < 0.4; theta_i += 0.005) {
				
				// Compute parameters of ellipse
				double semiMajorAxis = 0.5 * RichConst::aglTransmissionHeight() * (tan(theta + theta_i) - tan(theta - theta_i));
				double semiMinorAxis = RichConst::aglTransmissionHeight() / cos(theta) * tan(theta_i);
				
				// Distance from radCenter to ellipseCenter from top view.
				double rad2Center = RichConst::aglTransmissionHeight() * tan(theta - theta_i) + semiMajorAxis;
				
				vector < double > ellipseCenter {radA + rad2Center * cos(phi) , radB + rad2Center * sin(phi)};
				if ((ellipseCenter[0] - centerA) / (centerA - radA) < 0) cout << "Wrong\n";
				if ((ellipseCenter[1] - centerB) / (centerB - radB) < 0) cout << "Wrong\n";
				
				// Compute the likelihood of theta_i
				double weight = 1.0;
				double hitsAvailable = hits.size();
				double inliers = 0;
				vector < double > tempInlierDistance;
				vector < double > absTempInlierDistance;
				for (int i = 0; i < hits.size(); i ++) {
					
					double hit2Rad = distance(radCenter, hits[i]);
					double hit2Center = distance(ellipseCenter, hits[i]);
//					cout << fabs(distance(ellipseCenter, center)) << endl;
					
					if (hit2Center > 22.0) {
						hitsAvailable --;
						continue;
					}
					
					// Angle phi being radCenter-ellipseCenter-hit
					double cosPhi = (square(rad2Center) + square(hit2Center) - square(hit2Rad)) / 2.0 / rad2Center / hit2Center;
					
					// rExpected being the expected distance of the point on the ellipse to the center, which is on the direction of photon hit to center.
					double rExpected = sqrt(square(semiMajorAxis * cosPhi) + square(semiMinorAxis * sqrt(1.0 - cosPhi * cosPhi)));
					weight *= weightModel(hit2Center, rExpected);
					if (fabs(rExpected - hit2Center) < 3) {
						tempInlierDistance.push_back(rExpected - hit2Center);
						absTempInlierDistance.push_back(fabs(rExpected - hit2Center));
						inliers ++;
					}
				}
				
				if (inliers < 4) continue;
				if (hitsAvailable * 1.0 / hits.size() < 0.5 || hitsAvailable < 4) continue;
				
				if (maxWeight < weight) {
					maxWeight = weight;
					recTheta = theta_i;
					
					double mean1 = 0;
					double err1 = 0;
					for (double ii : tempInlierDistance) {
						mean1 += ii;
					}
					
					meanInlierDistance = mean1 / inliers;
					for (double ii : tempInlierDistance) {
						err1 += square(ii - meanInlierDistance);
					}
					
					meanErr = sqrt(err1) / inliers;
					
					double mean2 = 0;
					double err2 = 0;
					for (double ii : absTempInlierDistance) {
						mean2 += ii;
					}
					
					absMeanInlierDistance = mean2 / inliers;
					for (double ii : absTempInlierDistance) {
						err2 += square(ii - absMeanInlierDistance);
					}
					
					absMeanErr = sqrt(err2) / inliers;
				}
			}
			
			if (recTheta == 0) continue;
			
			
			

			double n2 = refractiveMap -> GetBinContent(refractiveMap -> FindBin(radA, radB));
//			double n2 = 1.055;
			tree.recBeta = 1.0 / n2 / cos(recTheta);
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
