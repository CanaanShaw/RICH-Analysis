//
//  main.cpp
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/4/14.
//

#include "../tooLib/sharedHeader.h"
#include "weightModel.h"
#include "recEllipse.h"
#include "../tooLib/tooLib.h"
#include "../tooLib/myTree.h"
#include "defs.h"

#define debugFlag false

int main(int argc, const char * argv[]) {
	
	// File initialization
#ifdef LOCAL_TEST
	
	string dataPath = "/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/data/*.root";
	myTree tree(dataPath.c_str(), "dailyData");
	
	string outPath	= "/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/rec.GaussFitEllipseFromXCode.root";
	TFile * outFile = new TFile(outPath.c_str(), "RECREATE");
	
	TFile * mapFile = new TFile("/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/RefractiveMap.root", "READ");
	TH2D * refractiveMap = (TH2D *) mapFile -> Get("a");
	
	TFile * modelFile = new TFile("/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/weightModel", "READ");
	TH1D * weightModelMap = (TH1D *) modelFile -> Get("weightModel");
#else
	
	cout << "Reading from file: " << argv[1] << endl;
	myTree tree(argv[1], "dailyData");
	
	cout << "Writing to file: " << argv[2] << endl;
	TFile * outFile = new TFile(argv[2], "RECREATE");
	
	TFile * mapFile = new TFile("/afs/cern.ch/user/j/jianan/private/betaRec/RefractiveMap.root", "READ");
	TH2D * refractiveMap = (TH2D *) mapFile -> Get("a");
	
	TFile * modelFile = new TFile("/afs/cern.ch/user/j/jianan/private/betaRec/weightModel", "READ");
	TH1D * weightModelMap = (TH1D *) modelFile -> Get("weightModel");
#endif
	
	tree.SetRichAddress();
	tree.SetTrackerAddress();
	tree.SetTofAddress();

	outFile -> cd();
	TTree * outTree = new TTree("treeRec", "Rec Result");
	tree.MakeRecAddress(outTree);
	tree.MakeRichAddress(outTree);
	tree.MakeTrackerAddress(outTree);
	double recTheta;		// Reconstructed cherenkov angle
	outTree -> Branch("recTheta", &recTheta, "recTheta/D");
	int inliers;
	outTree -> Branch("inliers", &inliers, "inliers/I");
	double inlierDistance[MAXIMUM];
	outTree -> Branch("inlierDistance", inlierDistance, "inlierDistance[inliers]/D");
	
	

	// Static correction of RICH spacial position
	const double dx = 0.09;
	const double dy = -0.075;
	
	// Start processing data
	for (int iEvent = 0; iEvent < tree.GetEntries(); iEvent ++) {
		
		tree.GetEntry(iEvent);
		if (iEvent % 10000 == 0)
			cout << "Current: " << iEvent << "	out of " << tree.GetEntries() << endl;
		if (iEvent > 500000) break;

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
//			int crossedCenterSize = 0;
//			double crossedCenterX = 0;
//			double crossedCenterY = 0;
			for (int iHit = 0; iHit < tree.nHits; iHit ++) {
				
				vector < double > tempHit {tree.pointX[iHit], tree.pointY[iHit], 0};
				if (!tree.isCrossed[iHit]) {
					hits.push_back(tempHit);
				}
//				else if(distance(center, tempHit) < 10) {
//					crossedCenterSize ++;
//					crossedCenterX += tree.pointX[iHit];
//					crossedCenterY += tree.pointY[iHit];
//				}
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
			
			recTheta = 0.0;
			inliers = 0;
			inlierDistance[0] = 0;
			Ellipse recEllipse;
			double maxWeight = 0.0;
			double theta = atan(distance(center, radCenter) / RichConst::aglTransmissionHeight());
			double phi = atan((centerB - radB) / (centerA - radA));
			if (centerA - radA < 0) phi += TMath::Pi();
			double maxBinContent = weightModelMap -> GetBinContent(weightModelMap -> GetMaximumBin());
			
			// Find the most possible Cherenkov angle theta_i
			vector < double > vecInlierDistace;
			for (double theta_i = 0.1; theta_i < 0.4; theta_i += 0.005) {
				
				// Compute parameters of ellipse
				Ellipse ellipse_i = GetEllipse(radCenter, theta, phi, theta_i);
				vector < double > ellipseCenter {ellipse_i.centerX, ellipse_i.centerY};
				
				// Compute the likelihood of theta_i
				double weight = 1.0;
				double hitsAvailable = hits.size();
				vector < double > tempInlierDistance;
				for (int i = 0; i < hits.size(); i ++) {
					
					double hit2Center = distance(ellipseCenter, hits[i]);
					if (hit2Center > 22.0) {
						hitsAvailable --;
						continue;
					}
					
					// rExpected being the expected distance of the point on the ellipse to the center, which is on the direction of photon hit to center.
					double rExpected = GetEllipseRadius(ellipse_i, hits[i]);
#ifdef WEIGHT_MODEL_CORRECTION
					weight *= exp(weightModelMap -> GetBinContent(weightModelMap -> FindBin(rExpected - hit2Center)) * 1.0 / maxBinContent);
#else
					weight *= weightModel(hit2Center, rExpected);
#endif
					if (fabs(rExpected - hit2Center) < 3) {
						tempInlierDistance.push_back(rExpected - hit2Center);
					}
				}
				
				if (tempInlierDistance.size() < 4) continue;
				if (hitsAvailable * 1.0 / hits.size() < 0.5 || hitsAvailable < 4) continue;
				
				if (maxWeight < weight) {
					maxWeight = weight;
					recTheta = theta_i;
					recEllipse = ellipse_i;
					vecInlierDistace = tempInlierDistance;
				}
			}
			
			if (recTheta == 0) continue;
			
			inliers = (int) vecInlierDistace.size();
			for (int ii = 0; ii < inliers; ii ++) {
				inlierDistance[ii] = vecInlierDistace[ii];
			}
			
#ifdef USING_N_CORRECTION
			double n2 = refractiveMap -> GetBinContent(refractiveMap -> FindBin(radA, radB));
#else
			double n2 = 1.055;
#endif
			
			// do refraction correction
#ifdef USING_REFRACTION_CORRECTION
			tree.recBeta = 1.0 / n2 / sqrt(1.0 - square(1.0 / n2 * sin(recTheta)));
#else
			tree.recBeta = 1.0 / n2 / cos(recTheta);
#endif
//			if (tree.recBeta > 1.0) {
//				tree.recBeta = max(1.0 / n2 / cos(recTheta - 0.005), 1.0);
//			}
			tree.recMass = tree.trRigidity * tree.trInnerCharge / tree.recBeta * sqrt(1.0 - square(tree.recBeta));

			// Correction for position changes from light guide to PMT pixel?
			
			outTree -> Fill();
			
#ifdef DRAW_ELLIPSE
			DrawEllipse(iEvent, recEllipse, hits);
#endif
		}
	}

	outFile -> cd();
	outTree -> Write();
	outFile -> Save();
	outFile -> Close();
	
	return 0;
}
