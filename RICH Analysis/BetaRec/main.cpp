//
//  main.cpp
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/4/14.
//

#include "../../tooLib/sharedHeader.h"
#include "../../tooLib/tooLib.h"
#include "../../tooLib/myTree.h"

//#include "weightModel.h"
#include "recEllipse.h"
#include "defs.h"
#include "testFunctions.h"


int main(int argc, const char * argv[]) {
	
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - // File Initialization
	
#ifdef LOCAL_TEST
	
	string dataPath = "/Users/canaanshaw/Desktop/CppFiles/RICHAnalysis.data/DataSample/out.1479168000.root";
	myTree tree(dataPath.c_str(), "dailyData");
	
	string outPath	= "/Users/canaanshaw/Desktop/CppFiles/RICHAnalysis.data/BetaRec/rec.1479168000.DirectThetaRec.root";
	TFile * outFile = new TFile(outPath.c_str(), "RECREATE");
	
	string nMapPath = "/Users/canaanshaw/Desktop/CppFiles/RICHAnalysis.data/BetaRec/RefractiveMap.root";
//	string nMapPath = "/Users/canaanshaw/Desktop/CppFiles/RefractiveMap.root";
	cout << "Loading refractive index file: " << nMapPath << endl;
//	TFile * mapFile = new TFile("/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/RefractiveMap.root", "READ");
	TFile * mapFile = new TFile(nMapPath.c_str(), "READ");
	TH2D * refractiveMap = (TH2D *) mapFile -> Get("RefractiveMap");
	
//	string weightModelPath = "/Users/canaanshaw/Desktop/CppFiles/WeightModel_Default.root";
//	cout << "Loading weight model file: " << weightModelPath << endl;
////	TFile * modelFile = new TFile("/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/weightModel", "READ");
//	TFile * modelFile = new TFile(weightModelPath.c_str(), "READ");
//	TH1D * weightModelMap = (TH1D *) modelFile -> Get("WeightModel");

	#ifdef USING_DYNAMIC_ALIGNMENT_CORRECTION
		string alignmentPath = "/Users/canaanshaw/Desktop/CppFiles/2023-7-1 HoughTransform/alignmentResult.root";
		cout << "Loading dynamic alignment file: " << alignmentPath << endl;
		TFile * alignmentFile = new TFile(alignmentPath.c_str(), "READ");
		TTree * alignmentTree = (TTree *) alignmentFile -> Get("alignment");
	
		int alignmentTime;
		alignmentTree -> SetBranchAddress("time", &alignmentTime);
		
		double biasX;
		double biasY;
		alignmentTree -> SetBranchAddress("biasX", &biasX);
		alignmentTree -> SetBranchAddress("biasY", &biasY);
	#endif
#else
	
	cout << "Reading from file: " << argv[1] << endl;
	myTree tree(argv[1], "dailyData");
	
	cout << "Writing to file: " << argv[2] << endl;
	TFile * outFile = new TFile(argv[2], "RECREATE");
	
	string nMapPath = "/afs/cern.ch/user/j/jianan/private/RICHAnalysis/RefractiveIndex/analysis/RefractiveMap.root";
	cout << "Loading refractive index file: " << nMapPath << endl;
	TFile * mapFile = new TFile(nMapPath.c_str(), "READ");
	TH2D * refractiveMap = (TH2D *) mapFile -> Get("RefractiveMap");
	
//	string weightModelPath = "/afs/cern.ch/user/j/jianan/private/RICHAnalysis/RefractiveIndex/analysis/WeightModel.root";
//	cout << "Loading weight model file: " << weightModelPath << endl;
//	TFile * modelFile = new TFile(weightModelPath.c_str(), "READ");
//	TH1D * weightModelMap = (TH1D *) modelFile -> Get("WeightModel");
	
	#ifdef USING_DYNAMIC_ALIGNMENT_CORRECTION
		string alignmentPath = "/afs/cern.ch/user/j/jianan/private/RICHAnalysis/Alignment/analysis/alignmentResult.root";
		cout << "Loading dynamic alignment file: " << alignmentPath << endl;
		TFile * alignmentFile = new TFile(alignmentPath.c_str(), "READ");
		TTree * alignmentTree = (TTree *) alignmentFile -> Get("alignment");

		int alignmentTime;
		alignmentTree -> SetBranchAddress("time", &alignmentTime);
		
		double biasX;
		double biasY;
		alignmentTree -> SetBranchAddress("biasX", &biasX);
		alignmentTree -> SetBranchAddress("biasY", &biasY);
	#endif
#endif
	
	tree.SetAMSAddress();
	tree.SetRichAddress();
	tree.SetTrackerAddress();
	tree.SetTofAddress();

	outFile -> cd();
	TTree * outTree = new TTree("treeRec", "Rec Result");
	tree.MakeAMSAddress(outTree);
	tree.MakeRecAddress(outTree);
	tree.MakeRichAddress(outTree);
	tree.MakeTrackerAddress(outTree);

	double emissionPointX;
	double emissionPointY;
	outTree -> Branch("emissionPointX", &emissionPointX, "emissionPointX/D");
	outTree -> Branch("emissionPointY", &emissionPointY, "emissionPointY/D");
	
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - // Standard Correction and Hits Selection

	// Static correction of RICH spacial position
	double dx = -0.07;
	double dy = 0.07;
	
	// Start processing data
	for (int iEvent = 0; iEvent < tree.GetEntries(); iEvent ++) {
		
		tree.GetEntry(iEvent);
		if (iEvent % 10000 == 0)
			cout << "Current: " << iEvent << "	out of " << tree.GetEntries() << endl;
		if (iEvent > 150000) break;

		if (tree.Select()) {
			
#ifdef USING_DYNAMIC_ALIGNMENT_CORRECTION
			
			alignmentTime = 0;
			int iAlignmentEvent = 0;
			while (alignmentTime < tree.AMSTime) {
				alignmentTree -> GetEntry(iAlignmentEvent ++);
			}
			
			dx = biasX;
			dy = biasY;
#endif
			
			// Tracker extrapolated photon hit position with correction
			// dx = recCenterX - tree.trPointX
			double centerA = tree.trPointX + dx;
			double centerB = tree.trPointY + dy;
			vector < double > center {centerA, centerB};
			
			// Tracker extrapolated radiation top plane point with correction
			double radA = tree.trRadX + dx;
			double radB = tree.trRadY + dy;
			vector < double > radCenter {radA, radB};
			
			// Assume that the cherenkov radiation is on average emitted at the vertical center of AGL.
			emissionPointX = radCenter[0] + RichConst::aglHeight / 2.0 * tan(tree.trTheta) * cos(tree.trPhi + TMath::Pi());
			emissionPointY = radCenter[1] + RichConst::aglHeight / 2.0 * tan(tree.trTheta) * sin(tree.trPhi + TMath::Pi());
			vector < double > emissionPoint {emissionPointX, emissionPointY};

			// Čerenkov radiation cone center extrapolated from top of the radiator plane to PMT plane.
			// This is more accurate than trPointX and trPointY.
			double radExtCenterX = tree.trRadX + cos(tree.trPhi + TMath::Pi()) * RichConst::richHeight() * sin(tree.trTheta);
			double radExtCenterY = tree.trRadY + sin(tree.trPhi + TMath::Pi()) * RichConst::richHeight() * sin(tree.trTheta);
			vector < double > radExtCenter {radExtCenterX, radExtCenterY};
			
			// Photon hits position without ones crossed by particles
			vector < vector < double > > hits;
			
			// hitsRaw is used to draw the original PMT signals.
			auto hitsRaw = hits;
			
			for (int iHit = 0; iHit < tree.nHits; iHit ++) {
				
				vector < double > tempHit {
					tree.pointX[iHit],
					tree.pointY[iHit],
					distance(center, vector < double > {tree.pointX[iHit], tree.pointY[iHit]})
				};
				
				hitsRaw.push_back(tempHit);
				
//				// in case the result is affected by reflected hits, we're trying to abandon the hits on the edge
//				if (distance(vector < double > {0, 0}, tempHit) > 65.0) continue;
//
//
//
				if (!tree.isCrossed[iHit] && tempHit[2] <= 22.0) {
					hits.push_back(tempHit);
				}
			}
			
			// Scan possible radius range
			int maxInliers = 0;
			double rInitial = 0;
			vector < vector < double > > hitsSelected;
			for (double r = 4.0; r < 20.0; r += 1.0) {
				int tempInliers = 0;
				vector < vector < double > > tempHit;
				for (int i = 0; i < hits.size(); i ++) {
					if (hits[i][2] > r - 2 && hits[i][2] < r + 2) {
						tempInliers ++;
						tempHit.push_back(hits[i]);
					}
				}

				if (tempInliers >= maxInliers) {
					maxInliers = tempInliers;
					rInitial = r;
					hitsSelected = tempHit;
				}
			}
			
//			if (hitsSelected.size() < 5) continue;
			if (hitsSelected.size() < 5 || hitsSelected.size() * 1.0 / hits.size() < 0.6) continue;
			
#ifdef USING_REFRACTION_CORRECTION
			
			for (int i = 0; i < hitsSelected.size(); i ++) {
				
				#ifdef USING_N_CORRECTION
					
					double nTemp = refractiveMap -> GetBinContent(refractiveMap -> FindBin(emissionPoint[0], emissionPoint[1]));
//					double nTemp = refractiveMap -> GetBinContent(refractiveMap -> FindBin(radA, radB));
					if (nTemp == 0) nTemp = RichConst::refractiveIndexDefault;
															  
					// Actually this corresponding refractive index n should be the value at the bottom of the radiator instead of rad center.
					hitsSelected[i] = RefractionCorrection(emissionPoint, hitsSelected[i], nTemp);
				#else
					hitsSelected[i] = RefractionCorrection(emissionPoint, hitsSelected[i]);
				#endif
			}
#endif
			
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - // Beta Reconstruction
			
			tree.recTheta = 0.0;
			
//			TCanvas * c = new TCanvas("c", "c", 800, 800);
			TH1D * thetaHist = new TH1D("thetaHist", "thetaHist", 100, 0, 0.5);
			for (vector < double > iHit : hitsSelected) {
				double a = sqrt(square(distance(emissionPoint, iHit)) + square(RichConst::aglTransmissionHeight()));
				double b = sqrt(square(distance(emissionPoint, radExtCenter)) + square(RichConst::aglTransmissionHeight()));
				double c = distance(iHit, radExtCenter);
				double theta_i = acos((a * a + b * b - c * c) / 2.0 / a / b);
				thetaHist -> Fill(theta_i);
				
			}
			
			TF1 * f = new TF1("f", "gaus", 0.0, 0.5);
			thetaHist -> Fit(f, "WW E M Q", "", 0.0, 0.5);
//			string path = to_string(iEvent);
//			string path2 = "/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/thetaHist/" + path + ".png";
//			thetaHist -> Draw();
//			c -> Print(path2.c_str());
//			delete c;
			
			tree.recTheta = f -> GetParameter(1);
			delete thetaHist;
			delete f;
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - // Fill Tree
			
#ifdef USING_N_CORRECTION
//			double n = refractiveMap -> GetBinContent(refractiveMap -> FindBin(radA, radB));
			double n = refractiveMap -> GetBinContent(refractiveMap -> FindBin(emissionPoint[0], emissionPoint[1]));
			if (n == 0) continue;
#else
			double n = RichConst::refractiveIndexDefault;
#endif
			
			tree.recBeta = 1.0 / n / cos(tree.recTheta);

#ifdef USING_BETA_BOUNDARY_CORRECTION
			if (tree.recBeta > 1.0) {
				tree.recBeta = 1.0 / n / sqrt(1.0 - square(1.0 / n * sin(tree.recTheta - 0.0025)));
			}
#endif
			tree.recMass = tree.trRigidity * tree.trInnerCharge / tree.recBeta * sqrt(1.0 - square(tree.recBeta));

			// Correction for position changes from light guide to PMT pixel?
			
			outTree -> Fill();
			
#ifdef DRAW_ELLIPSE
			if (tree.trTheta > 0.2)
			DrawEllipse(iEvent, recEllipse, hitsRaw, hitsSelected, center, radExtCenter, radCenter);
#endif
		}
	}

	outFile -> cd();
	outTree -> Write();
	outFile -> Save();
	outFile -> Close();
	
	return 0;
}
