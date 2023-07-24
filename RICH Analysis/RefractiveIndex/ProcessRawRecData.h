//
//  GetIndexMap.h
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/6/16.
//

#ifndef ProcessRawRecData_h
#define ProcessRawRecData_h

#include "../../tooLib/tooLib.h"
#include "../../tooLib/sharedHeader.h"
#include "../../tooLib/myTree.h"

void ProcessRawRecData(const char * readPath, const char * writePath) {
	// Generate refractive index map and weight model based on raw reconstruction data
	
	cout << "Reading from file: " << readPath << endl;
	myTree tree(readPath, "treeRec");
	tree.SetAMSAddress();
	tree.SetTrackerAddress();
	tree.SetRecAddress();
	int inliers;
	double inlierDistance[MAXIMUM];
	tree.pTree -> SetBranchAddress("inliers", &inliers);
	tree.pTree -> SetBranchAddress("inlierDistance", inlierDistance);
	
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
	
	TProfile2D * RefractiveMap = new TProfile2D("RefractiveMap", "RefractiveMap", 220, -63.75, 63.75, 220, -63.75, 63.75);
	TH1D * WeightModel = new TH1D("WeightModel", "WeightModel", 1500, -3.0, 3.0);
	long long entries = tree.GetEntries();
	for (int i = 0; i < entries; i++) {

		tree.GetEntry(i);
		if (!tree.Select()) continue;
		if (tree.recTheta == 0) continue;
		if (tree.trRigidity < 40) continue;
		if (tree.trInnerCharge < 0.8 || tree.trInnerCharge > 1.4) continue;
		if (tree.trTheta > 0.1) continue;
		if (tree.recBeta < 0.5) continue;
		
		double dx = -0.08;
		double dy = 0.075;
		
#ifdef USING_DYNAMIC_ALIGNMENT_CORRECTION
			
		alignmentTime = 0;
		int iAlignmentEvent = 0;
		while (alignmentTime < tree.AMSTime) {
			alignmentTree -> GetEntry(iAlignmentEvent ++);
		}
		
		dx = biasX;
		dy = biasY;
#endif
		
#ifdef USING_N_CORRECTION
		double n0 = refractiveMap -> GetBinContent(refractiveMap -> FindBin(tree.radX + dx, tree.radY + dy));
		if (n == 0) n = 1.055;
#else
		double n0 = 1.055;
#endif
		
		double n = tree.recBeta * n0;
		RefractiveMap -> Fill(tree.trRadX + dx, tree.trRadY + dy, n);
		
		for (int iInlier = 0; iInlier < inliers; iInlier ++) {
			WeightModel -> Fill(inlierDistance[iInlier]);
		}
	}

	TString mapPath = writePath;
	mapPath += "RefractiveMapItered.root";
	TFile * fOut = new TFile(mapPath, "RECREATE");
	fOut -> cd();
	RefractiveMap -> Write();
	fOut -> Save();
	fOut -> Close();
	
	TString modelPath = writePath;
	modelPath += "WeightModel.root";
	TFile * fOut2 = new TFile(modelPath, "RECREATE");
	fOut2 -> cd();
	WeightModel -> Write();
	fOut2 -> Save();
	fOut2 -> Close();
	delete RefractiveMap;
	delete WeightModel;
}

#endif /* ProcessRawRecData_h */
