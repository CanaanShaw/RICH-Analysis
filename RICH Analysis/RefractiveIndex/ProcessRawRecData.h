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
	tree.SetTrackerAddress();
	tree.SetRecAddress();
	int inliers;
	double inlierDistance[MAXIMUM];
	tree.pTree -> SetBranchAddress("inliers", &inliers);
	tree.pTree -> SetBranchAddress("inlierDistance", inlierDistance);
	
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
		
		double n = tree.recBeta * 1.055;
		RefractiveMap -> Fill(tree.trRadX + 0.09, tree.trRadY - 0.075, n);
		
		for (int iInlier = 0; iInlier < inliers; iInlier ++) {
			WeightModel -> Fill(inlierDistance[iInlier]);
		}
	}

	TString mapPath = writePath;
	mapPath += "RefractiveMap.root";
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
