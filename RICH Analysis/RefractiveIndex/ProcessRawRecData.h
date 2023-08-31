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
#include "GausProfile.h"
#include "defs.h"

void ProcessRawRecData(const char * readPath, const char * writePath) {
	// Generate refractive index map and weight model based on raw reconstruction data
	
#ifdef USING_N_CORRECTION
	string nMapPath = "/afs/cern.ch/user/j/jianan/private/RICHAnalysis/RefractiveIndex/analysis/RefractiveMap.root";
	cout << "Loading refractive index file: " << nMapPath << endl;
	TFile * mapFile = new TFile(nMapPath.c_str(), "READ");
	TH2D * refractiveMap = (TH2D *) mapFile -> Get("RefractiveMap");
#endif
	
	cout << "Reading from file: " << readPath << endl;
	myTree tree(readPath, "treeRec");
	tree.SetAMSAddress();
	tree.SetTrackerAddress();
	tree.SetRecAddress();
	tree.SetRichAddress();
	
	tree.pTree -> SetBranchAddress("emissionPointX", &emissionPointX);
	tree.pTree -> SetBranchAddress("emissionPointY", &emissionPointY);
	
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
	
	
	RadiatorProfile radProfile(220);
	long long entries = tree.GetEntries();
	for (int i = 0; i < entries; i++) {

		tree.GetEntry(i);
		if (!tree.Select()) continue;
		if (tree.recTheta == 0) continue;
		if (tree.trRigidity < 50) continue;
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
//		double n0 = refractiveMap -> GetBinContent(refractiveMap -> FindBin(tree.trRadX + dx, tree.trRadY + dy));
		double n0 = refractiveMap -> GetBinContent(refractiveMap -> FindBin(emissionPointX, emissionPointY));
		if (n0 == 0) n0 = RichConst::refractiveIndexDefault;
#else
		double n0 = RichConst::refractiveIndexDefault;
#endif
		
		double n = tree.recBeta * n0;
//		radProfile.Fill(tree.trRadX + dx, tree.trRadY + dy, n);
		radProfile.Fill(emissionPointX, emissionPointY, n);
	}
	
	radProfile.GausFit();
	TProfile2D * refMap = radProfile.GetResult();
	TProfile2D * refMapErr = radProfile.GetResultErr();

	TString mapPath = writePath;
	mapPath += "/RefractiveMap.root";
	TFile * fOut = new TFile(mapPath, "RECREATE");
	fOut -> cd();
	refMap -> Write();
	refMapErr -> Write();
	radProfile.Write(fOut);
	fOut -> Save();
	fOut -> Close();
}

#endif /* ProcessRawRecData_h */
