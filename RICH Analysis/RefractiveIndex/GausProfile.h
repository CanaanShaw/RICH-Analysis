//
//  GausProfile.h
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/7/26.
//

#ifndef GausProfile_h
#define GausProfile_h

#include "../../tooLib/sharedHeader.h"
#include "../../tooLib/tooLib.h"

struct RadiatorTile {
	double lowerBounds[3];
	double upperBounds[3];
	double boxCenterX;
	double boxCenterY;
	double boxCenterZ;
	
	TH1D * tileHist;
	double n;
	double nErr;
};

class RadiatorProfile {
	
private:
	int nBins;
	const double xMin = -63.75;
	const double xMax =  63.75;
	const double yMin = -63.75;
	const double yMax =  63.75;
	vector < vector < RadiatorTile > > tileArr;
	
public:
	
	RadiatorProfile () {};
	RadiatorProfile (int nBins);
	~RadiatorProfile ();
	vector < int > FindBin (double x, double y);
	void Fill (double x, double y, double value);
	void GausFit ();
	TProfile2D * GetResult ();
	TProfile2D * GetResultErr ();
	void Write (TFile * outFile);
};


RadiatorProfile::RadiatorProfile (int segment) {
	this -> nBins = segment;
	const double sideLengthX = (xMax - xMin) / nBins;
	const double sideLengthY = (yMax - yMin) / nBins;
	
	this -> tileArr.resize(nBins);
	for (int i = 0; i < nBins; i ++) {
		this -> tileArr[i].resize(nBins);
	}
	
	for (int x_i = 0; x_i < nBins; x_i ++) {
		for (int y_i = 0; y_i < nBins; y_i ++) {
			
			this -> tileArr[x_i][y_i].lowerBounds[0] = xMin + x_i * sideLengthX;
			this -> tileArr[x_i][y_i].lowerBounds[1] = yMin + y_i * sideLengthY;
			this -> tileArr[x_i][y_i].lowerBounds[2] = 0;
			
			this -> tileArr[x_i][y_i].upperBounds[0] = xMin + (x_i + 1) * sideLengthX;
			this -> tileArr[x_i][y_i].upperBounds[1] = yMin + (y_i + 1) * sideLengthY;
			this -> tileArr[x_i][y_i].upperBounds[2] = RichConst::aglHeight;
			
			this -> tileArr[x_i][y_i].boxCenterX = xMin + (x_i + 0.5) * sideLengthX;
			this -> tileArr[x_i][y_i].boxCenterY = yMin + (y_i + 0.5) * sideLengthY;
			this -> tileArr[x_i][y_i].boxCenterZ = 0.5 * RichConst::aglHeight;
			
			string nameX = convert < string > (x_i);
			string nameY = convert < string > (y_i);
			string name = nameX + ": " + nameY;
			this -> tileArr[x_i][y_i].tileHist = new TH1D(name.c_str(), name.c_str(), 90, 1.0, 1.1);
			this -> tileArr[x_i][y_i].n = 0;
			this -> tileArr[x_i][y_i].nErr = 0;
		}
	}
}

RadiatorProfile::~RadiatorProfile() {
	for (int x_i = 0; x_i < nBins; x_i ++) {
		for (int y_i = 0; y_i < nBins; y_i ++) {
			delete this -> tileArr[x_i][y_i].tileHist;
		}
	}
}

vector < int > RadiatorProfile::FindBin (double x, double y) {
	
	int x_i = 0;
	int y_i = 0;
	bool flag = 0;
	for (x_i = 0; x_i < this -> nBins && !flag; x_i ++) {
		for (y_i = 0; y_i < this -> nBins && !flag; y_i ++) {
			bool xBox = (x >= this -> tileArr[x_i][y_i].lowerBounds[0] && x < this -> tileArr[x_i][y_i].upperBounds[0]);
			bool yBox = (y >= this -> tileArr[x_i][y_i].lowerBounds[1] && y < this -> tileArr[x_i][y_i].upperBounds[1]);
			if (xBox && yBox) flag = 1;
		}
	}
	
	if (x_i == this -> nBins && y_i == this -> nBins) return vector < int > {0, 0};
	return vector < int > {x_i, y_i};
}

void RadiatorProfile::Fill (double x, double y, double value) {
	vector < int > index = this -> FindBin(x, y);
	if (index == vector < int > {0, 0}) cout << "RadiatorProfile::FindBin() - - - - Error\n";
	
	int xIndex = index[0];
	int yIndex = index[1];
	this -> tileArr[xIndex][yIndex].tileHist -> Fill(value);
	
	return;
}

void RadiatorProfile::GausFit () {
	for (int x_i = 0; x_i < this -> nBins; x_i ++) {
		for (int y_i = 0; y_i < this -> nBins; y_i ++) {
			if (this -> tileArr[x_i][y_i].tileHist -> GetEntries() != 0) {
				TF1 * f = new TF1("f", "gaus", 1.0, 1.1);
				this -> tileArr[x_i][y_i].tileHist -> Fit(f, "WWEM", "", 1.0, 1.1);
				this -> tileArr[x_i][y_i].n = f -> GetParameter(1);
				this -> tileArr[x_i][y_i].nErr = f -> GetParError(1);
				delete f;
				
				if (tileArr[x_i][y_i].n < 1.03 || tileArr[x_i][y_i].n > 1.06) {
					tileArr[x_i][y_i].n = 0;
					tileArr[x_i][y_i].nErr = 0;
				}
			}
		}
	}
}

TProfile2D * RadiatorProfile::GetResult () {
	TProfile2D * RefractiveMap = new TProfile2D("RefractiveMap",
												"RefractiveMap",
												this -> nBins,
												xMin,
												xMax,
												this -> nBins,
												yMin,
												yMax);
	for (int x_i = 0; x_i < this -> nBins; x_i ++) {
		for (int y_i = 0; y_i < this -> nBins; y_i ++) {
			RefractiveMap -> Fill(tileArr[x_i][y_i].boxCenterX, tileArr[x_i][y_i].boxCenterY, tileArr[x_i][y_i].n);
		}
	}
	
	return RefractiveMap;
}

TProfile2D * RadiatorProfile::GetResultErr () {
	TProfile2D * RefractiveMap = new TProfile2D("RefractiveMapErr",
												"RefractiveMapErr",
												this -> nBins,
												xMin,
												xMax,
												this -> nBins,
												yMin,
												yMax);
	for (int x_i = 0; x_i < this -> nBins; x_i ++) {
		for (int y_i = 0; y_i < this -> nBins; y_i ++) {
			RefractiveMap -> Fill(tileArr[x_i][y_i].boxCenterX, tileArr[x_i][y_i].boxCenterY, tileArr[x_i][y_i].nErr);
		}
	}
	
	return RefractiveMap;
}

void RadiatorProfile::Write (TFile * outFile) {
	for (int x_i = 0; x_i < this -> nBins; x_i ++) {
		for (int y_i = 0; y_i < this -> nBins; y_i ++) {
			outFile -> cd();
			tileArr[x_i][y_i].tileHist -> Write();
		}
	}
	
	outFile -> Save();
}

#endif /* GausProfile_h */
