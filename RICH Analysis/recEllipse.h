//
//  recEllipse.h
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/6/12.
//

#ifndef recEllipse_h
#define recEllipse_h

#include "tooLib.h"

struct Ellipse{
	// five parameters to discribe an ellipse
	double centerX;
	double centerY;
	double semiMajorAxis;
	double semiMinorAxis;
	double phi;
};

Ellipse GetEllipse(vector < double > radCenter, double trTheta, double trPhi, double theta_i){
	// Return these five parameters of a ellipse:
	// 		center position X, center position Y, semi-major axis, semi-minor axis, spin angle wrt X axis
	
	double semiMajorAxis = 0.5 * RichConst::aglTransmissionHeight() * (tan(trTheta + theta_i) - tan(trTheta - theta_i));
	double semiMinorAxis = RichConst::aglTransmissionHeight() / cos(trTheta) * tan(theta_i);

	// Distance from radCenter to ellipseCenter from top view.
	double rad2Center = RichConst::aglTransmissionHeight() * tan(trTheta - theta_i) + semiMajorAxis;
	
	double ellipseCenterX = radCenter[0] + rad2Center * cos(trPhi);
	double ellipseCenterY = radCenter[1] + rad2Center * sin(trPhi);
//	if ((ellipseCenter[0] - centerA) / (centerA - radA) < 0) cout << "Wrong\n";
//	if ((ellipseCenter[1] - centerB) / (centerB - radB) < 0) cout << "Wrong\n";
	return Ellipse {ellipseCenterX, ellipseCenterY, semiMajorAxis, semiMinorAxis, trPhi};
}

double GetEllipseRadius(const Ellipse & ellipse, vector < double > hit) {
	
	// Angle from hit to ellipse center wrt X axis
	double theta = atan((ellipse.centerY - hit[1]) / (ellipse.centerX - hit[0]));
	
	// Angle from hit to ellipse center wrt major axis
	double phi = theta - ellipse.phi;
	return sqrt(square(ellipse.semiMajorAxis * cos(phi)) + square(ellipse.semiMinorAxis * sin(phi)));
}

void DrawEllipse(const int & iEvent, const Ellipse & ellipse, const vector < vector < double > > & hits) {
	
	TCanvas * c = new TCanvas("c", "Ellipse", 800, 800);
	TGraph * g = new TGraph();
	for (vector < double > iHit : hits) {
		
		g -> AddPoint(iHit[0], iHit[1]);
	}
	
	gPad -> Modified();
	g -> GetXaxis() -> SetLimits(-67.0, 67.0);
	g -> SetMinimum(-67.0);
	g -> SetMaximum(67.0);
	gPad -> Update();

	g -> SetMarkerStyle(8);
	g -> SetMarkerSize(0.3);

	c -> cd();
	g -> Draw("AP");
	
	TEllipse * circle = new TEllipse(ellipse.centerX, ellipse.centerY, ellipse.semiMajorAxis, ellipse.semiMinorAxis, 0, 360, 180 * ellipse.phi / TMath::Pi());
	circle -> SetFillColor(0);
	circle -> SetFillStyle(0);
	circle -> SetLineColor(2);
	circle -> DrawClone("SAME");
	
	TEllipse * circle2 = new TEllipse(0, 0, 67.0);
	circle2 -> SetFillColor(0);
	circle2 -> SetFillStyle(0);
	circle2 -> SetLineColorAlpha(kBlack, 0.1);
	circle2 -> SetLineWidth(3);
	circle2 -> DrawClone("SAME");
	
	TBox * box = new TBox(-32.0, -32.0, 32.0, 32.0);
	box -> SetFillColor(0);
	box -> SetFillStyle(0);
	box -> SetLineWidth(3);
	box -> SetLineColorAlpha(kBlack, 0.3);
	box -> Draw();

	string path = to_string(iEvent);
	string path2 = "/Users/canaanshaw/Desktop/CppFiles/betaAnalysis/recEllipse/" + path + ".png";
	c -> Print(path2.c_str());
	
	delete c;
	delete g;
	delete circle;
	delete box;
	delete circle2;
}

#endif /* recEllipse_h */
