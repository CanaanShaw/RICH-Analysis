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

#endif /* recEllipse_h */
