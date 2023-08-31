//
//  testFunctions.h
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/6/8.
//

#ifndef testFunctions_h
#define testFunctions_h

#include "../../tooLib/sharedHeader.h"
vector < double > FindCrossPoint(vector < double > hits, vector < double > ringCenter,
								 vector < double > radCenter, double theta, double phi, double z = 47.17) {
	// Intersection point of each line between photon hit and emission point, with spin plane going through ringCenter.
	
	double t = -1.0 * (hits[0] * sin(theta) * cos(phi) + hits[1] * sin(theta) * sin(phi)) / ((radCenter[0] - hits[0]) * sin(theta) * cos(phi) + (radCenter[1] - hits[1]) * sin(theta) * sin(phi) + z * cos(theta));
	
	double xOut = hits[0] + (radCenter[0] - hits[0]) * t;
	double yOut = hits[1] + (radCenter[1] - hits[1]) * t;
	double zOut = z * t;
	
	return vector < double > {xOut, yOut, zOut};
}

vector < double > RefractionCorrection (vector < double > radCenter,
										vector < double > hit, double n = RichConst::refractiveIndexDefault,
										double radiatorThickness = RichConst::aglHeight / 2.0,
										//	   ^^^^^^^^^^^^^^^^^ should this term be 2.5 / 2.0 ?
										double z = RichConst::richHeight() - RichConst::aglHeight
										) {
	// Return the actual position the photon hits on the PMT plane without the refraction effect.
	double d = distance(radCenter, hit);
	const double eps = 1e-2;
	
	
	// Change forward propagation to backward propagation
	// in order to minimum the error
	double x1 = 0;
	double x2 = 0;
	double theta = 0;
	for (x1 = 0; x1 <= d / 2.0; x1 += 1e-5) {
		theta = atan(x1 / radiatorThickness);
		x2 = z * tan(asin(n * sin(theta)));
		if (fabs(x1 + x2 - d) < eps) break;
	}
	
	if (x1 > d / 2.0) {
		return hit;
	}
	
	double dActual = x1 * (z + radiatorThickness) / radiatorThickness;
	
	double phi = atan((hit[1] - radCenter[1]) / (hit[0] - radCenter[0]));
	if (hit[0] - radCenter[0] < 0) phi += TMath::Pi();
	
//	cout << radCenter[0] + dActual * cos(phi) << ":	" << hit[0] << endl;
	return vector < double > {radCenter[0] + dActual * cos(phi), radCenter[1] + dActual * sin(phi), hit[2]};
}



#endif /* testFunctions_h */
