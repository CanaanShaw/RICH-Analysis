//
//  testFunctions.h
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/6/8.
//

#ifndef testFunctions_h
#define testFunctions_h

#include "../tooLib/sharedHeader.h"
vector < double > FindCrossPoint(vector < double > hits, vector < double > ringCenter,
								 vector < double > radCenter, double theta, double phi, double z = 47.17) {
	// Intersection point of each line between photon hit and emission point, with spin plane going through ringCenter.
	
	double t = -1.0 * (hits[0] * sin(theta) * cos(phi) + hits[1] * sin(theta) * sin(phi)) / ((radCenter[0] - hits[0]) * sin(theta) * cos(phi) + (radCenter[1] - hits[1]) * sin(theta) * sin(phi) + z * cos(theta));
	
	double xOut = hits[0] + (radCenter[0] - hits[0]) * t;
	double yOut = hits[1] + (radCenter[1] - hits[1]) * t;
	double zOut = z * t;
	
	return vector < double > {xOut, yOut, zOut};
}

#endif /* testFunctions_h */
