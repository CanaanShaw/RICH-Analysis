//
//  recPack.h
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/6/12.
//

#ifndef recPack_h
#define recPack_h

#include "tooLib.h"
#include "sharedHeader.h"

double computeCubeWeight(double variable1, double variable2, double bandWidth = 3.0) {

	return 1.0 / (1.0 + fabs(cube((variable1 - variable2) / bandWidth)));
}

double computeFourthWeight(double variable1, double variable2, double bandWidth = 2.0) {

	return 1.0 / (1.0 + square(square((variable1 - variable2) / bandWidth)));
}

vector < double > recWeight(vector < vector < double > > hitsClean, vector < double > center, vector < double > emissionPoint) {

	// hitsClean are PMT signals excluding isCrossed ones.
	// returns: vector < double > {recR, weight, massCenterQuality, radiusQuality}.

	bool isNaF = (fabs(emissionPoint[0]) < 17.4 && fabs(emissionPoint[1]) < 17.4);

	double recR = 0;

	// Quality indicator for scatter degree of rings.
	// Higher value (with maximum of 1) with hits uniformly distribute on the ring.
	double massCenterX = 0;
	double massCenterY = 0;
	double massCenterQuality = 0;

	// Quality indicator for variation between primary recR and final recR.
	double radiusQuality = 0;

	int nHits = hitsClean.size();

	// Actual number of hits used for reconstruction.
	int hitsUsed = nHits;

	// Sloppy guess of Čerenkov ring radii range.
	const double aglMinR = 5.0;
	const double aglMaxR = 20.0;
	const double NaFMinR = 15.0;
	const double NaFMaxR = 85.0;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Select hits to do primary rec

	double meanR = 0;

	for (int iHit = 0; iHit < nHits; iHit ++) {

		double r_i = distance(center, vector < double > {hitsClean[iHit][0], hitsClean[iHit][1]});

		// Do boxing filter
		if (!isNaF) {

			if (r_i < aglMinR || r_i > aglMaxR) {

				hitsUsed --;
				continue;
			}
		} else {

			if (r_i < NaFMinR || r_i > NaFMaxR) {

				hitsUsed --;
				continue;
			}
		}

		meanR += r_i;
		massCenterX += hitsClean[iHit][0];
		massCenterY += hitsClean[iHit][1];
	}

	// In case number of hits being too low or noise being too strong.
	// This term have wasted too many entries!!
	if (hitsUsed < 3 || hitsUsed * 1.0 / nHits < 0.5) return vector < double > {-1.0, -1.0, -1.0, -1.0};

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Do rec with weight

	meanR /= hitsUsed;
	double weight = 0.3;

	// Initial guess of refined radius.
	double refinedR = meanR;
	for (int iHit = 0; iHit < nHits; iHit ++) {

		double r_i = distance(center, vector < double > {hitsClean[iHit][0], hitsClean[iHit][1]});
		if (!isNaF) {

			if (r_i < aglMinR || r_i > aglMaxR || r_i > meanR * 1.1 + 2 || r_i < meanR * 0.9 - 2) {

				continue;
			}
		} else {

			if (r_i < NaFMinR || r_i > NaFMaxR || r_i > meanR * 1.1 + 2 || r_i < meanR * 0.9 - 2) {

				continue;
			}
		}

		// This weight model is mainly affected by exponential term and denominator term around (refinedR - r_i).
		double w_i = computeCubeWeight(refinedR, r_i);
		if (fabs(refinedR - r_i) < 0.2) w_i = 1.0;
		// double w_i = 1.0 / (1.0 + square((refinedR - r_i) / (1.0 / 6.0 * refinedR + 0.5)));

		// Update refinedR and weight.
		refinedR = (weight * refinedR + w_i * r_i) / (weight + w_i);
		weight += w_i;
	}
	
	recR = refinedR;

	// Normalize weight.
	weight /= (hitsUsed + 1.0);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Estimate rec quality

	massCenterX /= hitsUsed;
	massCenterY /= hitsUsed;
	massCenterQuality = 1.0 - distance(center, vector < double > {massCenterX, massCenterY}) / recR;

	radiusQuality = 1.0 / (1.0 + 0.5 * fabs(recR - meanR));
	return vector < double > {recR, weight, massCenterQuality, radiusQuality};
}

vector < double > recWeightWithEdgeCorrection(vector < vector < double > > hitsClean, vector < double > center, vector < double > emissionPoint) {

	// hitsClean are PMT signals excluding isCrossed ones.
	// returns: vector < double > {recR, weight, massCenterQuality, radiusQuality}.

	bool isNaF = (fabs(emissionPoint[0]) < 17.4 && fabs(emissionPoint[1]) < 17.4);

	double recR = 0;

	// Quality indicator for scatter degree of rings.
	// Higher value (with maximum of 1) with hits uniformly distribute on the ring.
	double massCenterX = 0;
	double massCenterY = 0;
	double massCenterQuality = 0;

	// Quality indicator for variation between primary recR and final recR.
	double radiusQuality = 0;

	int nHits = hitsClean.size();

	// Actual number of hits used for reconstruction.
	int hitsUsed = nHits;

	// Sloppy guess of Čerenkov ring radii range.
	const double aglMinR = 3.0;
	const double aglMaxR = 20.0;
	const double NaFMinR = 15.0;
	const double NaFMaxR = 85.0;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Select hits to do primary rec

	double meanR = 0;

	for (int iHit = 0; iHit < nHits; iHit ++) {

		double r_i = distance(center, vector < double > {hitsClean[iHit][0], hitsClean[iHit][1]});
		double r_0 = distance(vector < double > {0.0, 0.0}, vector < double > {hitsClean[iHit][0], hitsClean[iHit][1]});
		if (r_0 > 52 && r_i < 10) {
			r_i = 67.0 - r_0 + 67.0 - distance(vector < double > {0.0, 0.0}, center);
		}

		// Do boxing filter
		if (!isNaF) {

			if (r_i < aglMinR || r_i > aglMaxR) {

				hitsUsed --;
				continue;
			}
		} else {

			if (r_i < NaFMinR || r_i > NaFMaxR) {

				hitsUsed --;
				continue;
			}
		}

		meanR += r_i;
		massCenterX += hitsClean[iHit][0];
		massCenterY += hitsClean[iHit][1];
	}

	// In case number of hits being too low or noise being too strong.
	// This term have wasted too many entries!!
	if (hitsUsed < 3 || hitsUsed * 1.0 / nHits < 0.5) return vector < double > {-1.0, -1.0, -1.0, -1.0};

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Do rec with weight

	meanR /= hitsUsed;
	double weight = 0.3;

	// Initial guess of refined radius.
	double refinedR = meanR;
	for (int iHit = 0; iHit < nHits; iHit ++) {

		double r_i = distance(center, vector < double > {hitsClean[iHit][0], hitsClean[iHit][1]});
		double r_0 = distance(vector < double > {0.0, 0.0}, vector < double > {hitsClean[iHit][0], hitsClean[iHit][1]});
		if (r_0 > 52) {
			r_i = 67.0 - r_0 + 67.0 - distance(vector < double > {0.0, 0.0}, center);
		}

		if (!isNaF) {

			if (r_i < aglMinR || r_i > aglMaxR || r_i > meanR * 1.1 + 2 || r_i < meanR * 0.9 - 2) {

				continue;
			}
		} else {

			if (r_i < NaFMinR || r_i > NaFMaxR || r_i > meanR * 1.1 + 2 || r_i < meanR * 0.9 - 2) {

				continue;
			}
		}

		// This weight model is mainly affected by exponential term and denominator term around (refinedR - r_i).
		double w_i = computeCubeWeight(refinedR, r_i);
		// double w_i = 1.0 / (1.0 + square((refinedR - r_i) / (1.0 / 6.0 * refinedR + 0.5)));

		// Update refinedR and weight.
		refinedR = (weight * refinedR + w_i * r_i) / (weight + w_i);
		weight += w_i;
	}
	
	recR = refinedR;

	// Normalize weight.
	weight /= (hitsUsed + 1.0);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Estimate rec quality

	massCenterX /= hitsUsed;
	massCenterY /= hitsUsed;
	massCenterQuality = 1.0 - distance(center, vector < double > {massCenterX, massCenterY}) / recR;

	radiusQuality = 1.0 / (1.0 + 0.5 * fabs(recR - meanR));
	return vector < double > {recR, weight, massCenterQuality, radiusQuality};
}

vector < double > recWeightWithEdgeCorrection2(vector < vector < double > > hitsClean, vector < double > center, vector < double > emissionPoint) {

	// hitsClean are PMT signals excluding isCrossed ones.
	// returns: vector < double > {recR, weight, massCenterQuality, radiusQuality}.

	bool isNaF = (fabs(emissionPoint[0]) < 17.4 && fabs(emissionPoint[1]) < 17.4);

	double centerDistance = distance(center, vector < double > {0.0, 0.0});
	bool isReflected = (centerDistance > 52.0 && centerDistance < 67.0);

	double recR = 0;

	// Quality indicator for scatter degree of rings.
	// Higher value (with maximum of 1) with hits uniformly distribute on the ring.
	double massCenterX = 0;
	double massCenterY = 0;
	double massCenterQuality = 0;

	// Quality indicator for variation between primary recR and final recR.
	double radiusQuality = 0;

	int nHits = hitsClean.size();

	// Actual number of hits used for reconstruction.
	int hitsUsed = nHits;

	// Sloppy guess of Čerenkov ring radii range.
	const double aglMinR = 3.0;
	const double aglMaxR = 20.0;
	const double NaFMinR = 15.0;
	const double NaFMaxR = 85.0;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Select hits to do primary rec

	double meanR = 0;

	for (int iHit = 0; iHit < nHits; iHit ++) {

		double r_i = distance(center, vector < double > {hitsClean[iHit][0], hitsClean[iHit][1]});

		// Do boxing filter
		if (!isNaF) {

			if (r_i < aglMinR || r_i > aglMaxR) {

				hitsUsed --;
				continue;
			}
		} else {

			if (r_i < NaFMinR || r_i > NaFMaxR) {

				hitsUsed --;
				continue;
			}
		}

		meanR += r_i;
		massCenterX += hitsClean[iHit][0];
		massCenterY += hitsClean[iHit][1];
	}

	// In case number of hits being too low or noise being too strong.
	// This term have wasted too many entries!!
	if (hitsUsed < 3 || hitsUsed * 1.0 / nHits < 0.5) return vector < double > {-1.0, -1.0, -1.0, -1.0};

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Do rec with weight

	meanR /= hitsUsed;
	double weight = 0.3;

	// Initial guess of refined radius.
	double refinedR = meanR;
	for (int iHit = 0; iHit < nHits; iHit ++) {

		double r_i = distance(center, vector < double > {hitsClean[iHit][0], hitsClean[iHit][1]});

		if (!isNaF) {

			if (r_i < aglMinR || r_i > aglMaxR || r_i > meanR * 1.1 + 2 || r_i < meanR * 0.9 - 2) {

				continue;
			}
		} else {

			if (r_i < NaFMinR || r_i > NaFMaxR || r_i > meanR * 1.1 + 2 || r_i < meanR * 0.9 - 2) {

				continue;
			}
		}

		// This weight model is mainly affected by exponential term and denominator term around (refinedR - r_i).
		double w_i = computeCubeWeight(refinedR, r_i);
		// double w_i = 1.0 / (1.0 + square((refinedR - r_i) / (1.0 / 6.0 * refinedR + 0.5)));

		// Update refinedR and weight.
		refinedR = (weight * refinedR + w_i * r_i) / (weight + w_i);
		weight += w_i;
	}
	
	recR = refinedR;

	// Normalize weight.
	weight /= (hitsUsed + 1.0);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Estimate rec quality

	massCenterX /= hitsUsed;
	massCenterY /= hitsUsed;
	massCenterQuality = 1.0 - distance(center, vector < double > {massCenterX, massCenterY}) / recR;

	radiusQuality = 1.0 / (1.0 + 0.5 * fabs(recR - meanR));

	if (isReflected) {

		double theta = atan(center[1] / center[0]);
		vector < double > reflectedCenter = vector < double > {2 * 67.0 * cos(theta) - center[0], 2 * 67.0 * sin(theta) - center[1]};
		vector < double > reflectedRing = recWeightWithEdgeCorrection2(hitsClean, reflectedCenter, emissionPoint);
		if (reflectedRing[0] > recR && reflectedRing[1] >= weight) {

			cout << "corrected\n";
			return reflectedRing;
		}
	}

	return vector < double > {recR, weight, massCenterQuality, radiusQuality};
}

vector < double > recFastRANSAC(vector < vector < double > > hitsClean, vector < double > center) {

	// returns: vector < double > {ringCenterX, ringCenterY, recR}

	// Sloppy guess of Čerenkov ring radii range.
	const double aglMinR = 4.0;
	const double aglMaxR = 18.0;
	const double NaFMinR = 15.0;
	const double NaFMaxR = 85.0;

	double massCenterX = 0;
	double massCenterY = 0;
	double massCenterQuality = 0;
	double radiusQuality = 0;

	int nHits = hitsClean.size();
	int hitsUsed = nHits;

	double recR = 0;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Select hits to do primary rec

	double meanR = 0;

	for (int iHit = 0; iHit < nHits; iHit ++) {

		vector < double > hit_i = {hitsClean[iHit][0], hitsClean[iHit][1]};
		double r_i = distance(center, hit_i);
		if (r_i < aglMinR || r_i > aglMaxR) {

			hitsUsed --;
			continue;
		}

		meanR += r_i;
		massCenterX += hitsClean[iHit][0];
		massCenterY += hitsClean[iHit][1];
	}

	if (hitsUsed < 3 || hitsUsed * 1.0 / nHits < 0.5) return vector < double > {-1.0, -1.0, -1.0, -1.0};

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - do RANSAC rec

	meanR /= hitsUsed;
	double maxWeight = 0;

	for (int i = 0; i < nHits; i ++) {

		vector < double > hit_i = {hitsClean[i][0], hitsClean[i][1]};
		double r_i = distance(center, hit_i);

		if (r_i < aglMinR || r_i > aglMaxR) {

			continue;
		}

		// Calculate inliers.
		double weightSum = 0;
		for (int iHit = 0; iHit < nHits; iHit ++) {

			vector < double > hit = {hitsClean[iHit][0], hitsClean[iHit][1]};
			double dist = distance(center, hit);
			if (fabs(dist - r_i) > 1.5) continue;
			double weight_i = computeFourthWeight(dist, r_i, 1.0);

			weightSum += weight_i;
		}

		// Update recRing.
		if (maxWeight < weightSum) {

			recR = r_i;
			maxWeight = weightSum;
		}
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Estimate rec quality

	// maxWeight /=hitsUsed;
	// massCenterX /= hitsUsed;
	// massCenterY /= hitsUsed;
	// massCenterQuality = 1.0 - distance(center, vector < double > {massCenterX, massCenterY}) / recR;
	// radiusQuality = 1.0 / (1.0 + 0.5 * fabs(recR - meanR));

	// return vector < double > {recR, maxWeight, massCenterQuality, radiusQuality};

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Post Monte-Carlo Based Refinement

	std::random_device rd;
	std::default_random_engine engine(rd());
	std::normal_distribution < double > gaussian(0.0, recR / 50);

	int maxTries = 50;
	double refinedR = recR;
	for (int i = 0; i < maxTries; i ++) {

		double rTemp = recR + gaussian(engine);
		double weightSum = 0;
		for (int iHit = 0; iHit < nHits; iHit ++) {

			vector < double > hit = {hitsClean[iHit][0], hitsClean[iHit][1]};
			double dist = distance(center, hit);
			if (fabs(dist - rTemp) > 1.5) continue;
			double weight_i = computeFourthWeight(dist, rTemp, 1.0);

			weightSum += weight_i;
		}

		if (maxWeight < weightSum) {

			// cout << "good\n";
			refinedR = rTemp;
			maxWeight = weightSum;
		}
	}

	recR = refinedR;
	maxWeight /=hitsUsed;
	massCenterQuality = 1.0 - distance(center, vector < double > {massCenterX, massCenterY}) / recR;
	radiusQuality = 1.0 / (1.0 + 0.5 * fabs(recR - meanR));
	
	return vector < double > {recR, maxWeight, massCenterQuality, radiusQuality};
}

vector < double > recFastRANSACWithCorrection(vector < vector < double > > hitsClean, vector < double > center) {

	// returns: vector < double > {ringCenterX, ringCenterY, recR}

	// Sloppy guess of Čerenkov ring radii range.
	const double aglMinR = 3.0;
	const double aglMaxR = 20.0;
	const double NaFMinR = 15.0;
	const double NaFMaxR = 85.0;

	double massCenterX = 0;
	double massCenterY = 0;
	double massCenterQuality = 0;
	double radiusQuality = 0;

	int nHits = hitsClean.size();
	int hitsUsed = nHits;

	double recR = 0;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Select hits to do primary rec

	double meanR = 0;

	for (int iHit = 0; iHit < nHits; iHit ++) {

		vector < double > hit_i = {hitsClean[iHit][0], hitsClean[iHit][1]};
		double r_i = distance(center, hit_i);
		if (r_i < aglMinR || r_i > aglMaxR) {

			hitsUsed --;
			continue;
		}

		meanR += r_i;
		massCenterX += hitsClean[iHit][0];
		massCenterY += hitsClean[iHit][1];
	}

	if (hitsUsed < 3 || hitsUsed * 1.0 / nHits < 0.5) return vector < double > {-1.0, -1.0, -1.0, -1.0};

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - do RANSAC rec

	meanR /= hitsUsed;
	double maxWeight = 0;

	for (int i = 0; i < nHits; i ++) {

		vector < double > hit_i = {hitsClean[i][0], hitsClean[i][1]};
		double r_i = distance(center, hit_i);

		if (r_i < aglMinR || r_i > aglMaxR) {

			continue;
		}

		// Calculate inliers.
		double weightSum = 0;
		for (int iHit = 0; iHit < nHits; iHit ++) {

			vector < double > hit = {hitsClean[iHit][0], hitsClean[iHit][1]};
			double dist = distance(center, hit);
			double weight_i = computeCubeWeight(dist, r_i);

			weightSum += weight_i;
		}

		// Update recRing.
		if (maxWeight < weightSum) {

			recR = r_i;
			maxWeight = weightSum;
		}
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Estimate rec quality

	// maxWeight /=hitsUsed;
	// massCenterX /= hitsUsed;
	// massCenterY /= hitsUsed;
	// massCenterQuality = 1.0 - distance(center, vector < double > {massCenterX, massCenterY}) / recR;
	// radiusQuality = 1.0 / (1.0 + 0.5 * fabs(recR - meanR));

	// return vector < double > {recR, maxWeight, massCenterQuality, radiusQuality};

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Post Monte-Carlo Based Refinement

	std::random_device rd;
	std::default_random_engine engine(rd());
	std::normal_distribution < double > gaussian(0.0, recR / 10);

	int maxTries = 50;
	double refinedR = recR;
	for (int i = 0; i < maxTries; i ++) {

		double rTemp = recR + gaussian(engine);
		double weightSum = 0;
		for (int iHit = 0; iHit < nHits; iHit ++) {

			vector < double > hit = {hitsClean[iHit][0], hitsClean[iHit][1]};
			double dist = distance(center, hit);
			double weight_i = computeFourthWeight(dist, rTemp, 2.0);

			weightSum += weight_i;
		}

		if (maxWeight < weightSum) {

			// cout << "good\n";
			refinedR = rTemp;
			maxWeight = weightSum;
		}
	}

	recR = refinedR;
	maxWeight /=hitsUsed;
	massCenterQuality = 1.0 - distance(center, vector < double > {massCenterX, massCenterY}) / recR;
	radiusQuality = 1.0 / (1.0 + 0.5 * fabs(recR - meanR));

	double centerDistance = distance(center, vector < double > {0.0, 0.0});
	bool isReflected = (centerDistance > 52.0 && centerDistance < 67.0);

	if (isReflected) {

		double theta = atan(center[1] / center[0]);
		vector < double > reflectedCenter = vector < double > {2 * 67.0 * cos(theta) - center[0], 2 * 67.0 * sin(theta) - center[1]};
		vector < double > reflectedRing = recFastRANSACWithCorrection(hitsClean, reflectedCenter);
		if (reflectedRing[0] > recR && reflectedRing[1] >= maxWeight) {

			cout << "corrected\n";
			return reflectedRing;
		}
	}
	
	return vector < double > {recR, maxWeight, massCenterQuality, radiusQuality};
}

vector < double > recIndependentRANSAC(vector < vector < double > > hitsClean) {

	// returns: vector < double > {ringCenterX, ringCenterY, recR}

	std::random_device rd;
	std::default_random_engine engine(rd());
	std::uniform_real_distribution < float > rand(0, 1);

	int nHits = hitsClean.size();

	int maxInliers = 0;
	int maxTries = 300;

	const double threshold = 1.0;

	vector < double > recRing {0, 0, 0};

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - do RANSAC rec

	for (int i = 0; i < maxTries; i ++) {

		// Randomly select 3 hits to fit a ring
		int idx1 = max < int > (round(nHits * rand(engine) - 1.0), 0.0);
		int idx2 = max < int > (round(nHits * rand(engine) - 1.0), 0.0);
		int idx3 = max < int > (round(nHits * rand(engine) - 1.0), 0.0);

		if (idx1 == idx2 || idx1 == idx3 || idx2 == idx3) {

			i --;
			continue;
		}

		vector < vector < double > > points = {hitsClean[idx1], hitsClean[idx2], hitsClean[idx3]};
		vector < double > recRingTemp = recCircle(points);

		// Calculate inliers.
		double inliers = 0;
		for (int iHit = 0; iHit < nHits; iHit ++) {

			vector < double > center = {recRingTemp[0], recRingTemp[1]};
			vector < double > hit 	 = {hitsClean[iHit][0], hitsClean[iHit][1]};
			double dist = distance(center, hit) - recRingTemp[2];

			if (fabs(dist) < threshold) {
				inliers ++;
			}
		}

		// Update recRing.
		if (maxInliers < inliers) {

			recRing = recRingTemp;
			maxInliers = inliers;
		}
	}

	return recRing;
}


vector < double > recTest(vector < vector < double > > hitsClean, vector < double > center, vector < double > emissionPoint) {

	// hitsClean are PMT signals excluding isCrossed ones.
	// returns: vector < double > {recR, weight, massCenterQuality, radiusQuality}.

	bool isNaF = (fabs(emissionPoint[0]) < 17.4 && fabs(emissionPoint[1]) < 17.4);

	double recR = 0;

	// Quality indicator for scatter degree of rings.
	// Higher value (with maximum of 1) with hits uniformly distribute on the ring.
	double massCenterX = 0;
	double massCenterY = 0;
	double massCenterQuality = 0;

	// Quality indicator for variation between primary recR and final recR.
	double radiusQuality = 0;

	int nHits = hitsClean.size();

	// Actual number of hits used for reconstruction.
	int hitsUsed = nHits;

	// Sloppy guess of Čerenkov ring radii range.
	const double aglMinR = 5.0;
	const double aglMaxR = 20.0;
	const double NaFMinR = 15.0;
	const double NaFMaxR = 85.0;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Select hits to do primary rec

	double meanR = 0;

	for (int iHit = 0; iHit < nHits; iHit ++) {

		double r_i = distance(center, vector < double > {hitsClean[iHit][0], hitsClean[iHit][1]});

		// Do boxing filter
		if (!isNaF) {

			if (r_i < aglMinR || r_i > aglMaxR) {

				hitsUsed --;
				continue;
			}
		} else {

			if (r_i < NaFMinR || r_i > NaFMaxR) {

				hitsUsed --;
				continue;
			}
		}

		meanR += r_i;
		massCenterX += hitsClean[iHit][0];
		massCenterY += hitsClean[iHit][1];
	}

	// In case number of hits being too low or noise being too strong.
	// This term have wasted too many entries!!
	if (hitsUsed < 3 || hitsUsed * 1.0 / nHits < 0.5) return vector < double > {-1.0, -1.0, -1.0, -1.0};

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Do rec with weight

	double maxWeight = 0;
	double rFinal = -1.0;
	for (double dr = 12.0; dr < 18.0; dr += 0.1) {

		double weight_i = 1;
		for (int iHit = 0; iHit < nHits; iHit ++) {

			double r_i = distance(center, vector < double > {hitsClean[iHit][0], hitsClean[iHit][1]});
			if (!isNaF) {

				if (r_i < aglMinR || r_i > aglMaxR) {

					continue;
				}
			} else {

				if (r_i < NaFMinR || r_i > NaFMaxR) {

					continue;
				}
			}

			// This weight model is mainly affected by exponential term and denominator term around (refinedR - r_i).
			double dw = computeCubeWeight(dr, r_i, 1);
			weight_i *= max < double > (dw, 0.1);
			// cout << weight_i;
			// if (fabs(refinedR - r_i) < 0.2) w_i = 1.0;
			// double w_i = 1.0 / (1.0 + square((refinedR - r_i) / (1.0 / 6.0 * refinedR + 0.5)));

			// Update refinedR and weight.
			// refinedR = (weight * refinedR + w_i * r_i) / (weight + w_i);
			// weight += w_i;
		}

		if (maxWeight < weight_i) {
			maxWeight = weight_i;
			rFinal = dr;
		}
	}

	recR = rFinal;
	// cout << recR << endl;
	// Normalize weight.
	// weight /= (hitsUsed + 1.0);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - Estimate rec quality

	// massCenterX /= hitsUsed;
	// massCenterY /= hitsUsed;
	// massCenterQuality = 1.0 - distance(center, vector < double > {massCenterX, massCenterY}) / recR;

	// radiusQuality = 1.0 / (1.0 + 0.5 * fabs(recR - meanR));
	return vector < double > {recR, maxWeight, massCenterQuality, radiusQuality};
}




















#endif /* recPack_h */
