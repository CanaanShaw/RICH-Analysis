//
//  functions.h
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/6/7.
//

#ifndef weightModel_h
#define weightModel_h

double computeWeight (double x, std::vector < double > r, double band1 = 1.45, double band2 = 2.25, double ratio = 0.9) {

	double weight = 1.0;
	unsigned long size = r.size();
	for (int i = 0; i < size; i ++) {

		double dist = fabs(x - r[i]);
		if (dist < 0.3) {
			weight *= exp(1.0);
			continue;
		}
		weight *= exp(ratio / (1.0 + fabs(pow(dist / band1, 4))) + (1.0 - ratio) / (1.0 + pow(dist / band2, 3)));
	}

	weight /= size;
	return weight;
}

#endif /* weightModel_h */