package com.jujutsu.tsne.barneshut;

import static java.lang.Math.sqrt;

public class EuclideanDistance implements Distance{

	public EuclideanDistance() {
	}

	@Override
	public double distance(DataPoint d1, DataPoint d2) {
	    double dd = .0;
	    double [] x1 = d1._x;
	    double [] x2 = d2._x;
	    double diff;
	    for(int d = 0; d < d1._D; d++) {
	        diff = (x1[d] - x2[d]);
	        dd += diff * diff;
	    }
	    return sqrt(dd);

	}

}
