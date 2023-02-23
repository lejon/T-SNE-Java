package com.jujutsu.tsne.barneshut;

import static java.lang.Math.sqrt;

public class EuclideanDistance implements Distance{

	public EuclideanDistance() {
	}

	@Override
	public double distance(DataPoint d1, DataPoint d2) {
	    double dd = .0;
	    double diff;
	    for(int d = 0; d < d1._D; d++) {
	        diff = (d1._x[d] - d2._x[d]);
	        dd += diff * diff;
	    }
	    return sqrt(dd);
	}

}
