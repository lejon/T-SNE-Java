package com.jujutsu.tsne.barneshut;

import static java.lang.Math.min;
import static java.lang.Math.sqrt;

public class DataPoint {
	
	int _ind;
	double [] _x;
	int _D;
	
	public DataPoint() {
        _D = 1;
        _ind = -1;
    }

	public DataPoint(int D, int ind, double [] x) {
		_D = D;
		_ind = ind;
		_x = x.clone();
	}
	
	@Override
	public String toString() {
		String xStr = "";
		for (int i = 0; i < min(20,_x.length); i++) {
			xStr += _x[i] + ", ";
		}
		return "DataPoint (index=" + _ind+ ", Dim=" + _D + ", point=" + xStr + ")"; 
	}

	public int index() { return _ind; }
	int dimensionality() { return _D; }
	double x(int d) { return _x[d]; }
	
	public double euclidean_distance( DataPoint t1 ) {
	    double dd = .0;
	    double [] x1 = t1._x;
	    double [] x2 = _x;
	    double diff;
	    for(int d = 0; d < t1._D; d++) {
	        diff = (x1[d] - x2[d]);
	        dd += diff * diff;
	    }
	    return sqrt(dd);
	}
	
	public static double euclidean_distance( DataPoint t1, DataPoint t2 ) {
	    double dd = .0;
	    double [] x1 = t1._x;
	    double [] x2 = t2._x;
	    double diff;
	    for(int d = 0; d < t1._D; d++) {
	        diff = (x1[d] - x2[d]);
	        dd += diff * diff;
	    }
	    return sqrt(dd);
	}
}
