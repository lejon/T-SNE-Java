package com.jujutsu.tsne.barneshut;

import java.util.Arrays;
import java.util.Comparator;

public class DistanceComparator implements Comparator<DataPoint> {
	DataPoint refItem; 
	Distance dist;
	
	DistanceComparator(DataPoint refItem) {
		this.refItem = refItem;
		this.dist = new EuclideanDistance();
	}
	
	DistanceComparator(DataPoint refItem, Distance dist) {
		this.refItem = refItem;
		this.dist = dist;
	}

	@Override
	public int compare(DataPoint o1, DataPoint o2) {
		//System.out.println("Reference " + refItem._ind + " " + Arrays.toString(refItem._x) + " ===> Comparing " + o1._ind + " " + Arrays.toString(o1._x) + " < " + o2._ind + " " + Arrays.toString(o2._x) + "### Distance: x=" + dist.distance(o1, refItem) + " <==> " + dist.distance(o2, refItem));
		return dist.distance(o1, refItem) < dist.distance(o2, refItem) ? -1 :
			(dist.distance(o1, refItem) > dist.distance(o2, refItem) ? 1 : 0);
	}
}