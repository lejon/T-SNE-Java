package com.jujutsu.tsne.barneshut;

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
		return dist.distance(o1, refItem) < dist.distance(o2, refItem) ? -1 :
			(dist.distance(o1, refItem) > dist.distance(o2, refItem) ? 1 : 0);
	}
}