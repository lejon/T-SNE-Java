package com.jujutsu.tsne.barneshut;

import java.util.Comparator;

public class DistanceComparator implements Comparator<DataPoint> {
	DataPoint refItem; 
	DistanceComparator(DataPoint refItem) {
		this.refItem = refItem;
	}

	@Override
	public int compare(DataPoint o1, DataPoint o2) {
		return o1.euclidean_distance(refItem) < o2.euclidean_distance(refItem) ? -1 :
			(o1.euclidean_distance(refItem) > o2.euclidean_distance(refItem) ? 1 : 0);
	}
}