package com.jujutsu.tsne.barneshut;

import java.util.List;

public class TreeSearchResult {
	int n;
	List<Double> distances;
	List<DataPoint> indices;

	public TreeSearchResult(List<DataPoint> indices, List<Double> distances, int n) {
		this.indices = indices;
		this.distances = distances;
		this.n = n;
	}

	public List<DataPoint> getIndices() {
		return indices;
	}

	public List<Double> getDistances() {
		return distances;
	}


	public int getIndex() {
		return n;
	}

}