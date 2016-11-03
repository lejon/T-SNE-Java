package com.jujutsu.tsne.barneshut;

import static java.lang.Math.sqrt;

public class ParallelSPTree extends SPTree {

	public ParallelSPTree(int D, double[] inp_data, int N) {
		super(D, inp_data, N);
	}

	public ParallelSPTree(int D, double[] inp_data, int N, double[] inp_corner, double[] inp_width) {
		super(D, inp_data, N, inp_corner, inp_width);
	}

	public ParallelSPTree(int D, double[] inp_data, double[] inp_corner, double[] inp_width) {
		super(D, inp_data, inp_corner, inp_width);
	}

	public ParallelSPTree(SPTree inp_parent, int D, double[] inp_data, double[] inp_corner, double[] inp_width) {
		super(inp_parent, D, inp_data, inp_corner, inp_width);
	}

	public ParallelSPTree(SPTree inp_parent, int D, double[] inp_data, int N, double[] inp_corner, double[] inp_width) {
		super(inp_parent, D, inp_data, N, inp_corner, inp_width);
	}
	
	@Override
	SPTree[] getTreeArray(int no_children) {
		return new ParallelSPTree[no_children];
	}
	
	@Override
	SPTree getNewTree(SPTree root, double[] new_corner, double[] new_width) {
		return new ParallelSPTree(root, dimension, data, new_corner, new_width);
	}

	// Compute non-edge forces using Barnes-Hut algorithm
	@Override
	void computeNonEdgeForces(int point_index, double theta, double [] neg_f, Object accumulator)
	{
		AtomicDouble sum_Q = (AtomicDouble) accumulator;
		// Make sure that we spend no time on empty nodes or self-interactions
		if(cum_size == 0 || (is_leaf && size == 1 && index[0] == point_index)) return;

		// Compute distance between point and center-of-mass
		double D = .0;
		int ind = point_index * dimension;
		for(int d = 0; d < dimension; d++) buff[d] = data[ind + d] - center_of_mass[d];
		for(int d = 0; d < dimension; d++) D += buff[d] * buff[d];

		// Check whether we can use this node as a "summary"
		double max_width = 0.0;
		double cur_width;
		for(int d = 0; d < dimension; d++) {
			cur_width = boundary.getWidth(d);
			max_width = (max_width > cur_width) ? max_width : cur_width;
		}
		if(is_leaf || max_width / sqrt(D) < theta) {
			// Compute and add t-SNE force between point and current node
			D = 1.0 / (1.0 + D);
			double mult = cum_size * D;
			sum_Q.addAndGet(mult);
			mult *= D;
			for(int d = 0; d < dimension; d++) neg_f[d] += mult * buff[d];
		}
		else {
			// Recursively apply Barnes-Hut to children
			for(int i = 0; i < no_children; i++) children[i].computeNonEdgeForces(point_index, theta, neg_f, sum_Q);
		}
	}

}
