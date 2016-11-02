package com.jujutsu.tsne.barneshut;

import static java.lang.Math.sqrt;

import java.util.concurrent.RecursiveAction;
import java.util.concurrent.atomic.AtomicLong;

import com.jujutsu.tsne.barneshut.ParalellBHTsne.RecursiveGradientCalculator;

public class ParalellSPTree extends SPTree {

	public ParalellSPTree(int D, double[] inp_data, int N) {
		super(D, inp_data, N);
	}

	public ParalellSPTree(int D, double[] inp_data, int N, double[] inp_corner, double[] inp_width) {
		super(D, inp_data, N, inp_corner, inp_width);
	}

	public ParalellSPTree(int D, double[] inp_data, double[] inp_corner, double[] inp_width) {
		super(D, inp_data, inp_corner, inp_width);
	}

	public ParalellSPTree(SPTree inp_parent, int D, double[] inp_data, double[] inp_corner, double[] inp_width) {
		super(inp_parent, D, inp_data, inp_corner, inp_width);
	}

	public ParalellSPTree(SPTree inp_parent, int D, double[] inp_data, int N, double[] inp_corner, double[] inp_width) {
		super(inp_parent, D, inp_data, N, inp_corner, inp_width);
	}
	
	class RecursiveForceCalculator extends RecursiveAction {
		final static long serialVersionUID = 1L;
		int startRow = -1;
		int endRow = -1;
		int limit = 100;
		SPTree tree;
		double[][] neg_f;
		double theta;
		AtomicLong sum_Q;

		public RecursiveForceCalculator(SPTree tree, double [][] neg_f , double theta, 
				AtomicLong sum_Q, int startRow, int endRow, int ll) {
			this.limit = ll;
			this.startRow = startRow;
			this.endRow = endRow;
			this.tree = tree;
			this.neg_f = neg_f; 
			this.theta = theta; 
			this.sum_Q = sum_Q;
		}

		@Override
		protected void compute() {
			if ( (endRow-startRow) <= limit ) {
				for (int row = startRow; row < endRow; row++) {
					tree.computeNonEdgeForces(row, theta, neg_f[row], sum_Q);
				}
			}
			else {
				int range = (endRow-startRow);
				int startDoc1 = startRow;
				int endDoc1 = startRow + (range / 2);
				int startDoc2 = endDoc1;
				int endDoc2 = endRow;
				invokeAll(new RecursiveForceCalculator(tree,neg_f, theta, sum_Q, startDoc1, endDoc1, limit),
						new RecursiveForceCalculator(tree,neg_f, theta, sum_Q, startDoc2, endDoc2, limit));
			}
		}
	}

	// Compute non-edge forces using Barnes-Hut algorithm
	void computeNonEdgeForces(int point_index, double theta, double [] neg_f, Object accumulator)
	{
		AtomicLong sum_Q = (AtomicLong) accumulator;
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
			sum_Q.addAndGet(Double.doubleToLongBits(mult));
			mult *= D;
			for(int d = 0; d < dimension; d++) neg_f[d] += mult * buff[d];
		}
		else {
			// Recursively apply Barnes-Hut to children
			for(int i = 0; i < no_children; i++) children[i].computeNonEdgeForces(point_index, theta, neg_f, sum_Q);
		}
	}

}
