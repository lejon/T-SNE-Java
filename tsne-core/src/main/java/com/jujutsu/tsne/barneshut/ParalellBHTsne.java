package com.jujutsu.tsne.barneshut;

import java.util.concurrent.RecursiveAction;
import java.util.concurrent.atomic.AtomicLong;

public class ParalellBHTsne extends BarnesHutTSne {
	
	class RecursiveGradientCalculator extends RecursiveAction {
		final static long serialVersionUID = 1L;
		int startRow = -1;
		int endRow = -1;
		int limit = 100;
		SPTree tree;
		double[][] neg_f;
		double theta;
		AtomicLong sum_Q;

		public RecursiveGradientCalculator(SPTree tree, double [][] neg_f , double theta, 
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
				invokeAll(new RecursiveGradientCalculator(tree,neg_f, theta, sum_Q, startDoc1, endDoc1, limit),
						new RecursiveGradientCalculator(tree,neg_f, theta, sum_Q, startDoc2, endDoc2, limit));
			}
		}
	}

	// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
	@Override
	void computeGradient(double [] P, int [] inp_row_P, int [] inp_col_P, double [] inp_val_P, 
			double [] Y, int N, int D, double [] dC, double theta)
	{
		// Construct space-partitioning tree on current map
		ParalellSPTree tree = new ParalellSPTree(D, Y, N);

		// Compute all terms required for t-SNE gradient
		AtomicLong sum_Q = new AtomicLong();
		double [] pos_f = new double[N * D];
		double [][] neg_f = new double[N][D];

		tree.computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, N, pos_f);
		for(int n = 0; n < N; n++) tree.computeNonEdgeForces(n, theta, neg_f[n], sum_Q);

		// Compute final t-SNE gradient
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < D; j++) {
				dC[i*D+j] = pos_f[i*D+j] - (neg_f[i][j] / Double.longBitsToDouble(sum_Q.get()));
			}
		}
	}

}
