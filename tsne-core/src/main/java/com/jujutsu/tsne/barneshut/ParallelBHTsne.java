package com.jujutsu.tsne.barneshut;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.RecursiveAction;

public class ParallelBHTsne extends BHTSne {
	
	private ForkJoinPool gradientPool;
	private ExecutorService	gradientCalculationPool;
	
	class RecursiveGradientCalculator extends RecursiveAction {
		final static long serialVersionUID = 1L;
		int startRow = -1;
		int endRow = -1;
		int limit = 100;
		SPTree tree;
		double[][] neg_f;
		double theta;
		AtomicDouble sum_Q;

		public RecursiveGradientCalculator(SPTree tree, double [][] neg_f , double theta, 
				AtomicDouble sum_Q, int startRow, int endRow, int ll) {
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
	
	class ParallelGradientCalculator implements Callable<Double> {
		final static long serialVersionUID = 1L;
		int row = -1;
		int limit = 100;
		ParallelSPTree tree;
		double[][] neg_f;
		double theta;
		
		public ParallelGradientCalculator(ParallelSPTree tree, double [][] neg_f , double theta, int row, int ll) {
			this.limit = ll;
			this.row = row;
			this.tree = tree;
			this.neg_f = neg_f; 
			this.theta = theta; 
		}
		@Override
		public Double call() {
			return tree.computeNonEdgeForces(row, theta, neg_f[row], 0.0);
		}   
	}
		
	@Override
	double[] run(double[] X, int N, int D, int no_dims, double perplexity, int max_iter, double theta) {
		gradientPool = new ForkJoinPool(Runtime.getRuntime().availableProcessors());
		double [] Y = super.run(X, N, D, no_dims, perplexity, max_iter, theta);
		gradientPool.shutdown();
		return Y;
	}

	// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
	@Override
	void computeGradient(double [] P, int [] inp_row_P, int [] inp_col_P, double [] inp_val_P, 
			double [] Y, int N, int D, double [] dC, double theta)
	{
		// Construct space-partitioning tree on current map
		ParallelSPTree tree = new ParallelSPTree(D, Y, N);

		// Compute all terms required for t-SNE gradient
		double [] pos_f = new double[N * D];
		double [][] neg_f = new double[N][D];

		tree.computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, N, pos_f);

		//AtomicDouble Asum_Q = new AtomicDouble();
		double sum_Q = 0;
		//RecursiveGradientCalculator dslr = new RecursiveGradientCalculator(tree, neg_f, theta, sum_Q, 0, N, 20);                
		//gradientPool.invoke(dslr);
		
		List<ParallelGradientCalculator> calculators = new ArrayList<>();
		for(int n = 0; n < N; n++) {
			calculators.add(new ParallelGradientCalculator(tree, neg_f, theta, n, 20));
		}
		List<Future<Double>> results;
		try {
			results = gradientCalculationPool.invokeAll(calculators);
			for (Future<Double> result : results) {
				double tmp = result.get();
				//System.out.println("Got Q: " + tmp);
				sum_Q += tmp;
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (ExecutionException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		//for(int n = 0; n < N; n++) tree.computeNonEdgeForces(n, theta, neg_f[n], sum_Q);
		//System.out.println("Final Q: " + sum_Q);
		// Compute final t-SNE gradient
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < D; j++) {
				dC[i*D+j] = pos_f[i*D+j] - (neg_f[i][j] / sum_Q);
			}
		}
	}

}
