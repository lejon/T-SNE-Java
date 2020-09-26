package com.jujutsu.tsne.barneshut;

import static java.lang.Math.exp;
import static java.lang.Math.log;

import java.util.List;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.jujutsu.utils.MatrixOps;

public class ParallelBHTsne extends BHTSne {

	@Override
	void updateGradient(int N, int no_dims, double[] Y, double momentum, double eta, double[] dY, double[] uY,
			double[] gains) {
		IntStream.range(0, N * no_dims).parallel().forEach(i -> {
			// Update gains
			gains[i] = (sign_tsne(dY[i]) != sign_tsne(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
			if(gains[i] < .01) gains[i] = .01;

			// Perform gradient update (with momentum and gains)
			Y[i] = Y[i] + uY[i];
			uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
		});
	}
	// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
	@Override
	void computeGradient(double [] P, int [] inp_row_P, 
			int [] inp_col_P, double [] inp_val_P, 	double [] Y, int N, int D, 
			double [] dC, double theta)
	{
		// Construct space-partitioning tree on current map
		ParallelSPTree tree = new ParallelSPTree(D, Y, N);

		double totalSum_Q = 0.0;
		double[] sum_Q = new double[N];
		double[] pos_f = new double[N * D];
		double[][] neg_f = new double[N][D];
		double[][] buff = new double[N][D];
		tree.computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, N, pos_f);

		IntStream.range(0, N).parallel().forEach(i -> {
			tree.computeNonEdgeForces(i, theta, neg_f[i], buff[i], sum_Q);
		});
		totalSum_Q = DoubleStream.of(sum_Q).sum();

		// Compute final t-SNE gradient
		for (int n = 0; n < N; n++)
		{
			for (int d = 0; d < D; d++)
			{
				dC[n * D + d] = pos_f[n * D + d] - (neg_f[n][d] / totalSum_Q);
			}
		}
	}

	@Override
	// Compute input similarities with a fixed perplexity using ball trees
	void computeGaussianPerplexity(double [] X, int N, int D, int [] _row_P, int [] _col_P, double [] _val_P, double perplexity, int K) {
		if(perplexity > K) System.out.println("Perplexity should be lower than K!");

		// Allocate the memory we need
		/**_row_P = (int*)    malloc((N + 1) * sizeof(int));
		 *_col_P = (int*)    calloc(N * K, sizeof(int));
		 *_val_P = (double*) calloc(N * K, sizeof(double));
			    if(*_row_P == null || *_col_P == null || *_val_P == null) { Rcpp::stop("Memory allocation failed!\n"); }*/
		int [] row_P = _row_P;
		int [] col_P = _col_P;
		double [] val_P = _val_P;
		double [] cur_P = new double[N - 1];

		row_P[0] = 0;
		for(int n = 0; n < N; n++) row_P[n + 1] = row_P[n] + K;    

		// Build ball tree on data set
		ParallelVpTree<DataPoint> tree = new ParallelVpTree<DataPoint>(distance);
		final DataPoint [] obj_X = new DataPoint [N];
		for(int n = 0; n < N; n++) {
			double [] row = MatrixOps.extractRowFromFlatMatrix(X,n,D);
			obj_X[n] = new DataPoint(D, n, row);
		}
		tree.create(obj_X);

		// VERIFIED THAT TREES LOOK THE SAME
		//System.out.println("Created Tree is: ");
		//			AdditionalInfoProvider pp = new AdditionalInfoProvider() {			
		//				@Override
		//				public String provideInfo(Node node) {
		//					return "" + obj_X[node.index].index();
		//				}
		//			};
		//			TreePrinter printer = new TreePrinter(pp);
		//			printer.printTreeHorizontal(tree.getRoot());

		// Loop over all points to find nearest neighbors
		List<TreeSearchResult> results = tree.searchMultiple(tree, obj_X, K+1);

		for (TreeSearchResult res : results) {
			List<Double> distances = res.getDistances();
			List<DataPoint> indices = res.getIndices();
			int n = res.getIndex();

			// Initialize some variables for binary search
			boolean found = false;
			double beta = 1.0;
			double min_beta = -Double.MAX_VALUE;
			double max_beta =  Double.MAX_VALUE;
			double tol = 1e-5;

			// Iterate until we found a good perplexity
			int iter = 0; 
			double sum_P = 0.;
			while(!found && iter < 200) {

				// Compute Gaussian kernel row and entropy of current row
				sum_P = Double.MIN_VALUE;
				double H = .0;
				for(int m = 0; m < K; m++) {
					cur_P[m] = exp(-beta * distances.get(m + 1));
					sum_P += cur_P[m];
					H += beta * (distances.get(m + 1) * cur_P[m]);
				}
				H = (H / sum_P) + log(sum_P);

				// Evaluate whether the entropy is within the tolerance level
				double Hdiff = H - log(perplexity);
				if(Hdiff < tol && -Hdiff < tol) {
					found = true;
				}
				else {
					if(Hdiff > 0) {
						min_beta = beta;
						if(max_beta == Double.MAX_VALUE || max_beta == -Double.MAX_VALUE)
							beta *= 2.0;
						else
							beta = (beta + max_beta) / 2.0;
					}
					else {
						max_beta = beta;
						if(min_beta == -Double.MAX_VALUE || min_beta == Double.MAX_VALUE)
							beta /= 2.0;
						else
							beta = (beta + min_beta) / 2.0;
					}
				}

				// Update iteration counter
				iter++;
			}

			// Row-normalize current row of P and store in matrix 
			for(int m = 0; m < K; m++) {
				cur_P[m] /= sum_P;
				col_P[row_P[n] + m] = indices.get(m + 1).index();
				val_P[row_P[n] + m] = cur_P[m];
			}
		}
	}

}
