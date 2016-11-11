/*
 *
 * This Java port of Barnes Hut t-SNE is Copyright (c) Leif Jonsson 2016 and 
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */
package com.jujutsu.tsne.barneshut;

import static java.lang.Math.exp;
import static java.lang.Math.log;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

import com.jujutsu.tsne.PrincipalComponentAnalysis;
import com.jujutsu.utils.MatrixOps;

public class BHTSne implements BarnesHutTSne {

	protected final Distance distance = new EuclideanDistance();

	@Override
	public double[][] tsne(double[][] X, int no_dims, int initial_dims, double perplexity) {
		return tsne(X,no_dims,initial_dims,perplexity,20000,true);
	}

	@Override
	public double[][] tsne(double[][] X, int no_dims, int initial_dims, double perplexity, int maxIterations) {
		return tsne(X,no_dims,initial_dims,perplexity,maxIterations,true);
	}

	@Override
	public double[][] tsne(double[][] X, int no_dims, int initial_dims, double perplexity, int max_iter, boolean use_pca) {
		return tsne(X,no_dims,initial_dims,perplexity,max_iter,use_pca, 0.5);
	}

	@Override
	public double[][] tsne(double[][] X, int no_dims, int initial_dims, double perplexity, int max_iter, boolean use_pca, double theta) {
		int N = X.length;
		int D = X[0].length;
		return run(X, N, D, no_dims, initial_dims, perplexity, max_iter, use_pca, theta);
	}

	private double[] flatten(double[][] x) {
		int noCols = x[0].length;
		double [] flat = new double[x.length*x[0].length];
		for (int i = 0; i < x.length; i++) {
			for (int j = 0; j < x[i].length; j++) {
				flat[i*noCols+j] = x[i][j];
			}
		}
		return flat;
	}

	private double [][] expand(double[]x, int N, int D) {
		double [][] expanded = new double[N][D];
		for (int row = 0; row < N; row++) {
			for (int col = 0; col < D; col++) {
				expanded[row][col] = x[row*D+col];
			}
		}
		return expanded;
	}

	static double sign_tsne(double x) { return (x == .0 ? .0 : (x < .0 ? -1.0 : 1.0)); }

	// Perform t-SNE
	double [][] run(double [][] Xin, int N, int D, int no_dims, int initial_dims, double perplexity, 
			int max_iter, boolean use_pca, double theta) {
		boolean exact = (theta == .0) ? true : false;
		if(exact) throw new IllegalArgumentException("The Barnes Hut implementation does not support exact inference yet (theta==0.0), if you want exact t-SNE please use one of the standard t-SNE implementations (FastTSne for instance)");
		
		if(use_pca && D > initial_dims && initial_dims > 0) {
			PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
			Xin = pca.pca(Xin, initial_dims);
			D = initial_dims;
			System.out.println("X:Shape after PCA is = " + Xin.length + " x " + Xin[0].length);
		}
		double [] X = flatten(Xin);	
		
		double [] Y = new double[N*no_dims];
		System.out.println("X:Shape is = " + N + " x " + D);
		// Determine whether we are using an exact algorithm
		if(N - 1 < 3 * perplexity) { throw new IllegalArgumentException("Perplexity too large for the number of data points!\n"); }
		System.out.printf("Using no_dims = %d, perplexity = %f, and theta = %f\n", no_dims, perplexity, theta);

		// Set learning parameters
		double total_time = 0;
		int stop_lying_iter = 250, mom_switch_iter = 250;
		double momentum = .5, final_momentum = .8;
		double eta = 200.0;

		// Allocate some memory
		double [] dY    = new double[N * no_dims];
		double [] uY    = new double[N * no_dims];
		double [] gains = new double[N * no_dims];
		for(int i = 0; i < N * no_dims; i++) gains[i] = 1.0;

		// Normalize input data (to prevent numerical problems)
		System.out.println("Computing input similarities...");
		long start = System.currentTimeMillis();
		//zeroMean(X, N, D);
		double max_X = .0;
		for(int i = 0; i < N * D; i++) {
			if(X[i] > max_X) max_X = X[i];
		}

		for(int i = 0; i < N * D; i++) X[i] /= max_X;

		double [] P = null;
		int K  = (int) (3 * perplexity);
		int [] row_P = new int[N+1];
		int [] col_P = new int[N*K];
		double [] val_P = new double[N*K];
		/**_row_P = (int*)    malloc((N + 1) * sizeof(int));
		 *_col_P = (int*)    calloc(N * K, sizeof(int));
		 *_val_P = (double*) calloc(N * K, sizeof(double));*/
		// Compute input similarities for exact t-SNE
		if(exact) {

			// Compute similarities
			P = new double[N * N];
			computeGaussianPerplexity(X, N, D, P, perplexity);

			// Symmetrize input similarities
			System.out.println("Symmetrizing...");
			int nN = 0;
			for(int n = 0; n < N; n++) {
				int mN = 0;
				for(int m = n + 1; m < N; m++) {
					P[nN + m] += P[mN + n];
					P[mN + n]  = P[nN + m];
					mN += N;
				}
				nN += N;
			}
			double sum_P = .0;
			for(int i = 0; i < N * N; i++) sum_P += P[i];
			for(int i = 0; i < N * N; i++) P[i] /= sum_P;
		}

		// Compute input similarities for approximate t-SNE
		else {

			// Compute asymmetric pairwise input similarities
			computeGaussianPerplexity(X, N, D, row_P, col_P, val_P, perplexity, K);

			// Verified that val_P,col_P,row_P is the same at this point

			// Symmetrize input similarities
			SymResult res = symmetrizeMatrix(row_P, col_P, val_P, N);
			row_P = res.sym_row_P;
			col_P = res.sym_col_P;
			val_P = res.sym_val_P;

			double sum_P = .0;
			for(int i = 0; i < row_P[N]; i++) sum_P += val_P[i];
			for(int i = 0; i < row_P[N]; i++) val_P[i] /= sum_P;
		}
		long end = System.currentTimeMillis();

		// Lie about the P-values
		if(exact) { for(int i = 0; i < N * N; i++)        P[i] *= 12.0; }
		else {      for(int i = 0; i < row_P[N]; i++) val_P[i] *= 12.0; }

		// Initialize solution (randomly)
		for(int i = 0; i < N * no_dims; i++) Y[i] = ThreadLocalRandom.current().nextDouble() * 0.0001;

		// Perform main training loop
		if(exact) System.out.printf("Done in %4.2f seconds!\nLearning embedding...\n", (end - start) / 1000.0);
		else System.out.printf("Done in %4.2f seconds (sparsity = %f)!\nLearning embedding...\n", (end - start) / 1000.0, (double) row_P[N] / ((double) N * (double) N));
		start = System.currentTimeMillis();
		for(int iter = 0; iter < max_iter; iter++) {
			
			if(exact) computeExactGradient(P, Y, N, no_dims, dY);
			// Compute (approximate) gradient
			else computeGradient(P, row_P, col_P, val_P, Y, N, no_dims, dY, theta);
			
			updateGradient(N, no_dims, Y, momentum, eta, dY, uY, gains);

			// Make solution zero-mean
			zeroMean(Y, N, no_dims);

			// Stop lying about the P-values after a while, and switch momentum
			if(iter == stop_lying_iter) {
				if(exact) { for(int i = 0; i < N * N; i++)        P[i] /= 12.0; }
				else      { for(int i = 0; i < row_P[N]; i++) val_P[i] /= 12.0; }
			}
			if(iter == mom_switch_iter) momentum = final_momentum;

			// Print out progress
			if((iter > 0 && iter % 50 == 0) || iter == max_iter - 1) {
				end = System.currentTimeMillis();
				double C = .0;
				if(exact) C = evaluateError(P, Y, N, no_dims);
				else      C = evaluateError(row_P, col_P, val_P, Y, N, no_dims, theta);  // doing approximate computation here!
				if(iter == 0)
					System.out.printf("Iteration %d: error is %f\n", iter + 1, C);
				else {
					total_time += (end - start) / 1000.0;
					System.out.printf("Iteration %d: error is %f (50 iterations in %4.2f seconds)\n", iter, C, (end - start) / 1000.0);
				}
				start = System.currentTimeMillis();
			}
		}
		end = System.currentTimeMillis(); total_time += (end - start) / 1000.0;

		System.out.printf("Fitting performed in %4.2f seconds.\n", total_time);
		return expand(Y,N,no_dims);
	}

	void updateGradient(int N, int no_dims, double[] Y, double momentum, double eta, double[] dY, double[] uY,
			double[] gains) {
		for(int i = 0; i < N * no_dims; i++)  {
			// Update gains
			gains[i] = (sign_tsne(dY[i]) != sign_tsne(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
			if(gains[i] < .01) gains[i] = .01;

			// Perform gradient update (with momentum and gains)
			Y[i] = Y[i] + uY[i];
			uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
		}
	}

	// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
	void computeGradient(double [] P, int [] inp_row_P, int [] inp_col_P, double [] inp_val_P, double [] Y, int N, int D, double [] dC, double theta)
	{
		// Construct space-partitioning tree on current map
		SPTree tree = new SPTree(D, Y, N);

		// Compute all terms required for t-SNE gradient
		double [] sum_Q = new double[1];
		double [] pos_f = new double[N * D];
		double [][] neg_f = new double[N][D];

		tree.computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, N, pos_f);
		for(int n = 0; n < N; n++) tree.computeNonEdgeForces(n, theta, neg_f[n], sum_Q);

		// Compute final t-SNE gradient
		for(int n = 0; n < N; n++) {
			for(int d = 0; d < D; d++) {
				dC[n*D+d] = pos_f[n*D+d] - (neg_f[n][d] / sum_Q[0]);
			}
		}
	}

	// Compute gradient of the t-SNE cost function (exact)
	void computeExactGradient(double [] P, double [] Y, int N, int D, double [] dC) {

		// Make sure the current gradient contains zeros
		for(int i = 0; i < N * D; i++) dC[i] = 0.0;

		// Compute the squared Euclidean distance matrix
		double [] DD = new double[N * N];
		computeSquaredEuclideanDistance(Y, N, D, DD);

		// Compute Q-matrix and normalization sum
		double [] Q    = new double[N * N];
		double sum_Q = .0;
		int nN = 0;
		for(int n = 0; n < N; n++) {
			for(int m = 0; m < N; m++) {
				if(n != m) {
					Q[nN + m] = 1 / (1 + DD[nN + m]);
					sum_Q += Q[nN + m];
				}
			}
			nN += N;
		}

		// Perform the computation of the gradient
		nN = 0;
		int nD = 0;
		for(int n = 0; n < N; n++) {
			int mD = 0;
			for(int m = 0; m < N; m++) {
				if(n != m) {
					double mult = (P[nN + m] - (Q[nN + m] / sum_Q)) * Q[nN + m];
					for(int d = 0; d < D; d++) {
						dC[nD + d] += (Y[nD + d] - Y[mD + d]) * mult;
					}
				}
				mD += D;
			}
			nN += N;
			nD += D;
		}
	}

	// Evaluate t-SNE cost function (exactly)
	double evaluateError(double [] P, double [] Y, int N, int D) {

		// Compute the squared Euclidean distance matrix
		double [] DD = new double[N * N];
		double [] Q = new double[N * N];
		computeSquaredEuclideanDistance(Y, N, D, DD);

		// Compute Q-matrix and normalization sum
		int nN = 0;
		double sum_Q = Double.MIN_VALUE;
		for(int n = 0; n < N; n++) {
			for(int m = 0; m < N; m++) {
				if(n != m) {
					Q[nN + m] = 1 / (1 + DD[nN + m]);
					sum_Q += Q[nN + m];
				}
				else Q[nN + m] = Double.MIN_VALUE;
			}
			nN += N;
		}
		for(int i = 0; i < N * N; i++) Q[i] /= sum_Q;

		// Sum t-SNE error
		double C = .0;
		for(int n = 0; n < N * N; n++) {
			C += P[n] * log((P[n] + Double.MIN_VALUE) / (Q[n] + Double.MIN_VALUE));
		}

		return C;
	}

	// Evaluate t-SNE cost function (approximately)
	double evaluateError(int [] row_P, int [] col_P, double [] val_P, double [] Y, int N, int D, double theta)
	{
		// Get estimate of normalization term
		SPTree tree = new SPTree(D, Y, N);
		double [] buff = new double[D];
		double [] sum_Q = new double[1];
		for(int n = 0; n < N; n++) tree.computeNonEdgeForces(n, theta, buff, sum_Q);

		// Loop over all edges to compute t-SNE error
		int ind1, ind2;
		double C = .0, Q;
		for(int n = 0; n < N; n++) {
			ind1 = n * D;
			for(int i = row_P[n]; i < row_P[n + 1]; i++) {
				Q = .0;
				ind2 = col_P[i] * D;
				for(int d = 0; d < D; d++) buff[d]  = Y[ind1 + d];
				for(int d = 0; d < D; d++) buff[d] -= Y[ind2 + d];
				for(int d = 0; d < D; d++) Q += buff[d] * buff[d];
				Q = (1.0 / (1.0 + Q)) / sum_Q[0];
				C += val_P[i] * log((val_P[i] + Double.MIN_VALUE) / (Q + Double.MIN_VALUE));
			}
		}

		return C;
	}

	// Compute input similarities with a fixed perplexity
	void computeGaussianPerplexity(double [] X, int N, int D, double [] P, double perplexity) {

		// Compute the squared Euclidean distance matrix
		double [] DD = new double[N * N];
		computeSquaredEuclideanDistance(X, N, D, DD);

		// Compute the Gaussian kernel row by row
		int nN = 0;
		for(int n = 0; n < N; n++) {

			// Initialize some variables
			boolean found = false;
			double beta = 1.0;
			double min_beta = -Double.MAX_VALUE;
			double max_beta =  Double.MAX_VALUE;
			double tol = 1e-5;
			double sum_P = Double.MIN_VALUE;

			// Iterate until we found a good perplexity
			int iter = 0;
			while(!found && iter < 200) {

				// Compute Gaussian kernel row
				for(int m = 0; m < N; m++) P[nN + m] = exp(-beta * DD[nN + m]);
				P[nN + n] = Double.MIN_VALUE;

				// Compute entropy of current row
				sum_P = Double.MIN_VALUE;
				for(int m = 0; m < N; m++) sum_P += P[nN + m];
				double H = 0.0;
				for(int m = 0; m < N; m++) H += beta * (DD[nN + m] * P[nN + m]);
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

			// Row normalize P
			for(int m = 0; m < N; m++) P[nN + m] /= sum_P;
			nN += N;
		}
	}

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
		VpTree<DataPoint> tree = new VpTree<DataPoint>(distance);
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
		System.out.println("Building tree...");
		List<DataPoint> indices = new ArrayList<>();
		List<Double> distances = new ArrayList<>();
		for(int n = 0; n < N; n++) {

			if(n % 10000 == 0) System.out.printf(" - point %d of %d\n", n, N);

			// Find nearest neighbors
			indices.clear();
			distances.clear();
			//System.out.println("Looking at: " + obj_X.get(n).index());
			tree.search(obj_X[n], K + 1, indices, distances);

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

	void computeGaussianPerplexity(double [] X, int N, int D, int [] _row_P, int [] _col_P, double [] _val_P, double perplexity, double threshold) {
		// Allocate some memory we need for computations
		double [] buff  = new double[D];
		double [] DD    = new double[N];
		double [] cur_P = new double[N];

		// Compute the Gaussian kernel row by row (to find number of elements in sparse P)
		int total_count = 0;
		for(int n = 0; n < N; n++) {

			// Compute the squared Euclidean distance matrix
			for(int m = 0; m < N; m++) {
				for(int d = 0; d < D; d++) buff[d]  = X[n * D + d];
				for(int d = 0; d < D; d++) buff[d] -= X[m * D + d];
				DD[m] = .0;
				for(int d = 0; d < D; d++) DD[m] += buff[d] * buff[d];
			}

			// Initialize some variables
			boolean found = false;
			double beta = 1.0;
			double min_beta = -Double.MAX_VALUE;
			double max_beta =  Double.MAX_VALUE;
			double tol = 1e-5;

			// Iterate until we found a good perplexity
			int iter = 0; 
			double sum_P = 0.;
			while(!found && iter < 200) {

				// Compute Gaussian kernel row
				for(int m = 0; m < N; m++) cur_P[m] = exp(-beta * DD[m]);
				cur_P[n] = Double.MIN_VALUE;

				// Compute entropy of current row
				sum_P = Double.MIN_VALUE;
				for(int m = 0; m < N; m++) sum_P += cur_P[m];
				double H = 0.0;
				for(int m = 0; m < N; m++) H += beta * (DD[m] * cur_P[m]);
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

			// Row-normalize and threshold current row of P
			for(int m = 0; m < N; m++) cur_P[m] /= sum_P;
			for(int m = 0; m < N; m++) {
				if(cur_P[m] > threshold / (double) N) total_count++;
			}
		}

		// Allocate the memory we need
		_row_P = new int[N + 1];
		_col_P = new int[total_count];
		_val_P = new double[total_count];
		int [] row_P = _row_P;
		int [] col_P = _col_P;
		double [] val_P = _val_P;
		row_P[0] = 0;

		// Compute the Gaussian kernel row by row (this time, store the results)
		int count = 0;
		for(int n = 0; n < N; n++) {

			// Compute the squared Euclidean distance matrix
			for(int m = 0; m < N; m++) {
				for(int d = 0; d < D; d++) buff[d]  = X[n * D + d];
				for(int d = 0; d < D; d++) buff[d] -= X[m * D + d];
				DD[m] = .0;
				for(int d = 0; d < D; d++) DD[m] += buff[d] * buff[d];
			}

			// Initialize some variables
			boolean found = false;
			double beta = 1.0;
			double min_beta = -Double.MAX_VALUE;
			double max_beta =  Double.MAX_VALUE;
			double tol = 1e-5;

			// Iterate until we found a good perplexity
			int iter = 0; 
			double sum_P = 0.;
			while(!found && iter < 200) {

				// Compute Gaussian kernel row
				for(int m = 0; m < N; m++) cur_P[m] = exp(-beta * DD[m]);
				cur_P[n] = Double.MIN_VALUE;

				// Compute entropy of current row
				sum_P = Double.MIN_VALUE;
				for(int m = 0; m < N; m++) sum_P += cur_P[m];
				double H = 0.0;
				for(int m = 0; m < N; m++) H += beta * (DD[m] * cur_P[m]);
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

			// Row-normalize and threshold current row of P
			for(int m = 0; m < N; m++) cur_P[m] /= sum_P;
			for(int m = 0; m < N; m++) {
				if(cur_P[m] > threshold / (double) N) {
					col_P[count] = m;
					val_P[count] = cur_P[m];
					count++;
				}
			}
			row_P[n + 1] = count;
		}

	}

	class SymResult {
		int []    sym_row_P;
		int []    sym_col_P;
		double [] sym_val_P;

		public SymResult(int[] sym_row_P, int[] sym_col_P, double[] sym_val_P) {
			super();
			this.sym_row_P = sym_row_P;
			this.sym_col_P = sym_col_P;
			this.sym_val_P = sym_val_P;
		}
	}

	SymResult symmetrizeMatrix(int [] _row_P, int [] _col_P, double [] _val_P, int N) {

		// Get sparse matrix
		int [] row_P = _row_P;
		int [] col_P = _col_P;
		double [] val_P = _val_P;

		// Count number of elements and row counts of symmetric matrix
		int [] row_counts = new int[N];
		for(int n = 0; n < N; n++) {
			for(int i = row_P[n]; i < row_P[n + 1]; i++) {

				// Check whether element (col_P[i], n) is present
				boolean present = false;
				for(int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
					if(col_P[m] == n) present = true;
				}
				if(present) row_counts[n]++;
				else {
					row_counts[n]++;
					row_counts[col_P[i]]++;
				}
			}
		}
		int no_elem = 0;
		for(int n = 0; n < N; n++) no_elem += row_counts[n];

		// Allocate memory for symmetrized matrix
		int []    sym_row_P = new int[N + 1];
		int []    sym_col_P = new int[no_elem];
		double [] sym_val_P = new double[no_elem];

		// Construct new row indices for symmetric matrix
		sym_row_P[0] = 0;
		for(int n = 0; n < N; n++) sym_row_P[n + 1] = sym_row_P[n] + row_counts[n];

		// Fill the result matrix
		int [] offset = new int[N];
		for(int n = 0; n < N; n++) {
			for(int i = row_P[n]; i < row_P[n + 1]; i++) {                                  // considering element(n, col_P[i])

				// Check whether element (col_P[i], n) is present
				boolean present = false;
				for(int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
					if(col_P[m] == n) {
						present = true;
						if(n <= col_P[i]) {                                                 // make sure we do not add elements twice
							sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
							sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
							sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i] + val_P[m];
							sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i] + val_P[m];
						}
					}
				}

				// If (col_P[i], n) is not present, there is no addition involved
				if(!present) {
					sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
					sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
					sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i];
					sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i];
				}

				// Update offsets
				if(!present || (present && n <= col_P[i])) {
					offset[n]++;
					if(col_P[i] != n) offset[col_P[i]]++;               
				}
			}
		}

		// Divide the result by two
		for(int i = 0; i < no_elem; i++) sym_val_P[i] /= 2.0;

		return new SymResult(sym_row_P, sym_col_P, sym_val_P);
	}

	// Compute squared Euclidean distance matrix (using BLAS)
	void computeSquaredEuclideanDistance(double [] X, int N, int D, double [] DD) {
		double [] dataSums = new double[N];
		for(int n = 0; n < N; n++) {
			for(int d = 0; d < D; d++) {
				dataSums[n] += (X[n * D + d] * X[n * D + d]);
			}
		}
		for(int n = 0; n < N; n++) {
			for(int m = 0; m < N; m++) {
				DD[n * N + m] = dataSums[n] + dataSums[m];
			}
		}
		//double a1 = -2.0;
		//double a2 = 1.0;
		//dgemm(String arg0, String arg1, int arg2, int arg3, int arg4, double arg5, double[] arg6, int arg7, int arg8, double[] arg9, int arg10, int arg11, double arg12, double[] arg13, int arg14, int arg15);
		//  org.netlib.blas.Dgemm.dgemm("T", "N", N, N, D, a1, X, D, X, D, a2, DD, N);

		/* DGEMM - perform one of the matrix-matrix operations    */
		/* C := alpha*op( A )*op( B ) + beta*C */
		/* BLAS_extern void
		 * R BLAS declaration:
		    F77_NAME(dgemm)(const char *transa, const char *transb, const int *m,
		    		const int *n, const int *k, const double *alpha,
		    		const double *a, const int *lda,
		    		const double *b, const int *ldb,
		    		const double *beta, double *c, const int *ldc);*/

		//org.netlib.blas.Dgemm.dgemm
	}


	// Makes data zero-mean
	void zeroMean(double [] X, int N, int D) {

		// Compute data mean
		double [] mean = new double[D];
		for(int n = 0; n < N; n++) {
			for(int d = 0; d < D; d++) {
				mean[d] += X[n * D + d];
			}
		}
		for(int d = 0; d < D; d++) {
			mean[d] /= (double) N;
		}

		// Subtract data mean
		for(int n = 0; n < N; n++) {
			for(int d = 0; d < D; d++) {
				X[n * D + d] -= mean[d];
			}
		}
	}

	// Makes data zero-mean
	void zeroMean(double [][] X, int N, int D) {

		// Compute data mean
		double [] mean = new double[D];
		for(int n = 0; n < N; n++) {
			for(int d = 0; d < D; d++) {
				mean[d] += X[n][d];
			}
		}
		for(int d = 0; d < D; d++) {
			mean[d] /= (double) N;
		}

		// Subtract data mean
		for(int n = 0; n < N; n++) {
			for(int d = 0; d < D; d++) {
				X[n][d] -= mean[d];
			}
		}
	}

}