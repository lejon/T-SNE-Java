package com.jujutsu.tsne;

/**
 *
 * Author: Leif Jonsson (leif.jonsson@gmail.com)
 * 
 * This is a port of van der Maaten and Hintons Python implementation of t-sne
 *
 */
public class SimpleTSne implements TSne {
	MatrixOps mo = new MatrixOps();

	public double [][] tsne(double[][] X, int k, int initial_dims, double perplexity) {
		return tsne(X,k,initial_dims, perplexity, 2000, true);
	}

	public double [][] tsne(double[][] X, int k, int initial_dims, double perplexity, int maxIterations) {
		return tsne(X,k,initial_dims, perplexity, maxIterations, true);
	}

	public double [][] tsne(double[][] X, int no_dims, int initial_dims, double perplexity, int max_iter, boolean use_pca) {
		System.out.println("X:Shape is = " + X.length + " x " + X[0].length);
		// Initialize variables
		if(use_pca && X[0].length > initial_dims) {
			PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
			X = pca.pca(X, initial_dims);
			System.out.println("X:Shape after PCA is = " + X.length + " x " + X[0].length);
		}
		int n = X.length;
		double momentum = .5;
		double initial_momentum = 0.5;
		double final_momentum   = 0.8;
		int eta                 = 500;
		double min_gain         = 0.01;
		double [][] Y           = mo.rnorm(n,no_dims);
		double [][] dY          = mo.fillMatrix(n,no_dims,0.0);
		double [][] iY          = mo.fillMatrix(n,no_dims,0.0);
		double [][] gains       = mo.fillMatrix(n,no_dims,1.0);
		
		// Compute P-values
		double [][] P = x2p(X, 1e-5, perplexity).P;
		P = mo.plus(P , mo.transpose(P));
		P = mo.scalarDivide(P ,mo.sum(P));
		P = mo.scalarMult(P , 4);					// early exaggeration
		P = mo.maximum(P, 1e-12);

		System.out.println("Y:Shape is = " + Y.length + " x " + Y[0].length);
		
		// Run iterations
		for (int iter = 0; iter < max_iter; iter++) {
			// Compute pairwise affinities
			double [][] sum_Y = mo.transpose(mo.sum(mo.square(Y), 1));
			double [][] num = mo.scalarInverse(mo.scalarPlus(mo.addRowVector(mo.transpose(mo.addRowVector(mo.scalarMult(
					mo.times(Y, mo.transpose(Y)),
					-2),
					sum_Y)),
					sum_Y),
					1));
			mo.assignAtIndex(num, mo.range(n), mo.range(n), 0);
			double [][] Q = mo.scalarDivide(num , mo.sum(num));

			Q = mo.maximum(Q, 1e-12);

			// Compute gradient
			double[][] L = mo.scalarMultiply(mo.minus(P , Q), num);
		    dY = mo.scalarMult(mo.times(mo.minus(mo.diag(mo.sum(L, 1)),L) , Y), 4);
			
			// Perform the update
			if (iter < 20)
				momentum = initial_momentum;
			else
				momentum = final_momentum;
			gains = mo.plus(mo.scalarMultiply(mo.scalarPlus(gains,.2), mo.abs(mo.negate(mo.equal(mo.biggerThan(dY,0.0),mo.biggerThan(iY,0.0))))),
					mo.scalarMultiply(mo.scalarMult(gains,.8), mo.abs(mo.equal(mo.biggerThan(dY,0.0),mo.biggerThan(iY,0.0)))));

			mo.assignAllLessThan(gains, min_gain, min_gain);
			iY = mo.minus(mo.scalarMult(iY,momentum) , mo.scalarMult(mo.scalarMultiply(gains , dY),eta));
			Y = mo.plus(Y , iY);
			//double [][] tile = tile(mean(Y, 0), n, 1);
			Y = mo.minus(Y , mo.tile(MatrixOps.mean(Y, 0), n, 1));

			// Compute current value of cost function
			if ((iter % 100 == 0))   {
				double [][] logdivide = MatrixOps.log(mo.scalarDivide(P , Q));
				logdivide = MatrixOps.replaceNaN(logdivide,0);
				double C = mo.sum(mo.scalarMultiply(P , logdivide));
				System.out.println("Iteration " + (iter + 1) + ": error is " + C);
			} else if((iter + 1) % 10 == 0) {
				System.out.println("Iteration " + (iter + 1));
			}

			// Stop lying about P-values
			if (iter == 100)
				P = mo.scalarDivide(P , 4);
		}

		// Return solution
		return Y;
	}

	public R Hbeta (double [][] D, double beta){
		double [][] P = mo.exp(mo.scalarMult(mo.scalarMult(D,beta),-1));
		double sumP = mo.sum(P);   // sumP confirmed scalar
		double H = Math.log(sumP) + beta * mo.sum(mo.scalarMultiply(D,P)) / sumP;
		P = mo.scalarDivide(P,sumP);
		R r = new R();
		r.H = H;
		r.P = P;
		return r;
	}

	public R x2p(double [][] X,double tol, double perplexity){
		int n               = X.length;
		double [][] sum_X   = mo.sum(mo.square(X), 1);
		double [][] times   = mo.scalarMult(mo.times(X, mo.transpose(X)), -2);
		double [][] prodSum = mo.addColumnVector(mo.transpose(times), sum_X);
		double [][] D       = mo.addRowVector(prodSum, mo.transpose(sum_X));
		// D seems correct at this point compared to Python version
		double [][] P       = mo.fillMatrix(n,n,0.0);
		double [] beta      = mo.fillMatrix(n,n,1.0)[0];
		double logU         = Math.log(perplexity);
		System.out.println("Starting x2p...");
		for (int i = 0; i < n; i++) {
			if (i % 500 == 0)
				System.out.println("Computing P-values for point " + i + " of " + n + "...");
			double betamin = Double.NEGATIVE_INFINITY;
			double betamax = Double.POSITIVE_INFINITY;
			double [][] Di = mo.getValuesFromRow(D, i,mo.concatenate(mo.range(0,i),mo.range(i+1,n)));

			R hbeta = Hbeta(Di, beta[i]);
			double H = hbeta.H;
			double [][] thisP = hbeta.P;

			// Evaluate whether the perplexity is within tolerance
			double Hdiff = H - logU;
			int tries = 0;
			while(Math.abs(Hdiff) > tol && tries < 50){
				if (Hdiff > 0){
					betamin = beta[i];
					if (Double.isInfinite(betamax))
						beta[i] = beta[i] * 2;
					else 
						beta[i] = (beta[i] + betamax) / 2;
				} else{
					betamax = beta[i];
					if (Double.isInfinite(betamin))  
						beta[i] = beta[i] / 2;
					else 
						beta[i] = ( beta[i] + betamin) / 2;
				}

				hbeta = Hbeta(Di, beta[i]);
				H = hbeta.H;
				thisP = hbeta.P;
				Hdiff = H - logU;
				tries = tries + 1;
			}
			mo.assignValuesToRow(P, i,mo.concatenate(mo.range(0,i),mo.range(i+1,n)),thisP[0]);
		}

		R r = new R();
		r.P = P;
		r.beta = beta;
		double sigma = mo.mean(mo.sqrt(mo.scalarInverse(beta)));

		System.out.println("Mean value of sigma: " + sigma);

		return r;
	}
}
