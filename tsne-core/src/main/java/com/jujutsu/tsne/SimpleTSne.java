package com.jujutsu.tsne;

import static com.jujutsu.utils.MatrixOps.abs;
import static com.jujutsu.utils.MatrixOps.addColumnVector;
import static com.jujutsu.utils.MatrixOps.addRowVector;
import static com.jujutsu.utils.MatrixOps.assignAllLessThan;
import static com.jujutsu.utils.MatrixOps.assignAtIndex;
import static com.jujutsu.utils.MatrixOps.assignValuesToRow;
import static com.jujutsu.utils.MatrixOps.biggerThan;
import static com.jujutsu.utils.MatrixOps.concatenate;
import static com.jujutsu.utils.MatrixOps.diag;
import static com.jujutsu.utils.MatrixOps.equal;
import static com.jujutsu.utils.MatrixOps.exp;
import static com.jujutsu.utils.MatrixOps.fillMatrix;
import static com.jujutsu.utils.MatrixOps.getValuesFromRow;
import static com.jujutsu.utils.MatrixOps.log;
import static com.jujutsu.utils.MatrixOps.maximum;
import static com.jujutsu.utils.MatrixOps.mean;
import static com.jujutsu.utils.MatrixOps.negate;
import static com.jujutsu.utils.MatrixOps.plus;
import static com.jujutsu.utils.MatrixOps.range;
import static com.jujutsu.utils.MatrixOps.replaceNaN;
import static com.jujutsu.utils.MatrixOps.rnorm;
import static com.jujutsu.utils.MatrixOps.scalarDivide;
import static com.jujutsu.utils.MatrixOps.scalarInverse;
import static com.jujutsu.utils.MatrixOps.scalarMult;
import static com.jujutsu.utils.MatrixOps.scalarPlus;
import static com.jujutsu.utils.MatrixOps.sqrt;
import static com.jujutsu.utils.MatrixOps.square;
import static com.jujutsu.utils.MatrixOps.sum;
import static com.jujutsu.utils.MatrixOps.tile;
import static com.jujutsu.utils.MatrixOps.times;

import com.jujutsu.utils.MatrixOps;

/**
*
* Author: Leif Jonsson (leif.jonsson@gmail.com)
* 
* This is a Java implementation of van der Maaten and Hintons t-sne 
* dimensionality reduction technique that is particularly well suited 
* for the visualization of high-dimensional datasets
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
		String IMPLEMENTATION_NAME = this.getClass().getSimpleName();
		System.out.println("X:Shape is = " + X.length + " x " + X[0].length);
		System.out.println("Running " + IMPLEMENTATION_NAME + ".");
		// Initialize variables
		if(use_pca && X[0].length > initial_dims && initial_dims > 0) {
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
		double [][] Y           = rnorm(n,no_dims);
		double [][] dY          = fillMatrix(n,no_dims,0.0);
		double [][] iY          = fillMatrix(n,no_dims,0.0);
		double [][] gains       = fillMatrix(n,no_dims,1.0);
		
		// Compute P-values
		double [][] P = x2p(X, 1e-5, perplexity).P;
		P = plus(P , mo.transpose(P));
		P = scalarDivide(P,sum(P));
		P = scalarMult(P , 4);					// early exaggeration
		P = maximum(P, 1e-12);

		System.out.println("Y:Shape is = " + Y.length + " x " + Y[0].length);
		
		// Run iterations
		for (int iter = 0; iter < max_iter; iter++) {
			// Compute pairwise affinities
			double [][] sum_Y = mo.transpose(sum(square(Y), 1));
			double [][] num = scalarInverse(scalarPlus(addRowVector(mo.transpose(addRowVector(scalarMult(
					times(Y, mo.transpose(Y)),
					-2),
					sum_Y)),
					sum_Y),
					1));
			assignAtIndex(num, range(n), range(n), 0);
			double [][] Q = scalarDivide(num , sum(num));

			Q = maximum(Q, 1e-12);

			// Compute gradient
			double[][] L = mo.scalarMultiply(mo.minus(P , Q), num);
		    dY = scalarMult(times(mo.minus(diag(sum(L, 1)),L) , Y), 4);
			
			// Perform the update
			if (iter < 20)
				momentum = initial_momentum;
			else
				momentum = final_momentum;
			gains = plus(mo.scalarMultiply(scalarPlus(gains,.2), abs(negate(equal(biggerThan(dY,0.0),biggerThan(iY,0.0))))),
					mo.scalarMultiply(scalarMult(gains,.8), abs(equal(biggerThan(dY,0.0),biggerThan(iY,0.0)))));

			assignAllLessThan(gains, min_gain, min_gain);
			iY = mo.minus(scalarMult(iY,momentum) , scalarMult(mo.scalarMultiply(gains , dY),eta));
			Y = plus(Y , iY);
			//double [][] tile = tile(mean(Y, 0), n, 1);
			Y = mo.minus(Y , tile(mean(Y, 0), n, 1));

			// Compute current value of cost function
			if ((iter % 100 == 0))   {
				double [][] logdivide = log(scalarDivide(P , Q));
				logdivide = replaceNaN(logdivide,0);
				double C = sum(mo.scalarMultiply(P , logdivide));
				System.out.println("Iteration " + (iter + 1) + ": error is " + C);
			} else if((iter + 1) % 10 == 0) {
				System.out.println("Iteration " + (iter + 1));
			}

			// Stop lying about P-values
			if (iter == 100)
				P = scalarDivide(P , 4);
		}

		// Return solution
		return Y;
	}

	public R Hbeta (double [][] D, double beta){
		double [][] P = exp(scalarMult(scalarMult(D,beta),-1));
		double sumP = sum(P);   // sumP confirmed scalar
		double H = Math.log(sumP) + beta * sum(mo.scalarMultiply(D,P)) / sumP;
		P = scalarDivide(P,sumP);
		R r = new R();
		r.H = H;
		r.P = P;
		return r;
	}

	public R x2p(double [][] X,double tol, double perplexity){
		int n               = X.length;
		double [][] sum_X   = sum(square(X), 1);
		double [][] times   = scalarMult(times(X, mo.transpose(X)), -2);
		double [][] prodSum = addColumnVector(mo.transpose(times), sum_X);
		double [][] D       = addRowVector(prodSum, mo.transpose(sum_X));
		// D seems correct at this point compared to Python version
		double [][] P       = fillMatrix(n,n,0.0);
		double [] beta      = fillMatrix(n,n,1.0)[0];
		double logU         = Math.log(perplexity);
		System.out.println("Starting x2p...");
		for (int i = 0; i < n; i++) {
			if (i % 500 == 0)
				System.out.println("Computing P-values for point " + i + " of " + n + "...");
			double betamin = Double.NEGATIVE_INFINITY;
			double betamax = Double.POSITIVE_INFINITY;
			double [][] Di = getValuesFromRow(D, i,concatenate(range(0,i),range(i+1,n)));

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
			assignValuesToRow(P, i,concatenate(range(0,i),range(i+1,n)),thisP[0]);
		}

		R r = new R();
		r.P = P;
		r.beta = beta;
		double sigma = mean(sqrt(scalarInverse(beta)));

		System.out.println("Mean value of sigma: " + sigma);

		return r;
	}
}
