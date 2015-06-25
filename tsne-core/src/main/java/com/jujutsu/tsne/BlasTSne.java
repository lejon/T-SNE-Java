package com.jujutsu.tsne;

import static com.jujutsu.utils.MatrixOps.addColumnVector;
import static com.jujutsu.utils.MatrixOps.addRowVector;
import static com.jujutsu.utils.MatrixOps.assignValuesToRow;
import static com.jujutsu.utils.MatrixOps.concatenate;
import static com.jujutsu.utils.MatrixOps.equal;
import static com.jujutsu.utils.MatrixOps.fillMatrix;
import static com.jujutsu.utils.MatrixOps.getValuesFromRow;
import static com.jujutsu.utils.MatrixOps.mean;
import static com.jujutsu.utils.MatrixOps.negate;
import static com.jujutsu.utils.MatrixOps.range;
import static com.jujutsu.utils.MatrixOps.scalarInverse;
import static com.jujutsu.utils.MatrixOps.scalarMult;
import static com.jujutsu.utils.MatrixOps.sqrt;
import static com.jujutsu.utils.MatrixOps.square;
import static com.jujutsu.utils.MatrixOps.sum;
import static com.jujutsu.utils.MatrixOps.times;

import org.jblas.DoubleMatrix;

import com.jujutsu.utils.BlasOps;
import com.jujutsu.utils.MatrixOps;

/**
 *
 * Author: Leif Jonsson (leif.jonsson@gmail.com)
 * 
 * This is a port of van der Maaten and Hintons Python implementation of t-sne
 *
 */
public class BlasTSne implements TSne {
	
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
		int n                    = X.length;
		double momentum          = .5;
		double initial_momentum  = 0.5;
		double final_momentum    = 0.8;
		int eta                  = 500;
		double min_gain          = 0.01;
		DoubleMatrix Y           = DoubleMatrix.randn(n,no_dims);
		DoubleMatrix dY          = DoubleMatrix.zeros(n,no_dims);
		DoubleMatrix iY          = DoubleMatrix.zeros(n,no_dims);
		DoubleMatrix gains       = DoubleMatrix.ones(n,no_dims);
		
		// Compute P-values
		double [][] Pt = x2p(X, 1e-5, perplexity).P;
		DoubleMatrix P = new DoubleMatrix(Pt);
		P = P.add(P.transpose());
		P = P.div(P.sum());
		P = P.mul(4);					// early exaggeration
		P = P.max(1e-12);

		System.out.println("Y:Shape is = " + Y.rows + " x " + Y.columns);
		
		// Run iterations
		for (int iter = 0; iter < max_iter; iter++) {
			// Compute pairwise affinities
			DoubleMatrix sum_Y = BlasOps.square(Y).rowSums().transpose();
			DoubleMatrix num = BlasOps.scalarInverse(Y.mmul(Y.transpose()).mul(-2).addRowVector(sum_Y).transpose().addRowVector(sum_Y).add(1));
			BlasOps.assignAtIndex(num, range(n), range(n), 0);
			
			DoubleMatrix Q = num.div(num.sum());
			Q = Q.max(1e-12);

			// Compute gradient
			DoubleMatrix L = P.sub(Q).mul(num);
			
		    dY = DoubleMatrix.diag(L.rowSums()).sub(L).mmul(Y).mul(4);

		    // Perform the update
			if (iter < 20)
				momentum = initial_momentum;
			else
				momentum = final_momentum;
			
			DoubleMatrix gainsSmall = new DoubleMatrix();
			gainsSmall.copy(gains);
			DoubleMatrix gainsBig = new DoubleMatrix();
			gainsBig.copy(gains);
			
			gainsSmall = gainsSmall.add(0.2);
			
			gainsBig = gainsBig.mul(0.8);
			
			DoubleMatrix btNeg = BlasOps.abs(negate(equal(BlasOps.biggerThan(dY,0.0),BlasOps.biggerThan(iY,0.0))));
			gainsSmall = gainsSmall.mul(btNeg);
			DoubleMatrix bt    = BlasOps.abs(       equal(BlasOps.biggerThan(dY,0.0),BlasOps.biggerThan(iY,0.0)));
			gainsBig = gainsBig.mul(bt);
			
			gains = gainsSmall.add(gainsBig);

			BlasOps.assignAllLessThan(gains, min_gain, min_gain);
			iY = iY.mul(momentum).sub(gains.mul(dY).mul(eta));
			Y = Y.add(iY);
			Y = Y.sub(BlasOps.tile(Y.columnMeans(), n, 1));

			// Compute current value of cost function
			if (iter % 100 == 0)   {
				DoubleMatrix logdivide = BlasOps.log(P.div(Q));
				logdivide = BlasOps.replaceNaN(logdivide,0);
				double C = P.mul(logdivide).sum();
				System.out.println("Iteration " + iter + ": error is " + C);
			} else if(iter % 10 == 0) {
				System.out.println("Iteration " + iter);
			}

			// Stop lying about P-values
			if (iter == 100)
				P = P.div(4);
		}

		// Return solution
		return Y.toArray2();
	}

	public R Hbeta (double [][] D, double beta){
		DoubleMatrix Dd = new DoubleMatrix(D);
		DoubleMatrix P = BlasOps.exp(Dd.mul(-beta));
		double sumP = P.sum();   // sumP confirmed scalar
		double H = Math.log(sumP) + beta * Dd.mul(P).sum() / sumP;
		P = P.div(sumP);
		R r = new R();
		r.H = H;
		r.P = P.toArray2();
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
