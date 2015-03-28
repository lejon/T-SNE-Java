package com.jujutsu.tsne;

import static org.ejml.ops.CommonOps.add;
import static org.ejml.ops.CommonOps.addEquals;
import static org.ejml.ops.CommonOps.divide;
import static org.ejml.ops.CommonOps.elementDiv;
import static org.ejml.ops.CommonOps.elementLog;
import static org.ejml.ops.CommonOps.elementMult;
import static org.ejml.ops.CommonOps.elementPower;
import static org.ejml.ops.CommonOps.elementSum;
import static org.ejml.ops.CommonOps.mult;
import static org.ejml.ops.CommonOps.multAddTransB;
import static org.ejml.ops.CommonOps.scale;
import static org.ejml.ops.CommonOps.subtract;
import static org.ejml.ops.CommonOps.subtractEquals;
import static org.ejml.ops.CommonOps.sumRows;
import static org.ejml.ops.CommonOps.transpose;

import org.ejml.data.DenseMatrix64F;
/**
 *
 * Author: Leif Jonsson (leif.jonsson@gmail.com)
 * 
 * This is a port of van der Maaten and Hintons Python implementation of t-sne
 *
 */
public class MemOptimizedTSne extends FastTSne {
	
	public double [][] tsne(double[][] X, int no_dims, int initial_dims, double perplexity, int max_iter, boolean use_pca) {
		System.out.println("X:Shape is = " + X.length + " x " + X[0].length);
		System.out.println("Running MemOptimized TSne.");
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
		DenseMatrix64F Y        = new DenseMatrix64F(mo.rnorm(n,no_dims));
		DenseMatrix64F Ysqlmul  = new DenseMatrix64F(Y.numRows,Y.numRows); // Ysqlmul = n x n
		DenseMatrix64F dY       = new DenseMatrix64F(mo.fillMatrix(n,no_dims,0.0));
		DenseMatrix64F iY       = new DenseMatrix64F(mo.fillMatrix(n,no_dims,0.0));
		DenseMatrix64F gains    = new DenseMatrix64F(mo.fillMatrix(n,no_dims,1.0));
		DenseMatrix64F btNeg    = new DenseMatrix64F(n,no_dims);
		DenseMatrix64F bt       = new DenseMatrix64F(n,no_dims);
		
		// Compute P-values
		DenseMatrix64F P        = new DenseMatrix64F(x2p(X, 1e-5, perplexity).P); // P = n x n
		DenseMatrix64F Psized   = new DenseMatrix64F(P.numRows,P.numCols);        // L = n x n
		DenseMatrix64F diag     = new DenseMatrix64F(mo.fillMatrix(Psized.numRows,Psized.numCols,0.0));
		
		transpose(P,Psized);
		addEquals(P,Psized);
		divide(P ,elementSum(P));
		replaceNaN(P,Double.MIN_VALUE);
		scale(4.0,P);					// early exaggeration
		maximize(P, 1e-12);
		
		System.out.println("Y:Shape is = " + Y.getNumRows() + " x " + Y.getNumCols());

		DenseMatrix64F sqed  = new DenseMatrix64F(Y.numRows,Y.numCols);  // sqed = n x n
		DenseMatrix64F sum_Y = new DenseMatrix64F(1,Y.numRows);
		DenseMatrix64F Q     = new DenseMatrix64F(P.numRows,P.numCols);  // Q = n x n
		
		for (int iter = 0; iter < max_iter; iter++) {
			// Compute pairwise affinities
			elementPower(Y, 2, sqed);
			sumRows(sqed, sum_Y);
			multAddTransB(-2.0, Y, Y, Ysqlmul);
			addRowVector(Ysqlmul, sum_Y);
			transpose(Ysqlmul);
			addRowVector(Ysqlmul, sum_Y);
			
			add(Ysqlmul, 1.0);
			divide(1.0,Ysqlmul);
			assignAtIndex(Ysqlmul, mo.range(n), mo.range(n), 0);
			divide(Ysqlmul , elementSum(Ysqlmul), Q);

			maximize(Q, 1e-12);
			
			// Compute gradient
			subtract(P, Q, Psized);
			elementMult(Psized, Ysqlmul);
			DenseMatrix64F rowsum = sumRows(Psized,null); // rowsum = nx1
			double [] rsum  = new double[rowsum.numRows];
			for (int i = 0; i < rsum.length; i++) {
				rsum[i] = rowsum.get(i,0);
			}
			setDiag(diag,rsum);
			subtract(diag, Psized, Psized);
			mult(Psized, Y, dY);
			scale(4.0, dY);
			
			// Perform the update
			if (iter < 20)
				momentum = initial_momentum;
			else
				momentum = final_momentum;
			
			boolean [][] boolMtrx = mo.equal(biggerThan(dY,0.0),biggerThan(iY,0.0));
			
			
			setData(btNeg, mo.abs(mo.negate(boolMtrx)));
			setData(bt, mo.abs(boolMtrx));
			
			DenseMatrix64F gainsSmall = new DenseMatrix64F(gains);
			DenseMatrix64F gainsBig   = new DenseMatrix64F(gains);
			add(gainsSmall,0.2);
			scale(0.8,gainsBig);
			
			elementMult(gainsSmall, btNeg);
			elementMult(gainsBig, bt);
			add(gainsSmall,gainsBig,gains);

			assignAllLessThan(gains, min_gain, min_gain);
			
			scale(momentum,iY);
			DenseMatrix64F gainsdY = new DenseMatrix64F(gains.numRows,dY.numCols);
			elementMult(gains , dY, gainsdY);
			scale(eta,gainsdY);
			subtractEquals(iY , gainsdY);
			addEquals(Y , iY);
			DenseMatrix64F colMeanY = colMean(Y, 0);
			DenseMatrix64F meanTile = tile(colMeanY, n, 1);
			subtractEquals(Y , meanTile);

			// Compute current value of the cost function
			if (iter % 100 == 0)   {
				DenseMatrix64F Pdiv = new DenseMatrix64F(P);
				elementDiv(Pdiv , Q);
				elementLog(Pdiv,Psized);
				replaceNaN(Psized,Double.MIN_VALUE);
				elementMult(Psized,P);
				replaceNaN(Psized,Double.MIN_VALUE);
				double C = elementSum(Psized);
				System.out.println("Iteration " + iter + ": error is " + C);
			} else if(iter % 10 == 0) {
				System.out.println("Iteration " + iter);
			}

			// Stop lying about P-values
			if (iter == 100)
				divide(P , 4);
		}

		// Return solution
		return extractDoubleArray(Y);
	}
}
