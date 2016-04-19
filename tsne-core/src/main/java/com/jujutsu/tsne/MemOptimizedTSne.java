package com.jujutsu.tsne;

import static com.jujutsu.utils.EjmlOps.addRowVector;
import static com.jujutsu.utils.EjmlOps.assignAllLessThan;
import static com.jujutsu.utils.EjmlOps.assignAtIndex;
import static com.jujutsu.utils.EjmlOps.biggerThan;
import static com.jujutsu.utils.EjmlOps.colMean;
import static com.jujutsu.utils.EjmlOps.maximize;
import static com.jujutsu.utils.EjmlOps.replaceNaN;
import static com.jujutsu.utils.EjmlOps.setData;
import static com.jujutsu.utils.EjmlOps.setDiag;
import static com.jujutsu.utils.EjmlOps.tile;
import static com.jujutsu.utils.MatrixOps.abs;
import static com.jujutsu.utils.MatrixOps.equal;
import static com.jujutsu.utils.MatrixOps.fillMatrix;
import static com.jujutsu.utils.MatrixOps.negate;
import static com.jujutsu.utils.MatrixOps.range;
import static com.jujutsu.utils.MatrixOps.rnorm;
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
public class MemOptimizedTSne extends FastTSne {
	
	public double [][] tsne(double[][] X, int no_dims, int initial_dims, double perplexity, int max_iter, boolean use_pca) {
		String IMPLEMENTATION_NAME = this.getClass().getSimpleName();
		System.out.println("X:Shape is = " + X.length + " x " + X[0].length);
		System.out.println("Running " + IMPLEMENTATION_NAME + ".");
		// Initialize variables
		if(use_pca && X[0].length > initial_dims && initial_dims > 0) {
			PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
			X = pca.pca(X, initial_dims);
			System.out.println("X:Shape after PCA is = " + X.length + " x " + X[0].length);
			System.out.println(MatrixOps.doubleArrayToPrintString(X,10,10));
		}
		int n = X.length;
		double momentum = .5;
		double initial_momentum = 0.5;
		double final_momentum   = 0.8;
		int eta                 = 500;
		double min_gain         = 0.01;
		DenseMatrix64F Y        = new DenseMatrix64F(rnorm(n,no_dims));
		DenseMatrix64F Ysqlmul  = new DenseMatrix64F(Y.numRows,Y.numRows); // Ysqlmul = n x n
		DenseMatrix64F dY       = new DenseMatrix64F(fillMatrix(n,no_dims,0.0));
		DenseMatrix64F iY       = new DenseMatrix64F(fillMatrix(n,no_dims,0.0));
		DenseMatrix64F gains    = new DenseMatrix64F(fillMatrix(n,no_dims,1.0));
		DenseMatrix64F btNeg    = new DenseMatrix64F(n,no_dims);
		DenseMatrix64F bt       = new DenseMatrix64F(n,no_dims);
		
		// Compute P-values
		DenseMatrix64F P        = new DenseMatrix64F(x2p(X, 1e-5, perplexity).P); // P = n x n
		DenseMatrix64F Psized   = new DenseMatrix64F(P.numRows,P.numCols);        // L = n x n
		DenseMatrix64F diag     = new DenseMatrix64F(fillMatrix(Psized.numRows,Psized.numCols,0.0));
		
		transpose(P,Psized);
		addEquals(P,Psized);
		divide(P ,elementSum(P));
		replaceNaN(P,Double.MIN_VALUE);
		scale(4.0,P);					// early exaggeration
		maximize(P, 1e-12);
		
		System.out.println("Using perplexity: " + perplexity);
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
			assignAtIndex(Ysqlmul, range(n), range(n), 0);
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
			
			boolean [][] boolMtrx = equal(biggerThan(dY,0.0),biggerThan(iY,0.0));
			
			
			setData(btNeg, abs(negate(boolMtrx)));
			setData(bt, abs(boolMtrx));
			
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
				if(C < 0) {
					System.err.println("Warning: Error is negative, this is usually a very bad sign!");
				}
			} else if(iter % 10 == 0) {
				System.out.println("Iteration " + iter);
			}

			// Stop lying about P-values
			if (iter == 100)
				divide(P , 4);
		}

		// Return solution
		return MatrixOps.extractDoubleArray(Y);
	}
}
