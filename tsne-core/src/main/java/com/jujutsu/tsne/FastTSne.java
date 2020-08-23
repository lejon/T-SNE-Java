package com.jujutsu.tsne;

import static com.jujutsu.utils.EjmlOps.addRowVector;
import static com.jujutsu.utils.EjmlOps.assignAllLessThan;
import static com.jujutsu.utils.EjmlOps.assignAtIndex;
import static com.jujutsu.utils.EjmlOps.biggerThan;
import static com.jujutsu.utils.EjmlOps.colMean;
import static com.jujutsu.utils.EjmlOps.extractDoubleArray;
import static com.jujutsu.utils.EjmlOps.maximize;
import static com.jujutsu.utils.EjmlOps.replaceNaN;
import static com.jujutsu.utils.EjmlOps.setData;
import static com.jujutsu.utils.EjmlOps.setDiag;
import static com.jujutsu.utils.EjmlOps.tile;
import static com.jujutsu.utils.MatrixOps.abs;
import static com.jujutsu.utils.MatrixOps.addColumnVector;
import static com.jujutsu.utils.MatrixOps.assignValuesToRow;
import static com.jujutsu.utils.MatrixOps.concatenate;
import static com.jujutsu.utils.MatrixOps.equal;
import static com.jujutsu.utils.MatrixOps.fillMatrix;
import static com.jujutsu.utils.MatrixOps.getValuesFromRow;
import static com.jujutsu.utils.MatrixOps.mean;
import static com.jujutsu.utils.MatrixOps.negate;
import static com.jujutsu.utils.MatrixOps.range;
import static com.jujutsu.utils.MatrixOps.rnorm;
import static com.jujutsu.utils.MatrixOps.scalarInverse;
import static com.jujutsu.utils.MatrixOps.scalarMult;
import static com.jujutsu.utils.MatrixOps.sqrt;
import static com.jujutsu.utils.MatrixOps.square;
import static com.jujutsu.utils.MatrixOps.sum;
import static com.jujutsu.utils.MatrixOps.times;
import static org.ejml.dense.row.CommonOps_DDRM.add;
import static org.ejml.dense.row.CommonOps_DDRM.addEquals;
import static org.ejml.dense.row.CommonOps_DDRM.divide;
import static org.ejml.dense.row.CommonOps_DDRM.elementDiv;
import static org.ejml.dense.row.CommonOps_DDRM.elementExp;
import static org.ejml.dense.row.CommonOps_DDRM.elementLog;
import static org.ejml.dense.row.CommonOps_DDRM.elementMult;
import static org.ejml.dense.row.CommonOps_DDRM.elementPower;
import static org.ejml.dense.row.CommonOps_DDRM.elementSum;
import static org.ejml.dense.row.CommonOps_DDRM.mult;
import static org.ejml.dense.row.CommonOps_DDRM.multAddTransB;
import static org.ejml.dense.row.CommonOps_DDRM.scale;
import static org.ejml.dense.row.CommonOps_DDRM.subtract;
import static org.ejml.dense.row.CommonOps_DDRM.subtractEquals;
import static org.ejml.dense.row.CommonOps_DDRM.sumRows;
import static org.ejml.dense.row.CommonOps_DDRM.transpose;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.ejml.data.DMatrixRMaj;

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
public class FastTSne implements TSne {
	MatrixOps mo = new MatrixOps();
	protected volatile boolean abort = false;

	public static double[][] readBinaryDoubleMatrix(int rows, int columns, String fn) throws FileNotFoundException, IOException {
		File matrixFile = new File(fn);
		double [][] matrix = new double[rows][columns];
		try (DataInputStream dis =
				new DataInputStream(new BufferedInputStream(new FileInputStream(matrixFile.getAbsolutePath())))) {
			for (int i = 0; i < matrix.length; i++) {
				for (int j = 0; j < matrix[0].length; j++) {
					matrix[i][j] = dis.readDouble();
				}
			}
		}
		return matrix;
	}
	
	@Override
	public double [][] tsne(TSneConfiguration config) {
		double[][] X      = config.getXin();
		int no_dims       = config.getOutputDims();
		int initial_dims  = config.getInitialDims(); 
		double perplexity = config.getPerplexity();
		int max_iter      = config.getMaxIter();
		boolean use_pca   = config.usePca();
		
		String IMPLEMENTATION_NAME = this.getClass().getSimpleName();
		System.out.println("X:Shape is = " + X.length + " x " + X[0].length);
		System.out.println("Running " + IMPLEMENTATION_NAME + ".");
		long end = System.currentTimeMillis();
		long start = System.currentTimeMillis();
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
		DMatrixRMaj Y        = new DMatrixRMaj(rnorm(n,no_dims));
		DMatrixRMaj Ysqlmul  = new DMatrixRMaj(Y.numRows,Y.numRows);
		DMatrixRMaj dY       = new DMatrixRMaj(fillMatrix(n,no_dims,0.0));
		DMatrixRMaj iY       = new DMatrixRMaj(fillMatrix(n,no_dims,0.0));
		DMatrixRMaj gains    = new DMatrixRMaj(fillMatrix(n,no_dims,1.0));
		DMatrixRMaj btNeg    = new DMatrixRMaj(n,no_dims);
		DMatrixRMaj bt       = new DMatrixRMaj(n,no_dims);
		
		// Compute P-values
		DMatrixRMaj P        = new DMatrixRMaj(x2p(X, 1e-5, perplexity).P); // P = n x n
		DMatrixRMaj Ptr      = new DMatrixRMaj(P.numRows,P.numCols);
		DMatrixRMaj L        = new DMatrixRMaj(P); // L = n x n
		DMatrixRMaj logdivide = new DMatrixRMaj(P.numRows,P.numCols);
		DMatrixRMaj diag     = new DMatrixRMaj(fillMatrix(L.numRows, L.numCols, 0.0));
		
		transpose(P,Ptr);
		addEquals(P,Ptr);
		divide(P ,elementSum(P));
		replaceNaN(P,Double.MIN_VALUE);
		scale(4.0,P);					// early exaggeration
		maximize(P, 1e-12);
		
		System.out.println("Y:Shape is = " + Y.getNumRows() + " x " + Y.getNumCols());

		DMatrixRMaj sqed  = new DMatrixRMaj(Y.numRows,Y.numCols);
		DMatrixRMaj sum_Y = new DMatrixRMaj(1,Y.numRows);
		DMatrixRMaj num   = new DMatrixRMaj(Y.numRows, Y.numRows);
		DMatrixRMaj Q     = new DMatrixRMaj(P.numRows,P.numCols);
		
		for (int iter = 0; iter < max_iter && !abort; iter++) {
			// Compute pairwise affinities
			elementPower(Y, 2, sqed);
			sumRows(sqed, sum_Y);
			multAddTransB(-2.0, Y, Y, Ysqlmul);
			addRowVector(Ysqlmul, sum_Y);
			transpose(Ysqlmul);
			addRowVector(Ysqlmul, sum_Y);
			
			add(Ysqlmul, 1.0);
			divide(1.0,Ysqlmul);
			num.set(Ysqlmul);
			assignAtIndex(num, range(n), range(n), 0);
			divide(num , elementSum(num), Q);

			maximize(Q, 1e-12);
			
			// Compute gradient
			subtract(P, Q, L);
			elementMult(L, num);
			DMatrixRMaj rowsum = sumRows(L,null); // rowsum = nx1
			double [] rsum  = new double[rowsum.numRows];
			for (int i = 0; i < rsum.length; i++) {
				rsum[i] = rowsum.get(i,0);
			}
			setDiag(diag,rsum);
			subtract(diag, L, L);
			mult(L, Y, dY);
			scale(4.0, dY);
			
			// Perform the update
			if (iter < 20)
				momentum = initial_momentum;
			else
				momentum = final_momentum;
			
			boolean [][] boolMtrx = equal(biggerThan(dY,0.0),biggerThan(iY,0.0));
			
			
			setData(btNeg, abs(negate(boolMtrx)));
			setData(bt, abs(boolMtrx));
			
			DMatrixRMaj gainsSmall = new DMatrixRMaj(gains);
			DMatrixRMaj gainsBig   = new DMatrixRMaj(gains);
			add(gainsSmall,0.2);
			scale(0.8,gainsBig);
			
			elementMult(gainsSmall, btNeg);
			elementMult(gainsBig, bt);
			add(gainsSmall,gainsBig,gains);

			assignAllLessThan(gains, min_gain, min_gain);
			
			scale(momentum,iY);
			DMatrixRMaj gainsdY = new DMatrixRMaj(gains.numRows,dY.numCols);
			elementMult(gains , dY, gainsdY);
			scale(eta,gainsdY);
			subtractEquals(iY , gainsdY);
			addEquals(Y , iY);
			DMatrixRMaj colMeanY = colMean(Y, 0);
			DMatrixRMaj meanTile = tile(colMeanY, n, 1);
			subtractEquals(Y , meanTile);

			// Compute current value of cost function
			if (iter % 100 == 0)   {
				DMatrixRMaj Pdiv = new DMatrixRMaj(P);
				elementDiv(Pdiv , Q);
				elementLog(Pdiv,logdivide);
				replaceNaN(logdivide,Double.MIN_VALUE);
				elementMult(logdivide,P);
				replaceNaN(logdivide,Double.MIN_VALUE);
				double C = elementSum(logdivide);
				end = System.currentTimeMillis();
				System.out.printf("Iteration %d: error is %f (50 iterations in %4.2f seconds)\n", iter, C, (end - start) / 1000.0);
				if(C < 0) {
					System.err.println("Warning: Error is negative, this is usually a very bad sign!");
				}
				start = System.currentTimeMillis();
			} else if(iter % 10 == 0) {
				end = System.currentTimeMillis();
				System.out.printf("Iteration %d: (10 iterations in %4.2f seconds)\n", iter, (end - start) / 1000.0);
				start = System.currentTimeMillis();
			}

			// Stop lying about P-values
			if (iter == 100)
				divide(P , 4);
		}

		// Return solution
		return extractDoubleArray(Y);
	}
	
	public R Hbeta (double [][] D, double beta){
    	DMatrixRMaj P  = new DMatrixRMaj(D);
    	scale(-beta,P);
    	elementExp(P,P);
		double sumP = elementSum(P);   // sumP confirmed scalar
		DMatrixRMaj Dd  = new DMatrixRMaj(D);
		elementMult(Dd, P);
		double H = Math.log(sumP) + beta * elementSum(Dd) / sumP;
		scale(1/sumP,P);
		R r = new R();
		r.H = H;
		r.P = extractDoubleArray(P);
		return r;
	}

	public R x2p(double [][] X,double tol, double perplexity){
		int n               = X.length;
		double [][] sum_X   = sum(square(X), 1);
		double [][] times   = scalarMult(times(X, mo.transpose(X)), -2);
		double [][] prodSum = addColumnVector(mo.transpose(times), sum_X);
		double [][] D       = com.jujutsu.utils.MatrixOps.addRowVector(prodSum, mo.transpose(sum_X));
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

	@Override
	public void abort() {
		abort = true;
	}
}
