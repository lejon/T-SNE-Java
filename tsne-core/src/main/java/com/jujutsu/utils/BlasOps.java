package com.jujutsu.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.ejml.data.DMatrixRMaj;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.netlib.util.intW;

import com.github.fommil.netlib.BLAS;
import com.github.fommil.netlib.LAPACK;

public class BlasOps {
	
	void benchmark(int size) {
		int m = (int) Math.sqrt(size);

		// random matrices are full rank (and can always be inverted if square)
		// http://www.sciencedirect.com/science/article/pii/S0096300306009040
		double[] a = new double[m * m];
		//for (int i = 0; i < a.length; i++) {
		//	a[i] = ThreadLocalRandom.current().nextDouble();
		//}
		a[0] = 2;
		a[4] = 2;
		a[8] = 2;
		DoubleMatrix pInv = new DoubleMatrix(a);

		System.out.println("Start is: " + pInv.reshape(m, m));

		double[] aOrig = Arrays.copyOf(a, a.length);
		double[] b = new double[1];
		int[] p = new int[m];
		intW info = new intW(0);

		LAPACK.getInstance().dgetri(m, a, m, p, b, -1, info);
		//log.info(m + " supposedly has optimal work of " + b[0]);
		b = new double[(int)b[0]];

		LAPACK.getInstance().dgetrf(m, m, a, m, p, info);
		if (info.val != 0)
			throw new IllegalArgumentException();
		LAPACK.getInstance().dgetri(m, a, m, p, b, b.length, info);
		if (info.val != 0)
			throw new IllegalArgumentException();

		// quick check
		double[] c = new double[m * m];
		BLAS.getInstance().dgemm("N", "N", m, m, m, 1, aOrig, m, a, m, 0, c, m);
		pInv = new DoubleMatrix(a);

		System.out.println("Result is: " + pInv.reshape(m, m));

		System.out.println("BlasInvert Result is: " + blasInvert(new DoubleMatrix(a)));
	}

	public static DoubleMatrix blasInvert(DoubleMatrix myPrecision) {
		double[] b = new double[myPrecision.columns];
		int[] p = new int[myPrecision.columns];
		intW info = new intW(0);
		int m = myPrecision.columns;
		double [] a = myPrecision.toArray();
		//LAPACK.getInstance().dgetri(m, a, m, p, b, -1, info);
		//log.info(m + " supposedly has optimal work of " + b[0]);
		//b = new double[(int)b[0]];
		LAPACK.getInstance().dgetrf(m, m, a, m, p, info);
		if (info.val != 0)
			throw new IllegalArgumentException();
		LAPACK.getInstance().dgetri(m, a, m, p, b, b.length, info);
		if (info.val != 0)
			throw new IllegalArgumentException();
		DoubleMatrix pInv = new DoubleMatrix(a);
		return pInv.reshape(myPrecision.rows, myPrecision.columns);
	}

	public static DMatrixRMaj blasInvertDense(DMatrixRMaj myPrecision) {
		DoubleMatrix mp = new DoubleMatrix(myPrecision.getData().clone());
		mp = mp.reshape(myPrecision.numRows, myPrecision.numCols);
		DoubleMatrix inv = blasInvert(mp);
		DMatrixRMaj res = new DMatrixRMaj(inv.toArray2());
		return res;
	}

	
	public static DoubleMatrix square(DoubleMatrix in) {
		DoubleMatrix res = in.dup();
		for (int i = 0; i < res.getLength(); i++) {
			double val = res.get(i);
			res.put(i,val*val);
		}
		return res;
	}

	public static DoubleMatrix scalarInverse(DoubleMatrix in) {
		DoubleMatrix res = in.dup();
		for (int i = 0; i < res.getLength(); i++) {
			double val = res.get(i);
			res.put(i,1/val);
		}
		return res;
	}

	public static void assignAtIndex(DoubleMatrix num, int[] range, int[] range1, double value) {
		for (int j = 0; j < range.length; j++) {
			num.put(range[j],range1[j],value);
		}
	}

	/**
	 * Returns a new matrix of booleans where true is set if the values to the two matrices are
	 * the same at that index
	 * @param matrix1
	 * @param matrix2
	 * @return new matrix with booelans with values matrix1[i,j] == matrix2[i,j]
	 */
	public static boolean [][] equal(DoubleMatrix matrix1, DoubleMatrix matrix2) {
		boolean [][] equals = new boolean[matrix1.rows][matrix1.columns];
		if( matrix1.length != matrix2.length) {
			throw new IllegalArgumentException("Dimensions does not match");
		}
		if( matrix1.columns != matrix2.columns) {
			throw new IllegalArgumentException("Dimensions does not match");
		}
		for (int i = 0; i < matrix1.rows; i++) {
			for (int j = 0; j < matrix1.columns; j++) {
				equals[i][j] = Double.compare(matrix1.get(i,j), matrix2.get(i,j)) == 0;
			}
		}
		return equals;
	}
	/**
	 * All values in matrix that is less than <code>lessthan</code> is assigned
	 * the value <code>assign</code>
	 * @param matrix
	 * @param lessthan
	 * @param assign
	 * @return
	 */
	public static void assignAllLessThan(DoubleMatrix matrix, double lessthan, double assign) {
		for (int i = 0; i < matrix.length; i++) {
			if( matrix.get(i) < lessthan) {
				matrix.put(i, assign);
			}
		}
	}

	/**
	 * Returns a new matrix with values that are the log of the input matrix
	 * @param matrix
	 * @return  same matrix with values log'ed
	 */
	public static DoubleMatrix log(DoubleMatrix m1) {
		DoubleMatrix matrix = m1.dup();
		for (int i = 0; i < matrix.length; i++) {
			matrix.put(i, Math.log(m1.get(i)));
		}
		return matrix;
	}

	/**
	 * Returns a new matrix with values that are the log of the input matrix
	 * @param matrix
	 * @param infAsZero treat +- Infinity as zero, i.e replaces Infinity with 0.0
	 * if set to true
	 * @return  same matrix with values log'ed
	 */
	public static DoubleMatrix log(DoubleMatrix m1, boolean infAsZero) {
		DoubleMatrix matrix = m1.dup();
		for (int i = 0; i < matrix.length; i++) {
			matrix.put(i, Math.log(m1.get(i)));
			if(infAsZero && Double.isInfinite(matrix.get(i)))
				matrix.put(i,0.0);
		}
		return matrix;
	}

	/** 
	 * Replaces NaN's with repl
	 * @param matrix
	 * @param repl
	 * @return
	 */
	public static DoubleMatrix replaceNaN(DoubleMatrix matrix, double repl) {
		DoubleMatrix result = matrix.dup();
		for (int i = 0; i < matrix.length; i++) {
			if(Double.isNaN(matrix.get(i))) {
				result.put(i,repl);
			} else {
				result.put(i,matrix.get(i));
			}
		}

		return result;
	}
	
	public static DoubleMatrix tile(DoubleMatrix matrix, int rowtimes, int coltimes) {
		DoubleMatrix result = new DoubleMatrix(matrix.rows*rowtimes,matrix.columns*coltimes);
		for (int i = 0, resultrow = 0; i < rowtimes; i++) {
			for (int j = 0; j < matrix.rows; j++) {
				for (int k = 0, resultcol = 0; k < coltimes; k++) {
					for (int l = 0; l < matrix.columns; l++) {
						result.put(resultrow,resultcol++,matrix.get(j,l));
					}
				}
				resultrow++;
			}
		}

		return result;
	}

	/**
	 * Returns a new matrix of booleans where true is set if the value in the matrix is
	 * bigger than value
	 * @param matrix
	 * @param value
	 * @return new matrix with booelans with values matrix1[i,j] == matrix2[i,j]
	 */
	public static boolean [][] biggerThan(DoubleMatrix matrix, double value) {
		boolean [][] equals = new boolean[matrix.rows][matrix.columns];
		for (int i = 0; i < matrix.rows; i++) {
			for (int j = 0; j < matrix.columns; j++) {
				equals[i][j] = Double.compare(matrix.get(i,j), value) == 1;
			}
		}
		return equals;
	}

	/**
	 * @param booleans
	 * @return
	 */
	public static DoubleMatrix abs(boolean [][] booleans) {
		DoubleMatrix absolutes = new DoubleMatrix(booleans.length,booleans[0].length);
		for (int i = 0; i < booleans.length; i++) {
			for (int j = 0; j < booleans[0].length; j++) {
				absolutes.put(i,j, booleans[i][j] ? 1 : 0);
			}
		}
		return absolutes;
	}

	// Unit Tested
	/**
	 * Returns a new matrix with values exponentiated
	 * @param matrix
	 * @return new matrix with values exponentiated
	 */
	public static DoubleMatrix exp(DoubleMatrix m1) {
		DoubleMatrix matrix = m1.dup();
		for (int i = 0; i < matrix.length; i++) {
			matrix.put(i, Math.exp(m1.get(i)));
		}
		return matrix;
	}
	
	public static boolean containsNaNs(DoubleMatrix matrix) {
		for (int i = 0; i < matrix.length; i++) {
			if(Double.isNaN(matrix.get(i))) {
				return true;
			}
		}
		return false;	
	}
	
	public static DoubleMatrix sign(DoubleMatrix matrix) {
		DoubleMatrix signs = new DoubleMatrix(matrix.rows,matrix.columns);
		for (int i = 0; i < matrix.length; i++) {
				signs.put(i, matrix.get(i) >= 0 ? 1 : -1);
		}
		return signs;
	}
	
	/**
	 * Returns new vetcor with vals sqrt'ed
	 * @param vector
	 * @return new vector with values sqrt'ed
	 */
	public static DoubleMatrix sqrt(DoubleMatrix m1) {
		DoubleMatrix matrix = m1.dup();
		for (int i = 0; i < matrix.length; i++) {
			matrix.put(i,Math.sqrt(m1.get(i)));
		}
		return matrix;
	}
	
	// Adapted from: http://www.programcreek.com/java-api-examples/index.php?source_dir=darks-learning-master/src/main/java/darks/learning/dimreduce/pca/PCA.java
	// with License: 
	
	/**
	 *  
	 * Copyright 2014 The Darks Learning Project (Liu lihua) 
	 * 
	 * Licensed under the Apache License, Version 2.0 (the "License"); 
	 * you may not use this file except in compliance with the License. 
	 * You may obtain a copy of the License at 
	 * 
	 *     http://www.apache.org/licenses/LICENSE-2.0 
	 * 
	 * Unless required by applicable law or agreed to in writing, software 
	 * distributed under the License is distributed on an "AS IS" BASIS, 
	 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
	 * See the License for the specific language governing permissions and 
	 * limitations under the License. 
	 */
	
	/**
	 * Reduce matrix dimension
	 * 
	 * @param source Source matrix
	 * @param dimension Target dimension
	 * @return Target matrix
	 */
	public static DoubleMatrix pca(DoubleMatrix source, int dimension)
	{
		//C=X*X^t / m
		DoubleMatrix covMatrix = source.mmul(source.transpose()).div(source.columns);
		ComplexDoubleMatrix eigVal = Eigen.eigenvalues(covMatrix);
		ComplexDoubleMatrix[] eigVectorsVal = Eigen.eigenvectors(covMatrix);
		ComplexDoubleMatrix eigVectors = eigVectorsVal[0];
		//Sort sigen vector from big to small by eigen values 
		List<PCABean> beans = new ArrayList<PCABean>();
		for (int i = 0; i < eigVectors.columns; i++)
		{
			beans.add(new PCABean(eigVal.get(i).real(), eigVectors.getColumn(i)));
		}
		Collections.sort(beans);
		DoubleMatrix newVec = new DoubleMatrix(dimension, beans.get(0).vector.rows);
		for (int i = 0; i < dimension; i++)
		{
			ComplexDoubleMatrix dm = beans.get(i).vector;
			DoubleMatrix real = dm.getReal();
			newVec.putRow(i, real);
		}
		return newVec.mmul(source);
	}
	
	static class PCABean implements Comparable<PCABean>
	{
		double eigenValue;
		
		ComplexDoubleMatrix vector;
		
		public PCABean(double eigenValue, ComplexDoubleMatrix vector)
		{
			super();
			this.eigenValue = eigenValue;
			this.vector = vector;
		}



		@Override
		public int compareTo(PCABean o)
		{
			return Double.compare(o.eigenValue, eigenValue);
		}



		@Override
		public String toString()
		{
			return "PCABean [eigenValue=" + eigenValue + ", vector=" + vector + "]";
		}
		
	}
	
}
