package com.jujutsu.utils;

import org.jblas.DoubleMatrix;

public class BlasOps {
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
	
}
