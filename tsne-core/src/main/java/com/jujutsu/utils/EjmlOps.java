package com.jujutsu.utils;

import static org.ejml.ops.CommonOps.divide;
import static org.ejml.ops.CommonOps.sumCols;

import org.ejml.data.DenseMatrix64F;

public class EjmlOps {

	public static void maximize(DenseMatrix64F p, double minval) {
		int rows = p.getNumRows();
		int cols = p.getNumCols();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				double val = p.get(i, j);
				if(val<minval) p.unsafe_set(j, j, minval);
			}
		}
	}

	/**
	 * Returns a new matrix of booleans where true is set if the value in the matrix is
	 * bigger than value
	 * @param fastTSne TODO
	 * @param matrix
	 * @param value
	 * @return new matrix with booelans with values matrix1[i,j] == matrix2[i,j]
	 */
	public static boolean [][] biggerThan(DenseMatrix64F matrix, double value) {
		boolean [][] equals = new boolean[matrix.numRows][matrix.numCols];
		for (int i = 0; i < matrix.numRows; i++) {
			for (int j = 0; j < matrix.numCols; j++) {
				equals[i][j] = Double.compare(matrix.get(i,j), value) == 1;
			}
		}
		return equals;
	}

	/**
	 * Sets the diagonal of 'diag' to the values of 'diagElements' as long 
	 * as possible (i.e while there are elements left in diag and the dim of 'diag'
	 * is big enough...
	 * Note: This method ONLY affect the diagonal elements the others are left as
	 * when passed in.
	 * @param fastTSne TODO
	 * @param diag Modified to contain the elements of 'diagElements' on its diagonal
	 * @param diagElems
	 */
	public static void setDiag(DenseMatrix64F diag, double[] diagElems) {
		int idx = 0; 
		while(idx<diag.numCols&&idx<diag.numRows&&idx<diagElems.length) {
			diag.set(idx, idx, diagElems[idx++]);
		}
	}

	/**
	 * <p>
	 * Sets the data of<code>target</code> to that of the input matrix with the values and shape defined by the 2D array 'data'.
	 * It is assumed that 'data' has a row-major formatting:<br>
	 *  <br>
	 * data[ row ][ column ]
	 * </p>
	 * @param fastTSne TODO
	 * @param target 2D DenseMatrix. Modified to contain the values in 'data'.
	 * @param data 2D array representation of the matrix. Not modified.
	 */
	public static void setData(DenseMatrix64F target, double[][] data) {
	    int numRows = data.length;
	    int numCols = data[0].length;
	
	    double [] targetData = new double[ numRows*numCols ];
	
	    int pos = 0;
	    for( int i = 0; i < numRows; i++ ) {
	        double []row = data[i];
	
	        if( row.length != numCols ) {
	            throw new IllegalArgumentException("All rows must have the same length");
	        }
	
	        System.arraycopy(row,0,targetData,pos,numCols);
	        pos += numCols;
	    }
	    
	    target.setData(targetData);
	}

	/** 
	 * Replaces NaN's with repl
	 * @param matrix
	 * @param repl
	 * @return
	 */
	public static void replaceNaN(DenseMatrix64F matrix, double repl) {
		for (int i = 0; i < matrix.numRows; i++) {
			for (int j = 0; j < matrix.numCols; j++) {
				if(Double.isNaN(matrix.get(i,j))) {
					matrix.set(i,j,repl);
				} 
			}
		}
	}

	public static DenseMatrix64F fillWithRow(DenseMatrix64F matrix, int setrow) {
		int rows = matrix.numRows;
		int cols = matrix.numCols;
		DenseMatrix64F result = new DenseMatrix64F(rows,cols);
		for (int row = 0; row < rows; row++) {
			for (int col = 0; col < cols; col++) {
				result.set(row,col, matrix.get(setrow,col));				
			}
		}
		return result;
	}

	public static DenseMatrix64F tile(DenseMatrix64F matrix, int rowtimes, int coltimes) {
		DenseMatrix64F result = new DenseMatrix64F(matrix.numRows*rowtimes,matrix.numCols*coltimes);
		for (int i = 0, resultrow = 0; i < rowtimes; i++) {
			for (int j = 0; j < matrix.numRows; j++) {
				for (int k = 0, resultcol = 0; k < coltimes; k++) {
					for (int l = 0; l < matrix.numCols; l++) {
						result.set(resultrow,resultcol++,matrix.get(j,l));
					}
				}
				resultrow++;
			}
		}
		return result;
	}

	/**
	 * All values in matrix that is less than <code>lessthan</code> is assigned
	 * the value <code>assign</code>
	 * @param matrix
	 * @param lessthan
	 * @param assign
	 * @return
	 */
	public static void assignAllLessThan(DenseMatrix64F matrix, double lessthan, double assign) {
		for (int i = 0; i < matrix.numRows; i++) {
			for (int j = 0; j < matrix.numCols; j++) {
				if( matrix.get(i,j) < lessthan) {
					matrix.set(i,j,assign);
				}
			}
		}
	}

	public static DenseMatrix64F colMean(DenseMatrix64F y, int i) {
		DenseMatrix64F colmean = new DenseMatrix64F(1,y.numCols);
		sumCols(y,colmean);
		divide(colmean, y.numRows);
		return colmean;
	}

	public static void addRowVector(DenseMatrix64F matrix, DenseMatrix64F rowvector) {
		for (int i = 0; i < matrix.numRows; i++) {
			for (int j = 0; j < matrix.numCols; j++) {
				matrix.set(i,j,matrix.get(i,j) + rowvector.get(0,j));
			}
		}
	}

	public static void assignAtIndex(DenseMatrix64F num, int[] range, int[] range1, double value) {
		for (int j = 0; j < range.length; j++) {
			num.set(range[j], range1[j], value);
		}
	}

	public static double [][] extractDoubleArray(DenseMatrix64F p) {
		int rows = p.getNumRows();
		int cols = p.getNumCols();
		double [][] result = new double[rows][cols];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				result[i][j] = p.get(i, j);
			}
		}
		return result;
	}

}
