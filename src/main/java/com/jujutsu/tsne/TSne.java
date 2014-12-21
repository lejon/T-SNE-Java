package com.jujutsu.tsne;

import java.util.Random;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.TimeUnit;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import Jama.QRDecomposition;

/**
 *
 * User: Leif Jonsson (leif.jonsson@ericsson.com)
 * 
 * This is a port of van der Maaten and Hintons Python implementation of t-sne
 *
 */
public class TSne {
	Random rnd = new Random();
	static ExecutorService	transposerPool = Executors.newFixedThreadPool(4);
	static BlockingQueue<Object> transposesQueue = new LinkedBlockingQueue<Object>();
	private static ForkJoinPool pool = new ForkJoinPool();

	/**
	 * Generate random draw from Normal with mean mu and std. dev sigma
	 * @param mu
	 * @param sigma
	 * @return random sample
	 */
	static double dnrom(double mu, double sigma) {
		return mu + (ThreadLocalRandom.current().nextGaussian() * sigma);
	}

	/**
	 * Returns a new matrix which is the transpose of input matrix
	 * @param matrix
	 * @return
	 */
	static double[][] origtranspose(double[][] matrix) {
		int cols = matrix[0].length;
		int rows = matrix.length;
		double[][] transpose = new double[cols][rows];
		for (int i = 0; i < cols; i++)
			for (int j = 0; j < rows; j++)
				transpose[i][j] = matrix[j][i];
		return transpose;
	}

	static double[][] transpose(double[][] matrix) {
		return transpose(matrix, 1000);
	}
	/**
	 * Returns a new matrix which is the transpose of input matrix
	 * @param matrix
	 * @return
	 */
	static double[][] transpose(double[][] matrix, int ll) {
		int cols = matrix[0].length;
		int rows = matrix.length;
		double[][] transpose = new double[cols][rows];
		if(rows < 100 ) {
			for (int i = 0; i < cols; i++)
				for (int j = 0; j < rows; j++)
					transpose[i][j] = matrix[j][i];
		} else {
			MatrixTransposer process = new MatrixTransposer(matrix, transpose,0,rows,ll);                
			pool.invoke(process);
//			try {
//				pool.awaitTermination(50L, TimeUnit.MILLISECONDS);
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//			}
//			//pool.shutdown();
		}
		return transpose;
	}

	static class MatrixTransposer extends RecursiveAction {
		static final long serialVersionUID = 1L;
		double [][] orig;
		double [][] transpose;
		int startRow = -1;
		int endRow = -1;
		int limit = 1000;

		public MatrixTransposer(double [][] orig, double [][] transpose, int startRow, int endRow, int ll) {
			this.limit = ll;
			this.orig = orig;
			this.transpose = transpose;
			this.startRow = startRow;
			this.endRow = endRow;
		}
		
		public MatrixTransposer(double [][] orig, double [][] transpose, int startRow, int endRow) {
			this.orig = orig;
			this.transpose = transpose;
			this.startRow = startRow;
			this.endRow = endRow;
		}

		@Override
		protected void compute() {
			try {
				if ( (endRow-startRow) <= limit ) {
					int cols = orig[0].length;
					for (int i = 0; i < cols; i++) {
						for (int j = startRow; j < endRow; j++) {
							transpose[i][j] = orig[j][i];
						}
					}
				}
				else {
					int range = (endRow-startRow);
					int startRow1 = startRow;
					int endRow1 = startRow + (range / 2);
					int startRow2 = endRow1;
					int endRow2 = endRow;
					invokeAll(new MatrixTransposer(orig, transpose, startRow1, endRow1, limit),
							new MatrixTransposer(orig, transpose, startRow2, endRow2, limit));
				}
			}
			catch ( Exception e ) {
				e.printStackTrace();
			}
		}
	}

	static void plainTranspose(double[][] matrix, double[][] transpose, int startRow, int endRow) {
		int cols = matrix[0].length;
		for (int i = 0; i < cols; i++)
			for (int j = startRow; j < endRow; j++)
				transpose[i][j] = matrix[j][i];
	}

	/**
	 * Destructively sets the values in matrix to its exponentiated value
	 * @param matrix
	 * @return same matrix with values exponentiated
	 */
	static double [][] exp(double [][] matrix) {
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				matrix[i][j] = Math.exp(matrix[i][j]);
			}
		}
		return matrix;
	}

	/**
	 * Destructively sets the values in vector to its square root
	 * @param vector
	 * @return same vector with values sqrt'ed
	 */
	static double [] sqrt(double [] vector) {
		for (int i = 0; i < vector.length; i++) {
			vector[i] = Math.sqrt(vector[i]);
		}
		return vector;
	}

	/**
	 * @param vector
	 * @return mean of values in vector
	 */
	static double mean(double [] vector) {
		double sum = 0.0;
		for (int i = 0; i < vector.length; i++) {
			sum +=vector[i];
		}
		return sum/vector.length;
	}

	/**
	 * Destructively sets the values in matrix to its log value
	 * @param matrix
	 * @return  same matrix with values log'ed
	 */
	static double [][] log(double [][] matrix) {
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				matrix[i][j] = Math.log(matrix[i][j]);
			}
		}
		return matrix;
	}

	/**
	 * @param matrix
	 * @return scalar inverse of matrix
	 */
	static double [][] scalarInverse(double [][] matrix) {
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				matrix[i][j] = 1/matrix[i][j];
			}
		}
		return matrix;
	}

	/**
	 * @param vector
	 * @return scalar inverse of vector
	 */
	static double [] scalarInverse(double [] vector) {
		for (int i = 0; i < vector.length; i++) {
			vector[i] = 1/vector[i];
		}
		return vector;
	}

	/**
	 * @param m
	 * @param n
	 * @return new 2D matrix with normal random values with mean 0 and std. dev 1
	 */
	static double[][] rnorm(int m, int n) {
		double[][] array = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < array[i].length; j++) {				
				array[i][j] = dnrom(0.0,1.0);
			}
		}
		return array;
	}

	/**
	 * Returns a new matrix of booleans where true is set if the values to the two matrices are
	 * the same at that index
	 * @param matrix1
	 * @param matrix2
	 * @return new matrix with booelans with values matrix1[i,j] == matrix2[i,j]
	 */
	static boolean [][] equal(double [][] matrix1, double [][] matrix2) {
		boolean [][] equals = new boolean[matrix1.length][matrix1[0].length];
		if( matrix1.length != matrix2.length) {
			throw new IllegalArgumentException("Dimensions does not match");
		}
		if( matrix1[0].length != matrix2[0].length) {
			throw new IllegalArgumentException("Dimensions does not match");
		}
		for (int i = 0; i < matrix1.length; i++) {
			for (int j = 0; j < matrix1[0].length; j++) {
				equals[i][j] = Double.compare(matrix1[i][j], matrix2[i][j]) == 0;
			}
		}
		return equals;
	}

	/**
	 * @param booleans
	 * @return new matrix with booleans which are the negations of the input
	 */
	static boolean [][] negate(boolean [][] booleans) {
		boolean [][] negates = new boolean[booleans.length][booleans[0].length];
		for (int i = 0; i < booleans.length; i++) {
			for (int j = 0; j < booleans[0].length; j++) {
				negates[i][j] = !booleans[i][j];
			}
		}
		return negates;
	}

	/**
	 * @param booleans
	 * @return
	 */
	static double [][] abs(boolean [][] booleans) {
		double [][] absolutes = new double[booleans.length][booleans[0].length];
		for (int i = 0; i < booleans.length; i++) {
			for (int j = 0; j < booleans[0].length; j++) {
				absolutes[i][j] = booleans[i][j] ? 1 : 0;
			}
		}
		return absolutes;
	}

	static double [][] sign(double [][] matrix) {
		double [][] signs = new double[matrix.length][matrix[0].length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				signs[i][j] = matrix[i][j] > 0 ? 1 : (matrix[i][j] < 0 ? -1 : 0);
			}
		}
		return signs;
	}

	static double [][] mean(double [][] matrix, int axis) {
		// Axis = 0 => sum columns
		// Axis = 1 => sum rows
		// Axis = 3 => global (returns a 1 element array with the result)
		double [][] result;
		if( axis == 0) {
			result = new double[1][matrix[0].length];
			for (int j = 0; j < matrix[0].length; j++) {
				double colsum = 0.0;
				for (int i = 0; i < matrix.length; i++) {
					colsum += matrix[i][j];
				}
				result[0][j] = colsum / matrix.length;
			}
		}   else if (axis == 1) {
			result = new double[matrix.length][1];
			for (int i = 0; i < matrix.length; i++) {
				double rowsum = 0.0;
				for (int j = 0; j < matrix[0].length; j++) {
					rowsum += matrix[i][j];
				}
				result[i][0] = rowsum / matrix[0].length;
			}
		}   else if (axis == 2) {
			result = new double[1][1];
			for (int j = 0; j < matrix[0].length; j++) {
				for (int i = 0; i < matrix.length; i++) {
					result[0][0] += matrix[i][j];
				}
			}
			result[0][0] /=  (matrix[0].length *  matrix.length);
		}else {
			throw  new IllegalArgumentException("Axes other than 0,1,2 is unsupported");
		}
		return result;
	}

	// Should be called dim-sum! :)
	static double [][] sum(double [][] matrix, int axis) {
		// Axis = 0 => sum columns
		// Axis = 1 => sum rows
		double [][] result;
		if( axis == 0) {
			result = new double[1][matrix[0].length];
			for (int j = 0; j < matrix[0].length; j++) {
				double rowsum = 0.0;
				for (int i = 0; i < matrix.length; i++) {
					rowsum += matrix[i][j];
				}
				result[0][j] = rowsum;
			}
		}   else if (axis == 1) {
			result = new double[matrix.length][1];
			for (int i = 0; i < matrix.length; i++) {
				double colsum = 0.0;
				for (int j = 0; j < matrix[0].length; j++) {
					colsum += matrix[i][j];
				}
				result[i][0] = colsum;
			}
		}   else {
			throw  new IllegalArgumentException("Axes other than 0,1 is unsupported");
		}
		return result;
	}


	/**
	 * @param matrix
	 * @return sum of all values in the matrix
	 */
	static double sum(double [][] matrix) {
		double sum = 0.0;
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				sum+=matrix[i][j];
			}
		}
		return sum;
	}

	/**
	 * Return a new matrix with the max value of either the value in the matrix 
	 * or maxval otherwise 
	 * @param matrix
	 * @param maxval
	 * @return
	 */
	static double [][] maximum(double [][] matrix, double maxval) {
		double [][] maxed = new double[matrix.length][matrix[0].length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				maxed[i][j] = matrix[i][j] > maxval ? matrix[i][j] : maxval;
			}
		}
		return maxed;
	}

	/**
	 * All values in matrix that is less than <code>lessthan</code> is assigned
	 * the value <code>assign</code>
	 * @param matrix
	 * @param lessthan
	 * @param assign
	 * @return
	 */
	static double [][] assignAllLessThan(double[][] matrix, double lessthan, double assign) {
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				if( matrix[i][j] < lessthan) {
					matrix[i][j] = assign;
				}
			}
		}
		return matrix;
	}

	/**
	 * @param matrix
	 * @return a new matrix with the values of matrix squared
	 */
	static double [][] square(double [][] matrix) {
		return scalarPow(matrix,2);
	}

	/** 
	 * Replaces NaN's with repl
	 * @param matrix
	 * @param repl
	 * @return
	 */
	static double [][] replaceNaN(double [][] matrix, double repl) {
		double [][] result = new double[matrix.length][matrix[0].length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				if(Double.isNaN(matrix[i][j])) {
					result[i][j] = repl;
				} else {
					result[i][j] = matrix[i][j];
				}
			}
		}
		return result;
	}

	static double [][] scalarPow(double [][] matrix, double power) {
		double [][] result = new double[matrix.length][matrix[0].length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				result[i][j] += Math.pow(matrix[i][j],power);
			}
		}
		return result;
	}

	static double [][] addColumnVector(double [][] matrix, double [][] colvector) {
		double [][] result = new double[matrix.length][matrix[0].length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				result[i][j] = matrix[i][j] + colvector[i][0];
			}
		}
		return result;
	}

	static double [][] addRowVector(double [][] matrix, double [][] rowvector) {
		double [][] result = new double[matrix.length][matrix[0].length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				result[i][j] = matrix[i][j] + rowvector[0][j];
			}
		}
		return result;
	}

	static double [][] fillWithRowOld(double [][] matrix, int row) {
		double [][] result = new double[matrix.length][matrix[0].length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				result[i][j] = matrix[row][j];
			}
		}
		return result;
	}

	static double [][] fillWithRow(double [][] matrix, int row) {
		int rows = matrix.length;
		int cols = matrix[0].length;
		double [][] result = new double[rows][cols];
		for (int i = 0; i < rows; i++) {
			System.arraycopy(matrix[row], 0, result[i], 0, cols);
		}
		return result;
	}

	static double [][] tile(double [][] matrix, int rowtimes, int coltimes) {
		double [][] result = new double[matrix.length*rowtimes][matrix[0].length*coltimes];
		for (int i = 0, resultrow = 0; i < rowtimes; i++) {
			for (int j = 0; j < matrix.length; j++) {
				for (int k = 0, resultcol = 0; k < coltimes; k++) {
					for (int l = 0; l < matrix[0].length; l++) {
						result[resultrow][resultcol++] = matrix[j][l];
					}
				}
				resultrow++;
			}
		}

		return result;
	}

	static double[][] normalize(double[][] x, double[] meanX, double[] stdevX) {
		double[][] y = new double[x.length][x[0].length];
		for (int i = 0; i < y.length; i++)
			for (int j = 0; j < y[i].length; j++)
				y[i][j] = (x[i][j] - meanX[j]) / stdevX[j];
		return y;
	}

	static int [] range(int n) {
		int [] result = new int[n];
		for (int i = 0; i < n; i++) {
			result[i] = i;
		}
		return result;
	}

	static int [] range(int a, int b) {
		if( b < a ) {
			throw new IllegalArgumentException("b has to be larger than a");
		}
		int val = a;
		int [] result = new int[b-a];
		for (int i = 0; i < (b-a); i++) {
			result[i] = val++;
		}
		return result;
	}

	static int [] concatenate(int [] v1,int [] v2) {
		int [] result = new int[v1.length+v2.length];
		int index = 0;
		for (int i = 0; i < v1.length; i++, index++) {
			result[index] = v1[index];
		}
		for (int i = 0; i < v2.length; i++, index++) {
			result[index] = v2[i];
		}
		return result;
	}

	static double [][] scalarMultiply(double [][] v1,double [][] v2) {
		if( v1.length != v2.length || v1[0].length != v2[0].length ) {
			throw new IllegalArgumentException("a and b has to be of equal dimensions");
		}
		double [][] result = new double[v1.length][v1[0].length];
		for (int i = 0; i < v1.length; i++) {
			for (int j = 0; j < v1[0].length; j++) {
				result[i][j] = v1[i][j] * v2[i][j];
			}
		}
		return result;
	}

	static void assignAtIndex(double[][] num, int[] range, int[] range1, double value) {
		for (int j = 0; j < range.length; j++) {
			num[range[j]][range1[j]] = value;
		}
	}

	static double [][] getValuesFromRow(double[][] matrix, int row, int[] indicies) {
		double [][] values = new double[1][indicies.length];
		for (int j = 0; j < indicies.length; j++) {
			values[0][j] = matrix[row][indicies[j]];
		}
		return values;
	}

	static void assignValuesToRow(double[][] matrix, int row, int[] indicies, double [] values) {
		if( indicies.length != values.length ) {
			throw new IllegalArgumentException("Length of indicies and values have to be equal");
		}
		for (int j = 0; j < indicies.length; j++) {
			matrix[row][indicies[j]] = values[j];
		}
	}

	static double[] stddev(double[][] v) {
		double[] var = variance(v);
		for (int i = 0; i < var.length; i++)
			var[i] = Math.sqrt(var[i]);
		return var;
	}

	static double[] variance(double[][] v) {
		int m = v.length;
		int n = v[0].length;
		double[] var = new double[n];
		int degrees = (m - 1);
		double c;
		double s;
		for (int j = 0; j < n; j++) {
			c = 0;
			s = 0;
			for (int k = 0; k < m; k++)
				s += v[k][j];
			s = s / m;
			for (int k = 0; k < m; k++)
				c += (v[k][j] - s) * (v[k][j] - s);
			var[j] = c / degrees;
		}
		return var;
	}

	/**
	 * @param matrix
	 * @return a new vector with the column means of matrix
	 */
	static double[] colMeans(double[][] matrix) {
		int rows = matrix.length;
		int cols = matrix[0].length;
		double[] mean = new double[cols];
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				mean[j] += matrix[i][j];
		for (int j = 0; j < cols; j++)
			mean[j] /= (double) rows;
		return mean;
	}

	static double[][] covariance(double[][] matrix) {
		int rows = matrix.length;
		int cols = matrix[0].length;
		double[][] covar = new double[cols][cols];
		int rank = (rows - 1);
		double c;
		double sum1;
		double sum2;
		for (int i = 0; i < cols; i++) {
			for (int j = 0; j < cols; j++) {
				c = 0;
				sum1 = 0;
				sum2 = 0;
				for (int k = 0; k < rows; k++) {
					sum1 += matrix[k][i];
					sum2 += matrix[k][j];
				}
				sum1 = sum1 / rows;
				sum2 = sum2 / rows;
				for (int k = 0; k < rows; k++)
					c += (matrix[k][i] - sum1) * (matrix[k][j] - sum2);
				covar[i][j] = c / rank;
			}
		}
		return covar;
	}

	static double[][] copyRows(double[][] input, int... indices) {
		double[][] matrix = new double[indices.length][input[0].length];
		for (int i = 0; i < indices.length; i++)
			System.arraycopy(input[indices[i]], 0, matrix[i], 0, input[indices[i]].length);
		return matrix;
	}

	/**
	 * @param rows
	 * @param cols
	 * @param fillvalue
	 * @return a new matrix filled with c's
	 */
	static double[][] fillMatrix(int rows, int cols, double fillvalue) {
		double[][] matrix = new double[rows][cols];
		for (int i = 0; i < matrix.length; i++)
			for (int j = 0; j < matrix[i].length; j++)
				matrix[i][j] = fillvalue;
		return matrix;
	}

	static double[][] plus(double[][] m1, double[][] m2) {
		double[][] matrix = new double[m1.length][m1[0].length];
		for (int i = 0; i < m1.length; i++)
			for (int j = 0; j < m1[0].length; j++)
				matrix[i][j] = m1[i][j] + m2[i][j];
		return matrix;
	}

	static double[][] scalarPlus(double[][] m1, double m2) {
		double[][] matrix = new double[m1.length][m1[0].length];
		for (int i = 0; i < m1.length; i++)
			for (int j = 0; j < m1[0].length; j++)
				matrix[i][j] = m1[i][j] + m2;
		return matrix;
	}

	static double[][] minus(double[][] m1, double[][] m2) {
		double[][] matrix = new double[m1.length][m1[0].length];
		for (int i = 0; i < m1.length; i++)
			for (int j = 0; j < m1[0].length; j++)
				matrix[i][j] = m1[i][j] - m2[i][j];
		return matrix;
	}

	static double[][] scalarDivide(double[][] numerator, double denom) {
		double[][] matrix = new double[numerator.length][numerator[0].length];
		for (int i = 0; i < numerator.length; i++)
			for (int j = 0; j < numerator[i].length; j++)
				matrix[i][j] = numerator[i][j] / denom;
		return matrix;
	}

	static double[][] divide(double[][] matrix1, double[]... matrix2) {
		return new QRDecomposition(Matrix.constructWithCopy(matrix2)).solve(Matrix.constructWithCopy(matrix1)).getArray();
	}

	static double[][] scalarMult(double[][] m1, double mul) {
		double[][] matrix = new double[m1.length][m1[0].length];
		for (int i = 0; i < m1.length; i++)
			for (int j = 0; j < m1[i].length; j++)
				matrix[i][j] = m1[i][j] * mul;
		return matrix;
	}

	static double[][] times(double[][] m1, double[][] m2) {
		double[][] array = new double[m1.length][m2[0].length];
		for (int i = 0; i < array.length; i++)
			for (int j = 0; j < array[i].length; j++) {
				double tmp = 0;
				for (int k = 0; k < m1[0].length; k++)
					tmp += m1[i][k] * m2[k][j];
				array[i][j] = tmp;
			}
		return array;
	}

	public static double [][] pca(double [][]X, int no_dims) {
		System.out.println("Preprocessing the data using PCA...");
		double[] stdevX = stddev(X);
		double[] meanX  = colMeans(X);
		double[][] Z    = normalize(X, meanX, stdevX);
		double[][] cov  = covariance(Z);
		EigenvalueDecomposition e = new EigenvalueDecomposition(new Matrix(cov));
		double[][] V    = e.getV().transpose().getArray();
		return copyRows(V,no_dims);
	}

	public static double [][] tsne(double[][] X, int k, int initial_dims, double perplexity) {
		return tsne(X,k,initial_dims, perplexity, 1500, false);
	}

	public static double [][] tsne(double[][] X, int no_dims, int initial_dims, double perplexity, int max_iter, boolean use_pca) {
		//Runs t-SNE on the dataset in the NxD array X to reduce its dimensionality to no_dims dimensions.
		//The syntax of the function is Y = tsne.tsne(X, no_dims, perplexity), where X is an NxD double [][] array."""

		// Initialize variables
		if(use_pca) X = pca(X, initial_dims);
		System.out.println("X:Shape is = " + X.length + " x " + X[0].length);
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
		P = plus(P , transpose(P));
		P = scalarDivide(P , sum(P));
		P = scalarMult(P , 4);					// early exaggeration
		P = maximum(P, 1e-12);

		// Run iterations
		for (int iter = 0; iter < max_iter; iter++) {
			// Compute pairwise affinities
			double [][] sum_Y = transpose(sum(square(Y), 1));
			double [][] num = scalarInverse(scalarPlus(addRowVector(transpose(addRowVector(scalarMult(
					times(Y, transpose(Y)),
					-2),
					sum_Y)),
					sum_Y),
					1));
			assignAtIndex(num, range(n), range(n), 0);
			double [][] Q = scalarDivide(num , sum(num));

			Q = maximum(Q, 1e-12);

			// Compute gradient
			double [][] PQ = minus(P , Q);
			for (int i = 0; i < n; i++) {
				double [][] PQrowi  = copyRows(PQ,i);
				double [][] numrowi = copyRows(num,i);
				dY[i] = sum(scalarMultiply(transpose(tile(scalarMultiply(PQrowi, numrowi), no_dims, 1)) , minus(fillWithRow(Y,i) , Y)), 0)[0];
			}
			// Perform the update
			if (iter < 20)
				momentum = initial_momentum;
			else
				momentum = final_momentum;
			gains = plus(scalarMultiply(scalarPlus(gains,.2),  abs(negate(equal(sign(gains),sign(gains))))),
					// Is the below correct? Above plus, below times...? It is the same in the orig Python impl.
					scalarMultiply(scalarMult(gains,.8), abs(equal(sign(gains),sign(gains)))));

			assignAllLessThan(gains, min_gain, min_gain);
			iY = minus(scalarMult(iY,momentum) , scalarMult(scalarMultiply(gains , dY),eta));
			Y = plus(Y , iY);
			//double [][] tile = tile(mean(Y, 0), n, 1);
			Y = minus(Y , tile(mean(Y, 0), n, 1));

			// Compute current value of cost function
			if (((iter + 1) % 100 == 0) && X.length<100)   {
				double [][] logdivide = log(divide(P , Q));
				logdivide = replaceNaN(logdivide,0);
				double C = sum(times(P , logdivide));
				System.out.println("Iteration " + (iter + 1) + ": error is " + C);
			} else System.out.println("Iteration " + (iter + 1));

			// Stop lying about P-values
			if (iter == 100)
				P = scalarDivide(P , 4);
		}

		// Return solution
		return Y;
	}

	static R Hbeta (double [][] D, double beta){
		double [][] P = exp(scalarMult(scalarMult(D,beta),-1));
		double sumP = sum(P);   // sumP confirmed scalar
		double H = Math.log(sumP) + beta * sum(times(D,transpose(P))) / sumP;
		P = scalarDivide(P,sumP);
		R r = new R();
		r.H = H;
		r.P = P;
		return r;
	}

	static R x2p(double [][] X,double tol, double perplexity){
		int n               = X.length;
		double [][] sum_X   = sum(square(X), 1);
		double [][] times   = scalarMult(times(X, transpose(X)), -2);
		double [][] prodSum = addColumnVector(transpose(times), sum_X);
		double [][] D       = addRowVector(prodSum, transpose(sum_X));
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
					if (Double.isInfinite(betamax)) beta[i] = beta[i] * 2;
					else beta[i] = (beta[i] + betamax)/2;
				} else{
					betamax = beta[i];
					if (Double.isInfinite(betamin))  beta[i] = beta[i]/ 2;
					else beta[i] = ( beta[i] + betamin) / 2;
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

		System.out.println("Mean of sigma summary: " + sigma);

		return r;
	}

	static class R {
		double [][] P;
		@SuppressWarnings("unused")
		double [] beta;
		double H;
	}
}
