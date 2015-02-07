package com.jujutsu.tsne;

import static org.junit.Assert.assertEquals;

import java.io.File;

import org.junit.Test;
import org.math.io.files.ASCIIFile;
import org.math.io.parser.ArrayString;

import Jama.Matrix;

public class MatrixOpTest {
		
	static double[][] multiply(double[][] matrix1, double[]... matrix2) {
		//return new QRDecomposition(Matrix.constructWithCopy(matrix2)).solve(Matrix.constructWithCopy(matrix1)).getArray();
		Matrix A = Matrix.constructWithCopy(matrix2);
		Matrix B = Matrix.constructWithCopy(matrix1); 
		return  A.times(B).getArray();
	}

	@Test
	public void testSum() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.sum(X) = " + TSne.sum(X));
	}

	@Test
	public void testMSum() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.msum(X) = " + TSne.sum(X,1));
	}

	@Test
	public void testTranspose() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.transpose(X) = \n" 
		+ ArrayString.printDoubleArray(TSne.transpose(X)));
	}

	@Test
	public void testSquare() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.square(X) = \n" 
		+ ArrayString.printDoubleArray(TSne.square(X)));
	}

	
	@Test
	public void testPCA() {
        double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/resources/datasets/nist_pca_2.txt")));
        System.out.println(ArrayString.printDoubleArray(X));
        //double [][] nmatrix = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/resources/datasets/mnist2500_X.txt")));
        double [][] matrix = TSne.rnorm(200,7);
        System.out.println(ArrayString.printDoubleArray(matrix));
        PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
        double [][] pcad = pca.pca(matrix, 2);
        System.out.println(ArrayString.printDoubleArray(pcad));
	}
	
	@Test
	public void testTimes() {
		double [][] A =  {{1,2,3},{3,2,1},{2,1,3}};
		double [][] B = {{4,5,6},{6,5,4},{4,6,5}};
		double [][] AxB = TSne.times(A, B);
		double [][] AxBM = Matrix.constructWithCopy(A).times(Matrix.constructWithCopy(B)).getArray();
		
		System.out.println("Matrices:");
		for (int i = 0; i < AxB.length; i++) {
			for (int j = 0; j < AxB[0].length; j++) {
				System.out.print(AxBM[i][j] + " => " + AxB[i][j] + ", ");
				assertEquals(AxBM[i][j],AxB[i][j],0.0000001);
				
			}
			System.out.println();
		}
	}

	@Test
	public void testTransposes() {
		int rows = 15302;
		int cols = 1143;
		double [][] matrix = new double [rows][cols];
		double [][] trmatrix = new double [cols][rows];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				matrix[i][j] = (double) j + (i*cols);
				trmatrix[j][i] = (double) j + (i*cols);
			}
		}
		double [][] tr1 = TSne.naivetranspose(matrix);
		assertEquals(tr1.length,cols);
		assertEquals(tr1[0].length,rows);
		double [][] tr2 = TSne.transpose(matrix,20);
		assertEquals(tr2.length,cols);
		assertEquals(tr2[0].length,rows);
		for (int i = 0; i < tr1.length; i++) {
			for (int j = 0; j < tr1[0].length; j++) {
				assertEquals(trmatrix[i][j],tr1[i][j],0.0000001);
				assertEquals("I: " + i + " J:" + j, trmatrix[i][j],tr2[i][j],0.0000001);
			}
		}
	}

	@Test
	public void testManyTransposes() {
		int rows = 15302;
		int cols = 1143;
		double [][] matrix = new double [rows][cols];
		double [][] trmatrix = new double [cols][rows];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				matrix[i][j] = (double) j + (i*cols);
				trmatrix[j][i] = (double) j + (i*cols);
			}
		}
		int noLaps = 10;
		long trtime = 0;
		long partrtime = 0;
		long time = 0;
		for (int laps = 0; laps < noLaps; laps++) {
			time = System.currentTimeMillis();
			double [][] tr1 = TSne.naivetranspose(matrix);
			trtime += (System.currentTimeMillis()-time);
			assertEquals(tr1.length,cols);
			assertEquals(tr1[0].length,rows);
			time = System.currentTimeMillis();
			double [][] tr2 = TSne.transpose(matrix,20);
			partrtime += (System.currentTimeMillis()-time);
			assertEquals(tr2.length,cols);
			assertEquals(tr2[0].length,rows);
			for (int i = 0; i < tr1.length; i++) {
				for (int j = 0; j < tr1[0].length; j++) {
					assertEquals(trmatrix[i][j],tr1[i][j],0.0000001);
					assertEquals("I: " + i + " J:" + j, trmatrix[i][j],tr2[i][j],0.0000001);
				}
			}
		}
		System.out.println("    Tr time: " + trtime);
		System.out.println("Par Tr time: " + partrtime);
	}

	public void timeTransposes() {
		int rows = 15302;
		int cols = 1143;
		double [][] matrix = new double [rows][cols];
		double [][] trmatrix = new double [cols][rows];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				matrix[i][j] = (double) j + (i*cols);
				trmatrix[j][i] = (double) j + (i*cols);
			}
		}
		int noLaps = 100;
		long trtime = 0;
		long partrtime = 0;
		long time = 0;
		for (int laps = 0; laps < noLaps; laps++) {
			time = System.currentTimeMillis();
			double [][] tr1 = TSne.naivetranspose(matrix);
			trtime += (System.currentTimeMillis()-time);
			assertEquals(tr1.length,cols);
			assertEquals(tr1[0].length,rows);
			time = System.currentTimeMillis();
			double [][] tr2 = TSne.transpose(matrix,20);
			partrtime += (System.currentTimeMillis()-time);
			assertEquals(tr2.length,cols);
			assertEquals(tr2[0].length,rows);
			for (int i = 0; i < tr1.length; i++) {
				for (int j = 0; j < tr1[0].length; j++) {
					assertEquals(trmatrix[i][j],tr1[i][j],0.0000001);
					assertEquals("I: " + i + " J:" + j, trmatrix[i][j],tr2[i][j],0.0000001);
				}
			}
		}
		System.out.println("    Tr time: " + trtime);
		System.out.println("Par Tr time: " + partrtime);
	}
	
	@Test
	public void timeTransposesNist() {
		double [][] matrix = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/resources/datasets/mnist2500_X.txt")));
		int rows = matrix.length;
		int cols = matrix[0].length;
		double [][] trmatrix = new double [cols][rows];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				trmatrix[j][i] = matrix[i][j];
			}
		}
		int noLaps = 500;
		long trtime = 0;
		long partrtime = 0;
		long time = 0;
		System.out.println("Size is" + rows + " x" + cols + "...");
		for (int laps = 0; laps < noLaps; laps++) {
			if((laps%100)==0) System.out.println("Iter " + laps + "...");
			time = System.currentTimeMillis();
			double [][] tr1 = TSne.naivetranspose(matrix);
			trtime += (System.currentTimeMillis()-time);
			assertEquals(tr1.length,cols);
			assertEquals(tr1[0].length,rows);
			time = System.currentTimeMillis();
			double [][] tr2 = TSne.transpose(matrix,20);
			partrtime += (System.currentTimeMillis()-time);
			assertEquals(tr2.length,cols);
			assertEquals(tr2[0].length,rows);
			for (int i = 0; i < tr1.length; i++) {
				for (int j = 0; j < tr1[0].length; j++) {
					assertEquals(trmatrix[i][j],tr1[i][j],0.0000001);
					assertEquals("I: " + i + " J:" + j, trmatrix[i][j],tr2[i][j],0.0000001);
				}
			}
		}
		System.out.println("    Tr time: " + trtime);
		System.out.println("Par Tr time: " + partrtime);
	}
	
	@Test
	public void timeScalarMultNist() {
		double [][] matrix1 = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/resources/datasets/mnist2500_X.txt")));
		double [][] matrix2 = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/resources/datasets/mnist2500_X.txt")));
		int rows = matrix1.length;
		int cols = matrix1[0].length;
		int noLaps = 1000;
		long trtime = 0;
		long partrtime = 0;
		long time = 0;
		System.out.println("Size is " + rows + " x " + cols + "...");
		for (int laps = 0; laps < noLaps; laps++) {
			if((laps%100)==0) System.out.println("Iter " + laps + "...");
			time = System.currentTimeMillis();
			double [][] tr1 = TSne.scalarMultiply(matrix1, matrix2);
			trtime += (System.currentTimeMillis()-time);
			time = System.currentTimeMillis();
			double [][] tr2 = TSne.parScalarMultiply(matrix1, matrix2);
			partrtime += (System.currentTimeMillis()-time);
			for (int i = 0; i < tr1.length; i++) {
				for (int j = 0; j < tr1[0].length; j++) {
					assertEquals("I: " + i + " J:" + j, tr1[i][j],tr2[i][j],0.0000001);
				}
			}
		}
		System.out.println("    Tr time: " + trtime);
		System.out.println("Par Tr time: " + partrtime);
	}

}
