package com.jujutsu.tsne;

import static org.junit.Assert.*;

import org.junit.Test;

public class MatrixOpTest {

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
		double [][] tr1 = TSne.origtranspose(matrix);
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
			double [][] tr1 = TSne.origtranspose(matrix);
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
			double [][] tr1 = TSne.origtranspose(matrix);
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


}
