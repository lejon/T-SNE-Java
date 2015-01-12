package com.jujutsu.tsne;

import java.io.File;

import org.junit.Test;
import org.math.io.files.ASCIIFile;
import org.math.io.parser.ArrayString;

public class MatrixOpVsPython {
	
	/*
	 * Notes:
	 * 	Java sum returns a column vector, but does Python perhaps only printing
	 */
		
	@Test
	public void testSum() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.sum(X) = " + TSne.sum(X));
	}

	@Test
	public void testMSum() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.sum(X,0) = \n" + ArrayString.printDoubleArray(TSne.sum(X,0)));
		System.out.println("TSne.sum(X,1) = \n" + ArrayString.printDoubleArray(TSne.sum(X,1)));
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
	public void testTimes() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.times(X,TSne.transpose(X)) = \n" 
		+ ArrayString.printDoubleArray(TSne.times(X,TSne.transpose(X))));
	}

	@Test
	public void testScaleTimes() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.scalarMult(X,-2) = \n" 
		+ ArrayString.printDoubleArray(TSne.scalarMult(X,-2)));
	}

	@Test
	public void testScalarPlus() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.scalarPlus(X,2) = \n" 
		+ ArrayString.printDoubleArray(TSne.scalarPlus(X,2)));
	}
	
	@Test
	public void testScalarInverse() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.scalarInverse(X) = \n" 
		+ ArrayString.printDoubleArray(TSne.scalarInverse(X)));
	}

	@Test
	public void testScalarDivide() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.scalarDivide(X, 2) = \n" 
		+ ArrayString.printDoubleArray(TSne.scalarDivide(X, 2)));
	}

	@Test
	public void testScalarMultiply() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.scalarMultiply(X, X) = \n" 
		+ ArrayString.printDoubleArray(TSne.scalarMultiply(X, X)));
	}

	@Test
	public void testRangeAssign() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		TSne.assignAtIndex(X, TSne.range(4), TSne.range(4), 0);
		System.out.println("assignAtIndex(num, range(n), range(n), 0) = \n" 
		+ ArrayString.printDoubleArray(X));
	}

	@Test
	public void testMinus() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.minus(X, X) = \n" 
		+ ArrayString.printDoubleArray(TSne.minus(X, X)));
	}

	@Test
	public void testTile() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		double [][] PQrowi  = TSne.copyCols(X,4);
		System.out.println("TSne.tile(X, 3, 1) = \n" 
		+ ArrayString.printDoubleArray(TSne.tile(PQrowi, 3, 1)));
		System.out.println("TSne.tile(X, 3, 1) = \n" 
		+ ArrayString.printDoubleArray(TSne.tile(PQrowi, 3, 2)));
	}
	
	@Test
	public void testAssignCol() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		X[3] = TSne.sum(X,0)[0];
		System.out.println("TSne.sum(X,0)[0] = \n" 
		+ ArrayString.printDoubleArray(X));
	}

	@Test
	public void testAssignAllLessThan() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		TSne.assignAllLessThan(X,3,-1);
		System.out.println("TSne.assignAllLessThan(X,3,-1) = \n" 
		+ ArrayString.printDoubleArray(X));
	}
	
	@Test
	public void testSign() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		TSne.assignAllLessThan(X,3,-1);
		System.out.println("TSne.sign(TSne.assignAllLessThan(X,3,-1)) = \n" 
		+ ArrayString.printDoubleArray(TSne.sign(X)));
	}
	
	@Test
	public void testEqual() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));		
		double [][] Y = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		TSne.assignAllLessThan(X,3,-1);
		TSne.assignAllLessThan(Y,2,-1);
		System.out.println("equal(sign(X),sign(Y) ="); 
		printBoolMtx(TSne.equal(TSne.sign(X),TSne.sign(Y)));
	}
	
	public void printBoolMtx(boolean [][]mtx) {
		for (int i = 0; i < mtx.length; i++) {
			for (int j = 0; j < mtx[0].length; j++) {
				System.out.print(mtx[i][j] + ", ");
			}
			System.out.println();
		}
	}

	@Test
	public void testMMean() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		System.out.println("TSne.mean(X,0) = \n" + ArrayString.printDoubleArray(TSne.mean(X,0)));
		System.out.println("TSne.mean(X,1) = \n" + ArrayString.printDoubleArray(TSne.mean(X,1)));
	}

	@Test
	public void testDivide() {
		double [][] X = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));		
		double [][] Y = TSneDemo.nistReadStringDouble(ASCIIFile.read(new File("src/main/python/small_ds.txt")));
		Y = TSne.scalarDivide(Y, 2);
		System.out.println("Y=\n" + ArrayString.printDoubleArray(Y));
		System.out.println(" X / (X/2) =" 
		+ ArrayString.printDoubleArray(TSne.divide(X, Y)));
	}
	
}
