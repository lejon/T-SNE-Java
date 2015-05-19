package com.jujutsu.tsne;

import static org.junit.Assert.assertEquals;

import java.io.File;

import org.junit.Test;

import com.jujutsu.utils.MatrixOps;
import com.jujutsu.utils.MatrixUtils;

public class MatrixOpVsPython {
	
	MatrixOps mo = new MatrixOps();
		
	double epsilon = 0.0000001;
	
	void assertEqualDoubleArrays(double[][] a1, double[][] a2, double tol) {
		for (int i = 0; i < a2.length; i++) {
			for (int j = 0; j < a2[i].length; j++) {
				assertEquals("At [" + i + "," + j + "]", a1[i][j], a2[i][j],tol);
			}
		}
	}

	void assertEqualDoubleVectors(double[] v1, double[] v2, double tol) {
		for (int i = 0; i < v2.length; i++) {
			assertEquals("At [" + i + "]", v1[i], v2[i],tol);
		}
	}

	void assertEqualIntVectors(int[] v1, int[] v2) {
		for (int i = 0; i < v2.length; i++) {
			assertEquals("At [" + i + "]", v1[i], v2[i]);
		}
	}
	
	@Test
	public void testSum() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		System.out.println("TSne.sum(X) = " + MatrixOps.sum(X));
		assertEquals(172.0,MatrixOps.sum(X), epsilon);
	}

	@Test
	public void testMSum() {
		double [] pysum0 = {30.,  30.,  38.,  42.,  32.};
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.sum(X,0) = \n" + ArrayString.printDoubleArray(tsne.sum(X,0)));
		double [][] sum0 = mo.sum(X,0);
		for (int i = 0; i < sum0.length; i++) {
			for (int j = 0; j < sum0[i].length; j++) {
				assertEquals(pysum0[j], sum0[i][j],epsilon);
			}
		}
		double [] pysum1 = {15.,  35.,  19.,  26.,  30.,  18.,  29.};
		//System.out.println("TSne.sum(X,1) = \n" + ArrayString.printDoubleArray(tsne.sum(X,1)));
		double [][] sum1 = mo.sum(X,1);
		for (int i = 0; i < sum1.length; i++) {
			for (int j = 0; j < sum1[i].length; j++) {
				assertEquals(pysum1[i], sum1[i][j],epsilon);
			}
		}
	}

	@Test
	public void testTranspose() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.transpose(X) = \n" 
		//+ ArrayString.printDoubleArray(tsne.transpose(X)));
		double [][] pytranspose = {
				{ 1.,  6.,  3.,  7.,  2.,  3.,  8.},
				{ 2.,  7.,  4.,  3.,  4.,  4.,  6.},
				{ 3.,  8.,  2.,  6.,  7.,  3.,  9.},
				{ 4.,  9.,  7.,  7.,  8.,  3.,  4.},
				{ 5.,  5.,  3.,  3.,  9.,  5.,  2.},};
		double [][] transpose = mo.transpose(X);
		assertEqualDoubleArrays(pytranspose, transpose, epsilon);
	}

	@Test
	public void testSquare() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.square(X) = \n" 
		//+ ArrayString.printDoubleArray(tsne.square(X)));
		double [][] pysquare = {
				{ 1.,   4.,   9.,  16.,  25.},
				{ 36.,  49.,  64.,  81.,  25.},
				{  9.,  16.,   4.,  49.,   9.},
				{ 49.,   9.,  36.,  49.,   9.},
				{  4.,  16.,  49.,  64.,  81.},
				{  9.,  16.,   9.,   9.,  25.},
				{ 64.,  36.,  81.,  16.,   4.},
				};
		double [][] square = mo.square(X);
		assertEqualDoubleArrays(pysquare, square, epsilon);
	}
	
	@Test
	public void testTimes() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.times(X,tsne.transpose(X)) = \n" 
		//+ ArrayString.printDoubleArray(tsne.times(X,tsne.transpose(X))));
		double [][] pydot = {
				{  55.,  105.,   60.,   74.,  108.,   57.,   73.},
				{ 105.,  255.,  140.,  189.,  213.,  122.,  208.},
				{  60.,  140.,   87.,  103.,  119.,   67.,  100.},
				{  74.,  189.,  103.,  152.,  151.,   87.,  162.},
				{ 108.,  213.,  119.,  151.,  214.,  112.,  153.},
				{  57.,  122.,   67.,   87.,  112.,   68.,   97.},
				{  73.,  208.,  100.,  162.,  153.,   97.,  201.}
				};
		double [][] times = mo.times(X,mo.transpose(X));
		assertEqualDoubleArrays(pydot, times, epsilon);
	}

	@Test
	public void testScaleTimes() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.scalarMult(X,-2) = \n" 
		//+ ArrayString.printDoubleArray(tsne.scalarMult(X,-2)));
		double [][] pyscle = 
				{{ -2.,  -4.,  -6.,  -8., -10.},
				 {-12., -14., -16., -18., -10.},
				 { -6.,  -8.,  -4., -14.,  -6.},
				 {-14.,  -6., -12., -14.,  -6.},
				 { -4.,  -8., -14., -16., -18.},
				 { -6.,  -8.,  -6.,  -6., -10.},
				 {-16., -12., -18.,  -8.,  -4.},
				 };
		double [][] scale = mo.scalarMult(X,-2);
		assertEqualDoubleArrays(pyscle, scale, epsilon);
	}

	@Test
	public void testScalarPlus() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.scalarPlus(X,2) = \n" 
		//+ ArrayString.printDoubleArray(tsne.scalarPlus(X,2)));
		double [][] pyplus = 
				{{  3.,   4.,   5.,   6.,   7.},
				 {  8.,   9.,  10.,  11.,   7.},
				 {  5.,   6.,   4.,   9.,   5.},
				 {  9.,   5.,   8.,   9.,   5.},
				 {  4.,   6.,   9.,  10.,  11.},
				 {  5.,   6.,   5.,   5.,   7.},
				 { 10.,   8.,  11.,   6.,   4.},
				 };
		double [][] plus = mo.scalarPlus(X,2);
		assertEqualDoubleArrays(pyplus, plus, epsilon);
	}
	
	@Test
	public void testScalarInverse() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.scalarInverse(X) = \n" 
		//+ ArrayString.printDoubleArray(tsne.scalarInverse(X)));
		double [][] pyinv = 
				{{ 1.,          0.5,         0.33333333,  0.25,        0.2       },
				 { 0.16666667,  0.14285714,  0.125,       0.11111111,  0.2       },
				 { 0.33333333,  0.25,        0.5,         0.14285714,  0.33333333},
				 { 0.14285714,  0.33333333,  0.16666667,  0.14285714,  0.33333333},
				 { 0.5,         0.25,        0.14285714,  0.125,       0.11111111},
				 { 0.33333333,  0.25,        0.33333333,  0.33333333,  0.2       },
				 { 0.125,       0.16666667,  0.11111111,  0.25,        0.5       }
				 };
		double [][] inv = mo.scalarInverse(X);
		assertEqualDoubleArrays(pyinv, inv, epsilon);
	}
	
	@Test
	public void testScalarInverseVector() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.scalarInverse(X[3,:]) = \n" 
		//+ ArrayString.printDoubleArray(tsne.scalarInverse(X[3])));
		double [] pyinv = { 0.14285714,  0.33333333,  0.16666667,  0.14285714,  0.33333333 };
		double [] inv = mo.scalarInverse(X[3]);
		assertEqualDoubleVectors(pyinv, inv, epsilon);
	}

	@Test
	public void testScalarDivide() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.scalarDivide(X, 2) = \n" 
		//+ ArrayString.printDoubleArray(tsne.scalarDivide(X, 2)));
		double [][] div = mo.scalarDivide(X, 2);
		double [][] pydiv = 
				{{ 0.5,  1.,   1.5,  2.,   2.5},
				 { 3.,   3.5,  4.,   4.5,  2.5},
				 { 1.5,  2.,   1.,   3.5,  1.5},
				 { 3.5,  1.5,  3.,   3.5,  1.5},
				 { 1.,   2.,   3.5,  4.,   4.5},
				 { 1.5,  2.,   1.5,  1.5,  2.5},
				 { 4.,   3.,   4.5,  2.,   1. },
				 };
		assertEqualDoubleArrays(pydiv, div, epsilon);
	}

	@Test
	public void testScalarMultiply() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.scalarMultiply(X, X) = \n" 
		//+ ArrayString.printDoubleArray(tsne.scalarMultiply(X, X)));
		double [][] pysm = 
				{{  1.,   4.,   9.,  16.,  25.},
				 { 36.,  49.,  64.,  81.,  25.},
				 {  9.,  16.,   4.,  49.,   9.},
				 { 49.,   9.,  36.,  49.,   9.},
				 {  4.,  16.,  49.,  64.,  81.},
				 {  9.,  16.,   9.,   9.,  25.},
				 { 64.,  36.,  81.,  16.,   4.},
				 };
		double [][] sm = mo.scalarMultiply(X, X);
		assertEqualDoubleArrays(pysm, sm, epsilon);
	}

	@Test
	public void testRangeAssign() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		mo.assignAtIndex(X, MatrixOps.range(4), MatrixOps.range(4), 0);
		//System.out.println("assignAtIndex(num, range(n), range(n), 0) = \n" 
		//+ ArrayString.printDoubleArray(X));
		double [][] pyasgn = 
				{{ 0.,  2.,  3.,  4.,  5.},
				 { 6.,  0.,  8.,  9.,  5.},
				 { 3.,  4.,  0.,  7.,  3.},
				 { 7.,  3.,  6.,  0.,  3.},
				 { 2.,  4.,  7.,  8.,  9.},
				 { 3.,  4.,  3.,  3.,  5.},
				 { 8.,  6.,  9.,  4.,  2.},
				 };
		assertEqualDoubleArrays(pyasgn, X, epsilon);
	}

	@Test
	public void testMinus() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.minus(X, X) = \n" 
		//+ ArrayString.printDoubleArray(tsne.sMinus(X, X)));
		double [][] pymin = 
				{{ 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 };
		double [][] min = mo.sMinus(X, X);
		assertEqualDoubleArrays(pymin, min, epsilon);
	}

	@Test
	public void testParMinus() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.minus(X, X) = \n" 
		//+ ArrayString.printDoubleArray(tsne.parScalarMinus(X, X)));
		double [][] pymin = 
				{{ 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 { 0.,  0.,  0.,  0.,  0.},
				 };
		double [][] min = mo.parScalarMinus(X, X);
		assertEqualDoubleArrays(pymin, min, epsilon);
	}

	@Test
	public void testTile() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		double [][] PQrowi  = mo.copyCols(X,4);
		//System.out.println("TSne.tile(X, 3, 1) = \n" 
		//+ ArrayString.printDoubleArray(tsne.tile(PQrowi, 3, 1)));
		double [][] pytile1 = 
				{{ 5.,  5.,  3.,  3.,  9.,  5.,  2.,},
				 { 5.,  5.,  3.,  3.,  9.,  5.,  2.,},
				 { 5.,  5.,  3.,  3.,  9.,  5.,  2.,},
				};
		double [][] tile1 = mo.tile(PQrowi, 3, 1);
		assertEqualDoubleArrays(pytile1, tile1, epsilon);
		//System.out.println("TSne.tile(X, 3, 2) = \n" 
		//+ ArrayString.printDoubleArray(tsne.tile(PQrowi, 3, 2)));
		double [][] pytile2 =
				{{ 5.,  5.,  3.,  3.,  9.,  5.,  2.,  5.,  5.,  3.,  3.,  9.,  5.,  2.},
				 { 5.,  5.,  3.,  3.,  9.,  5.,  2.,  5.,  5.,  3.,  3.,  9.,  5.,  2.},
				 { 5.,  5.,  3.,  3.,  9.,  5.,  2.,  5.,  5.,  3.,  3.,  9.,  5.,  2.},
				 };
		double [][] tile2 = mo.tile(PQrowi, 3, 2);
		assertEqualDoubleArrays(pytile2, tile2, epsilon);
	}
	
	@Test
	public void testAssignCol() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		X[3] = mo.sum(X,0)[0];
		//System.out.println("TSne.sum(X,0)[0] = \n" 
		//+ ArrayString.printDoubleArray(X));
		double [][] pyasgn = 
				{{  1.,   2.,   3.,   4.,   5.},
				 {  6.,   7.,   8.,   9.,   5.},
				 {  3.,   4.,   2.,   7.,   3.},
				 { 30.,  30.,  38.,  42.,  32.},
				 {  2.,   4.,   7.,   8.,   9.},
				 {  3.,   4.,   3.,   3.,   5.},
				 {  8.,   6.,   9.,   4.,   2.},
				 };
		assertEqualDoubleArrays(pyasgn, X, epsilon);
	}

	@Test
	public void testAssignAllLessThan() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		mo.assignAllLessThan(X,3,-1);
		//System.out.println("TSne.assignAllLessThan(X,3,-1) = \n" 
		//+ ArrayString.printDoubleArray(X));
		double [][] pylt =
				{{-1., -1.,  3.,  4.,  5.},
				 { 6.,  7.,  8.,  9.,  5.},
				 { 3.,  4., -1.,  7.,  3.},
				 { 7.,  3.,  6.,  7.,  3.},
				 {-1.,  4.,  7.,  8.,  9.},
				 { 3.,  4.,  3.,  3.,  5.},
				 { 8.,  6.,  9.,  4., -1.},
				 };
		assertEqualDoubleArrays(pylt, X, epsilon);
	}
	
	@Test
	public void testSign() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		mo.assignAllLessThan(X,3,-1);
		//System.out.println("TSne.sign(TSne.assignAllLessThan(X,3,-1)) = \n" 
		//+ ArrayString.printDoubleArray(tsne.sign(X)));
	}
	
	@Test
	public void testEqual() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");		
		double [][] Y = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		mo.assignAllLessThan(X,3,-1);
		mo.assignAllLessThan(Y,2,-1);
		System.out.println("equal(sign(X),sign(Y) ="); 
		printBoolMtx(mo.equal(mo.sign(X),mo.sign(Y)));
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
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		//System.out.println("TSne.mean(X,0) = \n" + ArrayString.printDoubleArray(tsne.mean(X,0)));
		double [] pymean0 = { 4.28571429,  4.28571429,  5.42857143,  6.,          4.57142857};
		double [][] mean0 = MatrixOps.mean(X,0);
		assertEqualDoubleVectors(pymean0, mean0[0], epsilon);
		double [] pymean1 = {3.,   7.,   3.8,  5.2,  6.,   3.6,  5.8};
		//System.out.println("TSne.mean(X,1) = \n" + ArrayString.printDoubleArray(tsne.mean(X,1)));
		double [][] mean1mtrx = MatrixOps.mean(X,1);
		double [] mean1 = new double [mean1mtrx.length];
		for (int i = 0; i < mean1mtrx.length; i++) {
			for (int j = 0; j < mean1mtrx[i].length; j++) {				
				mean1[i] = mean1mtrx[i][j];
			}
		}
		assertEqualDoubleVectors(pymean1, mean1, epsilon);
	}
	
	@Test
	public void testVMean() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		System.out.println("TSne.mean(X[3,:]) = \n" + mo.mean(X[3]));
		assertEquals(5.2, mo.mean(X[3]), epsilon);
	}

	@Test
	public void testElementWiseDivide() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");		
		double [][] Y = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		Y = mo.scalarDivide(Y, 2);
		//System.out.println("Y=\n" + ArrayString.printDoubleArray(Y));
		//System.out.println(" X / (X/2) =" 
		//+ ArrayString.printDoubleArray(tsne.scalarDivide(X, Y)));
		double [][] pydiv = 
			{{ 2.,  2.,  2.,  2.,  2.},
			 { 2.,  2.,  2.,  2.,  2.},
			 { 2.,  2.,  2.,  2.,  2.},
			 { 2.,  2.,  2.,  2.,  2.},
			 { 2.,  2.,  2.,  2.,  2.},
			 { 2.,  2.,  2.,  2.,  2.},
			 { 2.,  2.,  2.,  2.,  2.},
			};
		double [][] div = mo.scalarDivide(X, Y);
		assertEqualDoubleArrays(pydiv, div, epsilon);
	}
	
	@Test
	public void testSqrt() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");
		double [] pysqrt = { 2.64575131,  1.73205081,  2.44948974,  2.64575131,  1.73205081};
		//System.out.println("sqrt(X[3,:]) =\n" + ArrayString.printDoubleArray(tsne.sqrt(X[3])));
		assertEqualDoubleVectors(pysqrt, MatrixOps.sqrt(X[3]), epsilon);	
	}
	
	@Test
	public void testExp() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");		
		//System.out.println(" exp(X) =" 
		//+ ArrayString.printDoubleArray(tsne.exp(X)));
		double [][] pyexp = 
				{{  2.71828183e+00,   7.38905610e+00,   2.00855369e+01,   5.45981500e+01,	    1.48413159e+02},
				 {  4.03428793e+02,   1.09663316e+03,   2.98095799e+03,   8.10308393e+03,	    1.48413159e+02},
				 {  2.00855369e+01,   5.45981500e+01,   7.38905610e+00,   1.09663316e+03,	    2.00855369e+01},
				 {  1.09663316e+03,   2.00855369e+01,   4.03428793e+02,   1.09663316e+03,	    2.00855369e+01},
				 {  7.38905610e+00,   5.45981500e+01,   1.09663316e+03,   2.98095799e+03,	    8.10308393e+03},
				 {  2.00855369e+01,   5.45981500e+01,   2.00855369e+01,   2.00855369e+01,	    1.48413159e+02},
				 {  2.98095799e+03,   4.03428793e+02,   8.10308393e+03,   5.45981500e+01,	    7.38905610e+00}};
		double [][] jexp = mo.exp(X);
		assertEqualDoubleArrays(pyexp, jexp, 0.00001);
	}
	
	@Test
	public void testLog() {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/python/small_ds.txt"), " ");		
		//System.out.println(" log(X) =" 
		//+ ArrayString.printDoubleArray(tsne.log(X)));
		double [][] pylog = 
				{{ 0.,          0.69314718,  1.09861229,  1.38629436,  1.60943791},
				 { 1.79175947,  1.94591015,  2.07944154,  2.19722458,  1.60943791},
				 { 1.09861229,  1.38629436,  0.69314718,  1.94591015,  1.09861229},
				 { 1.94591015,  1.09861229,  1.79175947,  1.94591015,  1.09861229},
				 { 0.69314718,  1.38629436,  1.94591015,  2.07944154,  2.19722458},
				 { 1.09861229,  1.38629436,  1.09861229,  1.09861229,  1.60943791},
				 { 2.07944154,  1.79175947,  2.19722458,  1.38629436,  0.69314718},
				 };
		double [][] jlog = MatrixOps.log(X);
		assertEqualDoubleArrays(pylog, jlog, epsilon);
	}
		
	@Test 
	public void testConcatenate() {
		int [] v1 = {1,2,3,4};
		int [] v2 = {3,4,5,6};
		//System.out.println(ArrayString.printIntArray(v1) + " + " +  ArrayString.printIntArray(v2) + " = "
		//		+ ArrayString.printIntArray(tsne.concatenate(v1, v2)));
		int [] v3 = mo.concatenate(v1, v2);
		int [] expct = {1,2,3,4,3,4,5,6};
		assertEqualIntVectors(expct, v3);
	}
	
	@Test 
	public void testDiagWide() {
		double [][] diag = {{1,2,3,4,5}};
		double [][] dmatrix = mo.diag(diag);
		for (int i = 0; i < dmatrix.length; i++) {
			for (int j = 0; j < dmatrix[0].length; j++) {
				if(i==j) assertEquals(dmatrix[i][j],diag[0][i],epsilon);
			}
		}
	}

	@Test 
	public void testDiagLong() {
		double [][] diag = {{1},{2},{3},{4},{5}};
		double [][] dmatrix = mo.diag(diag);
		for (int i = 0; i < dmatrix.length; i++) {
			for (int j = 0; j < dmatrix[0].length; j++) {
				if(i==j) assertEquals(dmatrix[i][j],diag[i][0],epsilon);
			}
		}
	}

}
