package com.jujutsu.tsne;

/**
*
* Author: Leif Jonsson (leif.jonsson@gmail.com)
* 
* This is a Java implementation of van der Maaten and Hintons t-sne 
* dimensionality reduction technique that is particularly well suited 
* for the visualization of high-dimensional datasets
*
*/
public interface TSne {

	double [][] tsne(double[][] X, int k, int initial_dims, double perplexity);
	double [][] tsne(double[][] X, int k, int initial_dims, double perplexity, int maxIterations);

	double [][] tsne(double[][] X, int no_dims, int initial_dims, double perplexity, int max_iter, boolean use_pca);

	R Hbeta (double [][] D, double beta);
	
	R x2p(double [][] X,double tol, double perplexity);

	static class R {
		double [][] P;
		double [] beta;
		double H;
	}
}
