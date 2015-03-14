package com.jujutsu.tsne;

/**
 *
 * Author: Leif Jonsson (leif.jonsson@gmail.com)
 * 
 * This is a port of van der Maaten and Hintons Python implementation of t-sne
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
