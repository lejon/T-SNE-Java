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

	double [][] tsne(TSneConfiguration config);
	void abort();

	static class R {
		double [][] P;
		double [] beta;
		double H;
	}
}
