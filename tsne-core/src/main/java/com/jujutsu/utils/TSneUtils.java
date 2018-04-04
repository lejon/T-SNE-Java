package com.jujutsu.utils;

import com.jujutsu.tsne.CheckResult;
import com.jujutsu.tsne.PrincipalComponentAnalysis;
import com.jujutsu.tsne.TSneConfig;
import com.jujutsu.tsne.TSneConfiguration;

public class TSneUtils {
	
	public static TSneConfiguration buildConfig(double[][] xin, int outputDims, int initial_dims,
			double perplexity, int max_iter, boolean use_pca, double theta, boolean silent, boolean printError) {
		return new TSneConfig(xin, outputDims, initial_dims, perplexity, max_iter, use_pca, theta, silent, printError);
	}
	
	public static TSneConfiguration buildConfig(double[][] xin, int outputDims, int initial_dims,
			double perplexity, int max_iter, boolean use_pca, double theta, boolean silent) {
		return new TSneConfig(xin, outputDims, initial_dims, perplexity, max_iter, use_pca, theta, silent, true);
	}
	
	public static TSneConfiguration buildConfig(double[][] xin, int outputDims, int initial_dims,
			double perplexity, int max_iter) {
		return new TSneConfig(xin, outputDims, initial_dims, perplexity, max_iter, true, 0.5, false, true);
	}
	
	public CheckResult check(TSneConfiguration parameterObject) {
		int D = parameterObject.getXStartDim();
		double[][] Xin = parameterObject.getXin();
		boolean exact = (parameterObject.getTheta() == .0);
		
		if(exact)
			return new CheckResult(false,
					"The Barnes Hut implementation does not support exact inference yet (theta==0.0), if you want exact t-SNE please use one of the standard t-SNE implementations (FastTSne for instance)");

		if(parameterObject.usePca() && D > parameterObject.getInitialDims() && parameterObject.getInitialDims() > 0) {
			PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
			Xin = pca.pca(Xin, parameterObject.getInitialDims());
		}
		
		int N = parameterObject.getNrRows();
		
		double perplexity = parameterObject.getPerplexity();
		if(N - 1 < 3 * perplexity) { return new CheckResult(false,"Perplexity too large for the number of data points!\n"); }
		
		return new CheckResult(true, "Ok");
	}
}
