package com.jujutsu.tsne;

public class TSneConfig implements TSneConfiguration {
	protected double[][] xin;
	protected int outputDims;
	protected int initial_dims;
	protected double perplexity;
	protected int max_iter;
	protected boolean use_pca;
	protected double theta;
	protected boolean silent;
	protected boolean print_error;

	public TSneConfig(double[][] xin, int outputDims, int initial_dims, double perplexity, int max_iter,
			boolean use_pca, double theta, boolean silent, boolean print_error) {
		this.xin = xin;
		this.outputDims = outputDims;
		this.initial_dims = initial_dims;
		this.perplexity = perplexity;
		this.max_iter = max_iter;
		this.use_pca = use_pca;
		this.theta = theta;
		this.silent = silent;
		this.print_error = print_error;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#getXin()
	 */
	@Override
	public double[][] getXin() {
		return xin;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#setXin(double[][])
	 */
	@Override
	public void setXin(double[][] xin) {
		this.xin = xin;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#getOutputDims()
	 */
	@Override
	public int getOutputDims() {
		return outputDims;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#setOutputDims(int)
	 */
	@Override
	public void setOutputDims(int n) {
		this.outputDims = n;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#getInitialDims()
	 */
	@Override
	public int getInitialDims() {
		return initial_dims;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#setInitialDims(int)
	 */
	@Override
	public void setInitialDims(int initial_dims) {
		this.initial_dims = initial_dims;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#getPerplexity()
	 */
	@Override
	public double getPerplexity() {
		return perplexity;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#setPerplexity(double)
	 */
	@Override
	public void setPerplexity(double perplexity) {
		this.perplexity = perplexity;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#getMaxIter()
	 */
	@Override
	public int getMaxIter() {
		return max_iter;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#setMaxIter(int)
	 */
	@Override
	public void setMaxIter(int max_iter) {
		this.max_iter = max_iter;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#usePca()
	 */
	@Override
	public boolean usePca() {
		return use_pca;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#setUsePca(boolean)
	 */
	@Override
	public void setUsePca(boolean use_pca) {
		this.use_pca = use_pca;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#getTheta()
	 */
	@Override
	public double getTheta() {
		return theta;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#setTheta(double)
	 */
	@Override
	public void setTheta(double theta) {
		this.theta = theta;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#silent()
	 */
	@Override
	public boolean silent() {
		return silent;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#setSilent(boolean)
	 */
	@Override
	public void setSilent(boolean silent) {
		this.silent = silent;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#printError()
	 */
	@Override
	public boolean printError() {
		return print_error;
	}

	/* (non-Javadoc)
	 * @see com.jujutsu.tsne.barneshut.TSneConfiguration#setPrintError(boolean)
	 */
	@Override
	public void setPrintError(boolean print_error) {
		this.print_error = print_error;
	}

	@Override
	public int getXStartDim() {
		return xin[0].length;
	}

	@Override
	public int getNrRows() {
		return xin.length;
	}
	
	
}