package com.jujutsu.tsne.barneshut;

import static java.lang.Math.sqrt;

public class ParallelSPTree extends SPTree {

//	private static final int PAR_LIMIT = 100;
//
//	private static class ThreadPoolHolder {
//		public static final ExecutorService dimCalcExecutor = Executors.newFixedThreadPool(4);
//	}
//
//	public static ExecutorService getThreadPool() {
//		return ThreadPoolHolder.dimCalcExecutor;
//	}

	public ParallelSPTree(int D, double[] inp_data, int N) {
		super(D, inp_data, N);
	}

	public ParallelSPTree(int D, double[] inp_data, int N, double[] inp_corner, double[] inp_width) {
		super(D, inp_data, N, inp_corner, inp_width);
	}

	public ParallelSPTree(int D, double[] inp_data, double[] inp_corner, double[] inp_width) {
		super(D, inp_data, inp_corner, inp_width);
	}

	public ParallelSPTree(SPTree inp_parent, int D, double[] inp_data, double[] inp_corner, double[] inp_width) {
		super(inp_parent, D, inp_data, inp_corner, inp_width);
	}

	public ParallelSPTree(SPTree inp_parent, int D, double[] inp_data, int N, double[] inp_corner, double[] inp_width) {
		super(inp_parent, D, inp_data, N, inp_corner, inp_width);
	}

	@Override
	SPTree[] getTreeArray(int no_children) {
		return new ParallelSPTree[no_children];
	}

	@Override
	SPTree getNewTree(SPTree root, double[] new_corner, double[] new_width) {
		return new ParallelSPTree(root, dimension, data, new_corner, new_width);
	}

	// Compute non-edge forces using Barnes-Hut algorithm
//	@Override
//	double computeNonEdgeForces(int point_index, double theta, double [] neg_f, Object accumulator)
//	{
//		if(dimension < PAR_LIMIT) {
//			return computeNonEdgeForcesStd( point_index,  theta, neg_f,  accumulator);
//		} else {
//			return computeNonEdgeForcesPar( point_index,  theta, neg_f,  accumulator);
//		}
//	}

	//double computeNonEdgeForcesStd(int point_index, double theta, double [] neg_f, Object accumulator)
	@Override
	double computeNonEdgeForces(int point_index, double theta, double [] neg_f, Object accumulator)
	{
		Double sum_Q = (Double) accumulator;
		double input_sum_Q = sum_Q;
		double [] buff = new double[dimension];

		// Make sure that we spend no time on empty nodes or self-interactions
		if(cum_size == 0 || (is_leaf && size == 1 && index[0] == point_index)) return 0.0;

		// Compute distance between point and center-of-mass
		double D = .0;
		int ind = point_index * dimension;
		double max_width = 0.0;
		//double cur_width;
		for(int d = 0; d < dimension; d++) {
			buff[d] = data[ind + d] - center_of_mass[d];
			D += buff[d] * buff[d];
			// Check whether we can use this node as a "summary"
			double cur_width = boundary.getWidth(d);
			max_width = (max_width > cur_width) ? max_width : cur_width;
		}

		if(is_leaf || max_width / sqrt(D) < theta) {
			// Compute and add t-SNE force between point and current node
			D = 1.0 / (1.0 + D);
			double mult = cum_size * D;
			sum_Q += mult;
			mult *= D;
			for(int d = 0; d < dimension; d++) neg_f[d] += mult * buff[d];
		}
		else {
			// Recursively apply Barnes-Hut to children
			for(int i = 0; i < no_children; i++) sum_Q += children[i].computeNonEdgeForces(point_index, theta, neg_f, input_sum_Q);
		}

		return sum_Q;
	}

//	double computeNonEdgeForcesPar(int point_index, double theta, double [] neg_f, Object accumulator)
//	{
//		Double sum_Q = (Double) accumulator;
//		double input_sum_Q = sum_Q;
//		double [] buff = new double[dimension];
//
//		// Make sure that we spend no time on empty nodes or self-interactions
//		if(cum_size == 0 || (is_leaf && size == 1 && index[0] == point_index)) return 0.0;
//
//		// Compute distance between point and center-of-mass
//		double D = .0;
//		int ind = point_index * dimension;
//		double max_width = 0.0;
//		//double cur_width;
//		List<Callable<double []>> jobs = new ArrayList<>();
//		int start = 0;
//		int splitSize = dimension / PAR_LIMIT;
//		int splitLeftover = dimension % PAR_LIMIT;
//		int end  = splitSize + splitLeftover;
//		int threadCnt = 0;
//		System.out.println("Splitsize = " + splitSize + " leftover = " + splitLeftover);
//		while( end < dimension) {
//			final int tStart = start;
//			final int tEnd = end;
//			jobs.add((Callable<double []>)() -> {
//				return calcDim(buff, ind, tStart, tEnd);
//			});
//			start = end;
//			end = end + splitSize;
//			System.out.println("Started " + threadCnt++ + " calls...");
//		}
//
//		List<Future<double []>> results;
//		try {
//			results = getThreadPool().invokeAll(jobs);
//			for (Future<double []> result : results) {
//				double [] res = result.get();
//				D += res[0];
//				max_width = (max_width > res[1]) ? max_width : res[1];
//			}
//		} catch (InterruptedException e) {
//			e.printStackTrace();
//			System.exit(-1);
//		} catch (ExecutionException e) {
//			e.printStackTrace();
//			System.exit(-1);
//		}
//
//		if(is_leaf || max_width / sqrt(D) < theta) {
//			// Compute and add t-SNE force between point and current node
//			D = 1.0 / (1.0 + D);
//			double mult = cum_size * D;
//			sum_Q += mult;
//			mult *= D;
//			for(int d = 0; d < dimension; d++) neg_f[d] += mult * buff[d];
//		}
//		else {
//			// Recursively apply Barnes-Hut to children
//			for(int i = 0; i < no_children; i++) sum_Q += children[i].computeNonEdgeForces(point_index, theta, neg_f, input_sum_Q);
//		}
//
//		return sum_Q;
//	}
//
//	private double [] calcDim(double [] buff, int ind, int start, int end ) {
//		double [] res = new double[2]; 
//		for(int d = start; d < end; d++) {
//			buff[d] = data[ind + d] - center_of_mass[d];
//			res[0] += buff[d] * buff[d];
//			// Check whether we can use this node as a "summary"
//			double cur_width = boundary.getWidth(d);
//			res[1] = (res[1] > cur_width) ? res[1] : cur_width;
//		}
//		return res;	
//	}

}
