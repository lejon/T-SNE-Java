package com.jujutsu.tsne.barneshut;

public class ParallelSPTree extends SPTree {

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

	// Computes edge forces
    void computeEdgeForces(int[] row_P, int[] col_P, double[] val_P, int N, double[] pos_f)
    {
        // Loop over all edges in the graph
        int ind1 = 0;
        for (int n = 0; n < N; n++)
        {
            for (int i = row_P[n]; i < row_P[n + 1]; i++)
            {
                // Compute pairwise distance and Q-value
                double D = 1.0;
                int ind2 = col_P[i] * dimension;
                for (int d = 0; d < dimension; d++)
                {
                    buff[d] = data[ind1 + d] - data[ind2 + d];
                    D += buff[d] * buff[d];
                }
                D = val_P[i] / D;

                // Sum positive force
                for (int d = 0; d < dimension; d++)
                    pos_f[ind1 + d] += D * buff[d];
            }
            ind1 += dimension;
        }
    }
}
