package com.jujutsu.tsne.barneshut;

import static java.lang.Math.max;
import static java.lang.Math.sqrt;

import com.jujutsu.utils.MatrixOps;

public class SPTree {

	  // Fixed constants
    final static int QT_NODE_CAPACITY = 1;
        
	protected SPTree parent;
	protected int dimension;
	protected boolean is_leaf;
	protected int size;
	protected int cum_size;
	
	 // Axis-aligned bounding box stored as a center with half-dimensions to represent the boundaries of this quad tree
    Cell boundary;
    
    // Indices in this space-partitioning tree node, corresponding center-of-mass, and list of all children
    double[] data;
    double[] center_of_mass;
    int [] index = new int[QT_NODE_CAPACITY];
    
    // Children
    SPTree [] children;
    int no_children;
    
    protected double[] buff;

	public SPTree(int D, double[] inp_data, int N) {
		// Compute mean, width, and height of current map (boundaries of SPTree)
		int nD = 0;
		double [] mean_Y = new double [D];
		double []  min_Y = new double [D]; 
		double []  max_Y = new double [D]; 
		for(int d = 0; d < D; d++)  {
			min_Y[d] = Double.POSITIVE_INFINITY;
			max_Y[d] = Double.NEGATIVE_INFINITY;
		}
		for( int n = 0; n < N; n++) {
			for( int d = 0; d < D; d++) {
				mean_Y[d] += inp_data[n * D + d];
				if(inp_data[nD + d] < min_Y[d]) min_Y[d] = inp_data[nD + d];
				if(inp_data[nD + d] > max_Y[d]) max_Y[d] = inp_data[nD + d];
			}
			nD += D;
		}
		for(int d = 0; d < D; d++) mean_Y[d] /= (double) N;

		// Construct SPTree
		double [] width = new double [D];
		for(int d = 0; d < D; d++) width[d] = max(max_Y[d] - mean_Y[d], mean_Y[d] - min_Y[d]) + 1e-5;
		init(null, D, inp_data, mean_Y, width);
		fill(N);
	}

	// Main initialization function
	void init(SPTree inp_parent, int D, double [] inp_data, double [] inp_corner, double [] inp_width)
	{
		parent = inp_parent;
		dimension = D;
		no_children = 2;
		for(int d = 1; d < D; d++) no_children *= 2;
		data = inp_data;
		is_leaf = true;
		size = 0;
		cum_size = 0;

		center_of_mass = new double[D];
		boundary = new Cell(dimension);
		for(int d = 0; d < D; d++) {
			boundary.setCorner(d, inp_corner[d]);
			boundary.setWidth( d, inp_width[d]);
			center_of_mass[d] = .0;
		}

		children = getTreeArray(no_children);
		for(int i = 0; i < no_children; i++) children[i] = null;
		
		buff = new double[dimension];
	}
	
	// Constructor for SPTree with particular size and parent -- build the tree, too!
	SPTree(int D, double [] inp_data, int N, double [] inp_corner, double [] inp_width)
	{
		init(null, D, inp_data, inp_corner, inp_width);
		fill(N);
	}


	// Constructor for SPTree with particular size (do not fill the tree)
	SPTree(int D, double [] inp_data, double [] inp_corner, double [] inp_width)
	{
		init(null, D, inp_data, inp_corner, inp_width);
	}


	// Constructor for SPTree with particular size and parent (do not fill tree)
	SPTree(SPTree inp_parent, int D, double [] inp_data, double [] inp_corner, double [] inp_width) {
		init(inp_parent, D, inp_data, inp_corner, inp_width);
	}


	// Constructor for SPTree with particular size and parent -- build the tree, too!
	SPTree(SPTree inp_parent, int D, double [] inp_data, int N, double [] inp_corner, double [] inp_width)
	{
		init(inp_parent, D, inp_data, inp_corner, inp_width);
		fill(N);
	}

	// Update the data underlying this tree
	void setData(double [] inp_data)
	{
		data = inp_data;
	}


	// Get the parent of the current tree
	SPTree getParent()
	{
		return parent;
	}

	SPTree[] getTreeArray(int no_children) {
		return new SPTree[no_children];
	}

	// Insert a point into the SPTree
	boolean insert(int new_index)
	{
		// Ignore objects which do not belong in this quad tree
		double [] point = MatrixOps.extractRowFromFlatMatrix(data,new_index,dimension);

		if(!boundary.containsPoint(point))
			return false;

		// Online update of cumulative size and center-of-mass
		cum_size++;
		double mult1 = (double) (cum_size - 1) / (double) cum_size;
		double mult2 = 1.0 / (double) cum_size;
		for(int d = 0; d < dimension; d++) {
			center_of_mass[d] *= mult1;
			center_of_mass[d] += mult2 * point[d];
		}

		// If there is space in this quad tree and it is a leaf, add the object here
		if(is_leaf && size < QT_NODE_CAPACITY) {
			index[size] = new_index;
			size++;
			return true;
		}

		// Don't add duplicates for now (this is not very nice)
		boolean any_duplicate = false;
		for(int n = 0; n < size; n++) {
			boolean duplicate = true;
			for(int d = 0; d < dimension; d++) {
				if(point[d] != data[index[n] * dimension + d]) { duplicate = false; break; }
			}
			any_duplicate = any_duplicate || duplicate;
		}
		if(any_duplicate) return true;

		// Otherwise, we need to subdivide the current cell
		if(is_leaf) subdivide();

		// Find out where the point can be inserted
		for(int i = 0; i < no_children; i++) {
			if(children[i].insert(new_index)) return true;
		}

		// Otherwise, the point cannot be inserted (this should never happen)
		assert false;
		return false;
	}

	// Create four children which fully divide this cell into four quads of equal area
	void subdivide() {

		// Create new children
		double [] new_corner = new double[dimension];
		double [] new_width  = new double[dimension];
		for(int i = 0; i < no_children; i++) {
			int div = 1;
			for(int d = 0; d < dimension; d++) {
				new_width[d] = .5 * boundary.getWidth(d);
				if((i / div) % 2 == 1) new_corner[d] = boundary.getCorner(d) - .5 * boundary.getWidth(d);
				else                   new_corner[d] = boundary.getCorner(d) + .5 * boundary.getWidth(d);
				div *= 2;
			}
			children[i] = getNewTree(this, new_corner, new_width);
		}

		// Move existing points to correct children
		for(int i = 0; i < size; i++) {
			boolean success = false;
			for(int j = 0; j < no_children; j++) {
				if(!success) success = children[j].insert(index[i]);
			}
			index[i] = -1;
		}

		// Empty parent node
		size = 0;
		is_leaf = false;
	}

	SPTree getNewTree(SPTree root, double[] new_corner, double[] new_width) {
		return new SPTree(root, dimension, data, new_corner, new_width);
	}

	// Build SPTree on dataset
	void fill(int N)
	{
		for(int i = 0; i < N; i++) insert(i);
	}


	// Checks whether the specified tree is correct
	boolean isCorrect()
	{
		for(int n = 0; n < size; n++) {
			double [] point = MatrixOps.extractRowFromFlatMatrix(data, index[n], dimension);
			if(!boundary.containsPoint(point)) return false;
		}
		if(!is_leaf) {
			boolean correct = true;
			for(int i = 0; i < no_children; i++) correct = correct && children[i].isCorrect();
			return correct;
		}
		else return true;
	}



	// Build a list of all indices in SPTree
	void getAllIndices(int [] indices)
	{
		getAllIndices(indices, 0);
	}


	// Build a list of all indices in SPTree
	int getAllIndices(int [] indices, int loc)
	{

		// Gather indices in current quadrant
		for(int i = 0; i < size; i++) indices[loc + i] = index[i];
		loc += size;

		// Gather indices in children
		if(!is_leaf) {
			for(int i = 0; i < no_children; i++) loc = children[i].getAllIndices(indices, loc);
		}
		return loc;
	}


	int getDepth() {
		if(is_leaf) return 1;
		int depth = 0;
		for(int i = 0; i < no_children; i++) depth = max(depth, children[i].getDepth());
		return 1 + depth;
	}
	
	// Compute non-edge forces using Barnes-Hut algorithm
    double computeNonEdgeForces(int point_index, double theta, double[] neg_f, 
        double buff[], double[] sum_Q)
    {
        //double[] buff = new double[dimension];

        // Make sure that we spend no time on empty nodes or self-interactions
        if (cum_size == 0 || (is_leaf && size == 1 && index[0] == point_index))
            return 0.0;

        // Compute distance between point and center-of-mass
        double D = .0;
        int ind = point_index * dimension;
        // Check whether we can use this node as a "summary"
        double max_width = 0.0;
        double cur_width;
        for (int d = 0; d < dimension; d++)
        {
            buff[d] = data[ind + d] - center_of_mass[d];
            D += buff[d] * buff[d];
            cur_width = boundary.getWidth(d);
            max_width = (max_width > cur_width) ? max_width : cur_width;
        }

        if (is_leaf || max_width / sqrt(D) < theta)
        {
            // Compute and add t-SNE force between point and current node
            D = 1.0 / (1.0 + D);
            double mult = cum_size * D;
            sum_Q[point_index] += mult;
            mult *= D;
            for (int d = 0; d < dimension; d++)
                neg_f[d] += mult * buff[d];
        }
        else
        {

            // Recursively apply Barnes-Hut to children
            for (int i = 0; i < no_children; i++)
                children[i].computeNonEdgeForces(point_index, theta, neg_f,
                    buff, sum_Q);
        }
        return sum_Q[point_index];
    }



	// Computes edge forces
	void computeEdgeForces(int [] row_P, int [] col_P, double [] val_P, int N, double [] pos_f)
	{
		// Loop over all edges in the graph
		double [] buff = new double[dimension];
		int ind1 = 0;
		int ind2 = 0;
		double D;
		for(int n = 0; n < N; n++) {
			for(int i = row_P[n]; i < row_P[n + 1]; i++) {

				// Compute pairwise distance and Q-value
				D = 1.0;
				ind2 = col_P[i] * dimension;
				for(int d = 0; d < dimension; d++) { 
					buff[d] = data[ind1 + d] - data[ind2 + d];
					D += buff[d] * buff[d];
				} 
				D = val_P[i] / D;

				// Sum positive force
				for(int d = 0; d < dimension; d++) pos_f[ind1 + d] += D * buff[d];
			}
			ind1 += dimension;
		}
	}


	// Print out tree
	void print() 
	{
		if(cum_size == 0) {
			System.out.printf("Empty node\n");
			return;
		}

		if(is_leaf) {
			System.out.printf("Leaf node; data = [");
			for(int i = 0; i < size; i++) {
				double [] point = MatrixOps.extractRowFromFlatMatrix(data, index[i], dimension);
				for(int d = 0; d < dimension; d++) System.out.printf("%f, ", point[d]);
				System.out.printf(" (index = %d)", index[i]);
				if(i < size - 1) System.out.printf("\n");
				else System.out.printf("]\n");
			}        
		}
		else {
			System.out.printf("Intersection node with center-of-mass = [");
			for(int d = 0; d < dimension; d++) System.out.printf("%f, ", center_of_mass[d]);
			System.out.printf("]; children are:\n");
			for(int i = 0; i < no_children; i++) children[i].print();
		}
	}
	
	class Cell {		
		int dimension;
		double [] corner;
		double [] width;
		    
		// Constructs cell
		Cell(int inp_dimension) {
			dimension = inp_dimension;
			corner = new double[dimension];
			width  = new double[dimension];
		}

		Cell(int inp_dimension, double [] inp_corner, double [] inp_width) {
			dimension = inp_dimension;
			corner = new double[dimension];
			width  = new double[dimension];
			for(int d = 0; d < dimension; d++) setCorner(d, inp_corner[d]);
			for(int d = 0; d < dimension; d++) setWidth( d,  inp_width[d]);
		}

		double getCorner(int d) {
			return corner[d];
		}

		double getWidth(int d) {
			return width[d];
		}

		void setCorner(int d, double val) {
			corner[d] = val;
		}

		void setWidth(int d, double val) {
			width[d] = val;
		}

		// Checks whether a point lies in a cell
		boolean containsPoint(double point[])
		{
			for(int d = 0; d < dimension; d++) {
				if(corner[d] - width[d] > point[d]) return false;
				if(corner[d] + width[d] < point[d]) return false;
			}
			return true;
		}
	}

}
