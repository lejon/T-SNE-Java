package com.jujutsu.tsne.barneshut;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.concurrent.ThreadLocalRandom;

public class VpTree<StorageType, DistanceType> {

	double _tau;
	DataPoint [] _items;
	Node _root;

	public void create(DataPoint [] items) {
		_items = items.clone();
		_root = buildFromPoints(0,items.length);
	}

	public void search(DataPoint target, int k, List<DataPoint> results, List<Double> distances) {
		// Use a priority queue to store intermediate results on
		// Javas prio heap is by default in ascending order, we want descending... 
		PriorityQueue<HeapItem> heap = new PriorityQueue<HeapItem>(k,new Comparator<HeapItem>() {
			@Override
			public int compare(HeapItem o1, HeapItem o2) {
				return -1 * o1.compareTo(o2);
			}
		}); 
        
        // Variable that tracks the distance to the farthest point in our results
        _tau = Double.MAX_VALUE;
        
        // Perform the search
       _root.search(_root, target, k, heap);
        
        // Gather final results
        results.clear(); 
        distances.clear();
        while(!heap.isEmpty()) {
            results.add(_items[heap.peek().index]);
            distances.add(heap.peek().dist);
            heap.remove();
        }
        
        // Results are in reverse order 
        Collections.reverse(results);
        Collections.reverse(distances);
	}

	// Function that (recursively) fills the tree
	public Node buildFromPoints( int lower, int upper )
	{
		if (upper == lower) {     // indicates that we're done here!
			return null;
		}

		// Lower index is center of current node
		Node node = new Node();
		node.index = lower;

		if (upper - lower > 1) {      // if we did not arrive at leaf yet

			// Choose an arbitrary point and move it to the start
			int i = (int) (ThreadLocalRandom.current().nextDouble() * (upper - lower - 1)) + lower;
			swap(_items, lower, i);

			// Partition around the median distance
			int median = (upper + lower) / 2;
			nth_element(_items, lower + 1,	median,	upper, new DistanceComparator(_items[lower]));

			// Threshold of the new node will be the distance to the median
			node.threshold = distance(_items[lower], _items[median]);

			// Recursively build tree
			node.index = lower;
			node.left = buildFromPoints(lower + 1, median);
			node.right = buildFromPoints(median, upper);
		}

		// Return result
		return node;
	}
	
	public Node getRoot() {
		return _root;
	}
	
	// Quick and dirty... optimize later :D
	static void nth_element(DataPoint [] array, int low, int mid, int high,
			DistanceComparator distanceComparator) {
		DataPoint [] tmp = new DataPoint[high-low];
		for (int i = 0; i < tmp.length; i++) {
			tmp[i] = array[low+i];
		}
		Arrays.sort(tmp, distanceComparator);
		for (int i = 0; i < tmp.length; i++) {
			array[low+i] = tmp[i];
		}
	}
	
	static void nth_element(int [] array, int low, int mid, int high) {
		int [] tmp = new int[high-low];
		for (int i = 0; i < tmp.length; i++) {
			tmp[i] = array[low+i];
		}
		Arrays.sort(tmp);
		for (int i = 0; i < tmp.length; i++) {
			array[low+i] = tmp[i];
		}
	}

	private double distance(DataPoint dataPoint, DataPoint dataPoint2) {
		return dataPoint.euclidean_distance(dataPoint2);
	}

	private void swap(DataPoint [] items, int idx1,int idx2) {
		DataPoint dp = items[idx1];
		items[idx1] = items[idx2];
		items[idx2] = dp;
	}

	// An item on the intermediate result queue
    static class HeapItem implements Comparable<HeapItem> {
    	int index;
    	double dist;
    	HeapItem( int index, double dist) {
    		this.index = index;
    		this.dist = dist; 
    	}

    	@Override
		public int compareTo(HeapItem o) {
			return dist < o.dist ? -1 : (dist > o.dist ? 1 : 0);
		}
    	
    	@Override
    	public String toString() {
    		return "HeapItem (index=" + index + ",dist=" + dist + ")"; 
    	}
    };

	class Node {
		int index;
		double threshold;
		private Node left;
		private Node right;
		
		@Override
		public String toString() {
			return "Node(id=" + index + ")";
		}
		
		public Node getLeft() {
			return left;
		}

		public Node getRight() {
			return right;
		}

		// Helper function that searches the tree    
		void search(Node node, DataPoint target, int k, PriorityQueue<HeapItem> heap)
		{
			if(node == null) return;     // indicates that we're done here

			// Compute distance between target and current node
			double dist = distance(_items[node.index], target);

			// If current node within radius tau
			if(dist < _tau) {
				if(heap.size() == k) heap.remove();           // remove furthest node from result list (if we already have k results)
				heap.add(new HeapItem(node.index, dist));     // add current node to result list
				if(heap.size() == k) _tau = heap.peek().dist; // update value of tau (farthest point in result list)
			}

			// Return if we arrived at a leaf
			if(node.left == null && node.right == null) {
				return;
			}

			// If the target lies within the radius of ball
			if(dist < node.threshold) {
				if(dist - _tau <= node.threshold) {         // if there can still be neighbors inside the ball, recursively search left child first
					search(node.left, target, k, heap);
				}

				if(dist + _tau >= node.threshold) {         // if there can still be neighbors outside the ball, recursively search right child
					search(node.right, target, k, heap);
				}

				// If the target lies outsize the radius of the ball
			} else {
				if(dist + _tau >= node.threshold) {         // if there can still be neighbors outside the ball, recursively search right child first
					search(node.right, target, k, heap);
				}

				if (dist - _tau <= node.threshold) {         // if there can still be neighbors inside the ball, recursively search left child
					search(node.left, target, k, heap);
				}
			}
		}
	}
}
