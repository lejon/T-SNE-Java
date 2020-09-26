package com.jujutsu.tsne.barneshut;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ParallelVpTree<StorageType> extends VpTree<StorageType> {
	
	public ParallelVpTree(Distance distance) {
		super(distance);
	}
	
	public List<TreeSearchResult> searchMultiple(ParallelVpTree<StorageType> tree, DataPoint [] targets, int k) {
		VpTree<StorageType>.Node node = tree.getRoot();
		
		List<TreeSearchResult> results =  
		IntStream.range(0, targets.length).parallel().mapToObj(n -> {
			DataPoint target = targets[n];
			List<DataPoint> indices = new ArrayList<>();
			List<Double> distances = new ArrayList<>();
			PriorityQueue<HeapItem> heap = new PriorityQueue<HeapItem>(k,new Comparator<HeapItem>() {
				@Override
				public int compare(HeapItem o1, HeapItem o2) {
					return -1 * o1.compareTo(o2);
				}
			}); 

			double tau = Double.MAX_VALUE;
			// Perform the search
			node.search(node, target, k, heap, tau);

			// Gather final results
			while(!heap.isEmpty()) {
				indices.add(_items[heap.peek().index]);
				distances.add(heap.peek().dist);
				heap.remove();
			}
			
			// Results are in reverse order 
			Collections.reverse(indices);
			Collections.reverse(distances);

			return new TreeSearchResult(indices, distances,n);
		}).collect(Collectors.toList());
		
		return results;
	}

}
