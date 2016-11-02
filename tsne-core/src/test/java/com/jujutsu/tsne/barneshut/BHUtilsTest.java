package com.jujutsu.tsne.barneshut;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

import org.junit.Test;

import com.jujutsu.tsne.barneshut.VpTree.HeapItem;

public class BHUtilsTest {

	@Test
	public void testHeapItemOrdering() {
		DataPoint [] _items = new DataPoint [3];
		double [] p1 = {1.0,2.0,3.0};
		_items[0] = new DataPoint(3,0,p1);
		double [] p2 = {4.0,5.0,6.0};
		_items[1] = new DataPoint(3,1,p2);
		double [] p3 = {7.0,8.0,9.0};
		_items[2] = new DataPoint(3,2,p3);
		
		//System.out.println("Items=" + Arrays.toString(_items));
		int k = 10;
		PriorityQueue<HeapItem> heap = new PriorityQueue<HeapItem>(k,new Comparator<HeapItem>() {
			@Override
			public int compare(HeapItem o1, HeapItem o2) {
				return -1 * o1.compareTo(o2);
			}
		}); 
		heap.add(new HeapItem(0, 7.2));
		heap.add(new HeapItem(1, 4.1));
		heap.add(new HeapItem(2, 5.1));

		//System.out.println("Heap Is:");
		//System.out.println(heap);
        //System.out.println();
        
        // // Gather final results
		List<DataPoint> results  = new ArrayList<DataPoint>();
		List<Double> distances = new ArrayList<Double>();

        while(!heap.isEmpty()) {
            results.add(_items[heap.peek().index]);
            distances.add(heap.peek().dist);
            heap.remove();
        }
        
        
        Collections.reverse(results);
        DataPoint [] expectedRes = { _items[1], _items[2], _items[0] };
        DataPoint [] actualRes = new DataPoint[expectedRes.length];
        for (int i = 0; i < actualRes.length; i++) {
        	actualRes[i] = results.get(i);
		}

        assertArrayEquals(expectedRes, actualRes);
        
        Collections.reverse(distances);
        double [] expected = {4.1,5.1,7.2};
        double [] actual = new double[expected.length];
        for (int i = 0; i < actual.length; i++) {
			actual[i] = distances.get(i);
		}
        assertArrayEquals(expected, actual, 0.000000001);
	}
	
	@Test
	public void testDistanceComparator() {
		DataPoint [] _items = new DataPoint [3];
		double [] p1 = {1.0,2.0,3.0};
		_items[0] = new DataPoint(3,0,p1);
		double [] p2 = {4.0,5.0,6.0};
		_items[1] = new DataPoint(3,1,p2);
		double [] p3 = {10.0,11.0,12.0};
		_items[2] = new DataPoint(3,2,p3);
		DataPoint [] expectedRes = { _items[1], _items[0], _items[2] };
		Arrays.sort(_items, new DistanceComparator(_items[1]));
		assertArrayEquals(expectedRes, _items);
	}
	
	@Test
	public void testEuclidDistance() {
		double [] p1 = {1.0,2.0,3.0};
		DataPoint d1 = new DataPoint(3,0,p1);
		double [] p2 = {4.0,5.0,6.0};
		DataPoint d2 = new DataPoint(3,1,p2);
		assertEquals(5.196152,d1.euclidean_distance(d2),0.000001);
	}
	
	@Test
	public void testNthElementSecond() {
		DataPoint [] _items = new DataPoint [6];
		int idx = 0;
		double [] p6 = {1000.0,2000.0,3000.0};
		_items[idx++] = new DataPoint(3,5,p6);
		double [] p1 = {1.0,2.0,3.0};
		_items[idx++] = new DataPoint(3,0,p1);
		double [] p3 = {10.0,11.0,12.0};
		_items[idx++] = new DataPoint(3,2,p3);
		double [] p4 = {20.0,40.0,60.0};
		_items[idx++] = new DataPoint(3,3,p4);
		double [] p5 = {100.0,200.0,300.0};
		_items[idx++] = new DataPoint(3,4,p5);
		double [] p2 = {4.0,5.0,6.0};
		_items[idx++] = new DataPoint(3,1,p2);

		VpTree.nth_element(_items, 0, 1, 2, new DistanceComparator(_items[1]));
		//System.out.println("Items:" + Arrays.toString(_items));
	}
	
	@Test
	public void testNthElementIntMiddle() {
		int [] array = {5, 6, 4, 3, 2, 6, 7, 9, 3};
		int pivot = 4;
		VpTree.nth_element(array, 0, pivot, 9);
		//System.out.println("Items:" + Arrays.toString(array));
		assertEquals(5,array[pivot]);
	}

	@Test
	public void testNthElementInt1() {
		int [] array = {5, 6, 4, 3, 2, 6, 7, 9, 3};
		int pivot = 1;
		VpTree.nth_element(array, 0, pivot, 7);
		//System.out.println("Items:" + Arrays.toString(array));
		assertEquals(3,array[pivot]);
	}
	
	@Test
	public void testNthElement1() {
		int [] array = {5, 6, 4, 3, 2, 6, 7, 9, 3};
		int pivot = 1;
		VpTree.nth_element(array, 0, pivot, 2);
		//System.out.println("Items:" + Arrays.toString(array));
		assertEquals(6,array[pivot]);
	}
	
	@Test
	public void testNthElementFull() {
		int [] array = {5, 6, 4, 3, 2, 6, 7, 9, 3};
		int pivot = 4;
		VpTree.nth_element(array, 0, pivot, 9);
		//System.out.println("Items:" + Arrays.toString(array));
		assertEquals(5,array[pivot]);
	}
	
	@Test
	public void testNthElementInt2() {
		int [] array = {5, 6, 4, 3, 2, 6, 7, 9, 3};
		int pivot = 2;
		VpTree.nth_element(array, 0, pivot, 3);
		//System.out.println("Items:" + Arrays.toString(array));
		assertEquals(6,array[pivot]);
	}

	@Test
	public void testNthElementInt3() {
		int [] array = {5, 6, 4, 3, 2, 6, 7, 9, 3};
		int pivot = 6;
		VpTree.nth_element(array, 0, pivot, 8);
		//System.out.println("Items:" + Arrays.toString(array));
		assertEquals(7,array[pivot]);
	}
	
	@Test
	public void testNthElementInt4() {
		int [] array = {5, 6, 4, 3, 2, 6, 7, 9, 3};
		int pivot = 6;
		VpTree.nth_element(array, 5, pivot, 8);
		//System.out.println("Items:" + Arrays.toString(array));
		assertEquals(7,array[pivot]);
	}

	@Test
	public void testPrioHeap() {
		int k = 10;
		PriorityQueue<Integer> heap = new PriorityQueue<Integer>(k,Collections.reverseOrder());
		int [] array = {5, 6, 4, 3, 2, 6, 7, 9, 3};
		
		for (int i = 0; i < array.length; i++) {
			heap.add(array[i]);			
		}
//		System.out.println(heap);
		int cnt = 0;		
		int [] result = new int[array.length]; 
		while(!heap.isEmpty()) {
			result[cnt++] = heap.remove();
//			System.out.print(heap.remove() + ", ");
		}
		int [] expected = {9, 7, 6, 6, 5, 4, 3, 3, 2};
		assertArrayEquals(expected, result);
	}
	
}
