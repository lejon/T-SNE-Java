package com.jujutsu.tsne.barneshut;

import static org.junit.Assert.assertEquals;

import java.util.concurrent.ThreadLocalRandom;

import org.junit.Test;

public class ParallelBHUtilsTest {

	@Test
	public void testAtomicIntAddition() {
		double d1 = ThreadLocalRandom.current().nextDouble();
		double d2 = ThreadLocalRandom.current().nextDouble();
		
		AtomicDouble al = new AtomicDouble();
		al.addAndGet(d1);
		al.addAndGet(d2);
		
		assertEquals(d1+d2,al.get(),0.000000001);
	}
	
}
