package com.jujutsu.tsne.barneshut;

import java.util.concurrent.atomic.AtomicReference;

class AtomicDouble {
    private AtomicReference<Double> value = new AtomicReference<Double>(Double.valueOf(0.0));
    double addAndGet(double delta) {
        while (true) {
            Double currentValue = value.get();
            Double newValue = Double.valueOf(currentValue.doubleValue() + delta);
            if (value.compareAndSet(currentValue, newValue))
                return currentValue.doubleValue();
        }
    }
    
    double get() {
    	return Double.valueOf(value.get());
    }
}