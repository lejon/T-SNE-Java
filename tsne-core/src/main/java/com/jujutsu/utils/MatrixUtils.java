package com.jujutsu.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MatrixUtils {
	public static double[][] simpleRead2DMatrix(File file) {
		return simpleRead2DMatrix(file, " ");
	}
	
	public static double[][] simpleRead2DMatrix(File file, String columnDelimiter) {
		List<double[]> rows = new ArrayList<>();
		
        try (FileReader fr = new FileReader(file)) {
            BufferedReader b = new BufferedReader(fr);
            String line;
            while ((line = b.readLine()) != null && !line.matches("\\s*")) {
                String[] cols = line.trim().split(columnDelimiter);
                double [] row = new double[cols.length];
                for (int j = 0; j < cols.length; j++) {
                	if(!(cols[j].length()==0)) {
                		row[j] = Double.parseDouble(cols[j].trim());
                    }
                }
                rows.add(row);
            }
            b.close();
        } catch (IOException e) {
            throw new IllegalArgumentException(e);
        }
        
        double[][] array = new double[rows.size()][];
        int currentRow = 0;
        for (double[] ds : rows) {
			array[currentRow++] = ds;
		}
                
        return array;
    }
	
	public static String[] simpleReadLines(File file) {
		List<String> rows = new ArrayList<>();
		
        try (FileReader fr = new FileReader(file);
        	BufferedReader b = new BufferedReader(fr)) {
        	String line;
        	while ((line = b.readLine()) != null) {
        		rows.add(line.trim());
        	}
        } catch (IOException e) {
            throw new IllegalArgumentException(e);
        }
        
        String[] lines = new String[rows.size()];
        int currentRow = 0;
        for (String line : rows) {
        	lines[currentRow++] = line;
		}
                
        return lines;
    }
}
