package com.jujutsu.tsne.demos;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;

import org.math.io.files.ASCIIFile;
import org.math.plot.FrameView;
import org.math.plot.Plot2DPanel;
import org.math.plot.PlotPanel;
import org.math.plot.plots.ColoredScatterPlot;
import org.math.plot.plots.ScatterPlot;

import com.jujutsu.tsne.SimpleTSne;
import com.jujutsu.tsne.TSne;
import com.jujutsu.utils.MatrixOps;
import com.jujutsu.utils.MatrixUtils;

public class TSneASCIIDemo {

	static double perplexity = 20.0;
	private static int initial_dims = 50;

	public TSneASCIIDemo() {}
	
	public static double[][] stringToMatrix(String s, String columnDelimiter, String rowDelimiter) {
        String[] rows = s.split(rowDelimiter);
        List<double []> mRows = new ArrayList<>();
        row:
        for (int i = 0; i < rows.length; i++) {
            List<Double> colvals = new ArrayList<Double>();
            String[] cols = rows[i].split(columnDelimiter);
            for (int j = 0; j < cols.length; j++) {
                if(!(cols[j].length()==0)) {
                	if(cols[j].equals("NA")) {
                		continue row;
                	} else {
                		colvals.add(Double.parseDouble(cols[j]));
                	}
                }
            }
            double [] mRow = new double[colvals.size()];
            for (int j = 0; j < colvals.size(); j++) {
            	mRow[j] = colvals.get(j);
            }
            mRows.add(mRow);
        }
        
        double[][] array = new double[mRows.size()][];
        for (int i = 0; i < array.length; i++) {
			array[i] = mRows.get(i);
		}

        return array;
    }
	
    public static double[][] nistReadStringDouble(String s) {
        return stringToMatrix(s, " ", "\n");
    }

    public static double[][] nistReadStringDouble(String s, String columnDelimiter) {
        return stringToMatrix(s, columnDelimiter, "\n");
    }

    public static double [][] runTSne(double [][] matrix) {
    	System.out.println("Shape is: " + matrix.length + " x " + matrix[0].length);
    	TSne tsne = new SimpleTSne();
    	double [][] Y = tsne.tsne(matrix, 2, initial_dims, perplexity);
    	System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + MatrixOps.doubleArrayToPrintString(Y));
    	return Y;
    }

	static void plot2D(String[] labels, double[][] Y) {
		Plot2DPanel plot = new Plot2DPanel();
    	if(labels != null) {
    		for (int i = 0; i < labels.length; i++) {
    			labels[i] = labels[i].trim().substring(0, 1);
    		}
    		ColoredScatterPlot setosaPlot = new ColoredScatterPlot("setosa", Y, labels);
    		plot.plotCanvas.addPlot(setosaPlot);
    	} else {
    		ScatterPlot dataPlot = new ScatterPlot("Data", PlotPanel.COLORLIST[0], Y);
    		plot.plotCanvas.addPlot(dataPlot);

    	}
    	plot.plotCanvas.setNotable(true);
    	plot.plotCanvas.setNoteCoords(true);
    	
    	FrameView plotframe = new FrameView(plot);
    	plotframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	plotframe.setVisible(true);
	}

	
	public static LoadResult loadData(String [] args) {
		String usage = "Usage: TSneASCIIDemo [options] datafile \nAvailable options:\n	-labelFile filename <default = 'null'> -colDelim <default = ','>\n";
		if(args.length<1) System.out.println("No arguments given.\n" + usage);
		File file = null;
		File labelFile = null;
		String colDelim = ",";
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-colDelim")) {
				colDelim = args[i+1].trim(); 
				if(colDelim.equals("space")) colDelim = " "; 
				i++;
			} else if (args[i].equals("-labelFile")) {
				labelFile = new File(args[i+1].trim());
				i++;
			} else {
				file = new File(args[i].trim());
				if (!file.exists()) {
					System.out.println("File " + file + " doesn't exists.\n" + usage);
					return null;
				} 
				i++;
			}
		}
		
		String [] labels = null;
		if(labelFile!=null) labels = MatrixUtils.simpleReadLines(labelFile);
		LoadResult lr = null;
		if(file!=null) {
			System.out.println("Running T-SNE on " + file.getAbsolutePath() + " ...");
			lr = new LoadResult(MatrixUtils.simpleRead2DMatrix(file,colDelim),labels);
		}
		return lr;
	}
	
	static void printDataSample(int printSampleSize, double[][] matrix) {
		if(matrix!=null) {
			System.out.println("Loaded ("+printSampleSize+" samples): ");
			for (int i = 0; i < Math.min(printSampleSize, matrix.length); i++) {
				for (int j = 0; j < matrix[0].length; j++) {
					System.out.print(matrix[i][j] + ",");
				}
				System.out.println();
			}
			System.out.println("...");
		}
	}
	
	static class LoadResult {
		public LoadResult(double[][] matrix, String[] labels) {
			super();
			this.matrix = matrix;
			this.labels = labels;
		}
		public double [][] matrix;
		public String [] labels;
	}

	public static void main(String[] args) {
		LoadResult lr = loadData(args);
		printDataSample(10,lr.matrix);
		if(lr!=null)  {
			double [][] Y = runTSne(lr.matrix);
			plot2D(lr.labels, Y);
			File out = new File("tsne-result.txt");
			ASCIIFile.writeDoubleArray(out, Y);
		}
	}
}
