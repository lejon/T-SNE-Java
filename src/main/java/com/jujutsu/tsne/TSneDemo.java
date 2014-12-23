package com.jujutsu.tsne;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;

import org.math.io.files.ASCIIFile;
import org.math.io.parser.ArrayString;
import org.math.plot.FrameView;
import org.math.plot.Plot2DPanel;
import org.math.plot.PlotPanel;
import org.math.plot.plots.ColoredScatterPlot;
import org.math.plot.plots.IconScatterPlot;
import org.math.plot.plots.ScatterPlot;

public class TSneDemo {
	
	static double perplexity = 40.0;
	private static int initial_dims = 50;

	public TSneDemo() {}
	
    public static double[][] nistReadStringDouble(String s, String columnDelimiter, String rowDelimiter) {
        double[][] array;
        String[] rows = s.split(rowDelimiter);
        array = new double[rows.length][];
        for (int i = 0; i < rows.length; i++) {
            List<Double> colvals = new ArrayList<Double>();
            String[] cols = rows[i].split(columnDelimiter);
            for (int j = 0; j < cols.length; j++) {
                if(!(cols[j].length()==0)) {
                    colvals.add(Double.parseDouble(cols[j]));
                }
            }
            array[i] = new double[colvals.size()];
            for (int j = 0; j < colvals.size(); j++) {
                array[i][j] = colvals.get(j);
            }
        }

        return array;
    }
	
    public static double[][] nistReadStringDouble(String s) {
        return nistReadStringDouble(s, " ", "\n");
    }

    public static double[][] nistReadStringDouble(String s, String columnDelimiter) {
        return nistReadStringDouble(s, columnDelimiter, "\n");
    }
	
	public static void pca_iris() {
    	double [][] X = nistReadStringDouble(ASCIIFile.read(new File("src/main/resources/datasets/iris_X.txt")), ",");
    	System.out.println("Input is = " + X.length + " x " + X[0].length + " => \n" + ArrayString.printDoubleArray(X));
    	double [][] Y = TSne.pca(X,2);
    	System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + ArrayString.printDoubleArray(Y));
    }
    
    public static void tsne_iris() {
    	double [][] X = nistReadStringDouble(ASCIIFile.read(new File("src/main/resources/datasets/iris_X.txt")), ",");
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
		double [][] Y = TSne.tsne(X, 2, initial_dims, perplexity);
        System.out.println("Shape is: " + Y.length + " x " + Y[0].length);
        
        double [][] setosa = new double [initial_dims][2];
        String [] setosaNames = new String[initial_dims];
        double [][] versicolor = new double [initial_dims][2];
        String [] versicolorNames = new String[initial_dims];
        double [][] virginica = new double [initial_dims][2];
        String [] virginicaNames = new String[initial_dims];
        
        int cnt = 0;
        for (int i = 0; i < initial_dims; i++, cnt++) {
        	for (int j = 0; j < 2; j++) {
            	setosa[i][j] = Y[cnt][j];
            	setosaNames[i] = "setosa";
			}
        }
        for (int i = 0; i < initial_dims; i++, cnt++) {
        	for (int j = 0; j < 2; j++) {
        		versicolor[i][j] = Y[cnt][j];
        		versicolorNames[i] = "versicolor";
			}
        }
        for (int i = 0; i < initial_dims; i++, cnt++) {
        	for (int j = 0; j < 2; j++) {
        		virginica[i][j] = Y[cnt][j];
        		virginicaNames[i] = "virginica";
			}
        }
        
        System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + ArrayString.printDoubleArray(Y));
        Plot2DPanel plot = new Plot2DPanel();
        
        ScatterPlot setosaPlot = new ScatterPlot("setosa", PlotPanel.COLORLIST[0], setosa);
        setosaPlot.setTags(setosaNames);
        
        ScatterPlot versicolorPlot = new ScatterPlot("versicolor", PlotPanel.COLORLIST[1], versicolor);
        versicolorPlot.setTags(versicolorNames);
        ScatterPlot virginicaPlot = new ScatterPlot("versicolor", PlotPanel.COLORLIST[2], virginica);
        virginicaPlot.setTags(virginicaNames);
        
        plot.plotCanvas.setNotable(true);
        plot.plotCanvas.setNoteCoords(true);
        plot.plotCanvas.addPlot(setosaPlot);
        plot.plotCanvas.addPlot(versicolorPlot);
        plot.plotCanvas.addPlot(virginicaPlot);
        
        //int setosaId = plot.addScatterPlot("setosa", setosa);
        //int versicolorId = plot.addScatterPlot("versicolor", versicolor);
        //int virginicaId = plot.addScatterPlot("virginica", virginica);
        
        FrameView plotframe = new FrameView(plot);
        plotframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        plotframe.setVisible(true);
    }
    
    public static void tsne_mnist(int nistSize) {
        System.out.println("Running example on " + nistSize + " MNIST digits...");
        double [][] X = nistReadStringDouble(ASCIIFile.read(new File("src/main/resources/datasets/mnist" + nistSize + "_X.txt")));
    	String [] labels = new ASCIIFile(new File("src/main/resources/datasets/mnist2500_labels.txt")).readLines();
    	for (int i = 0; i < labels.length; i++) {
			labels[i] = labels[i].trim().substring(0, 1);
		}
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
        double [][] Y = TSne.tsne(X, 2, initial_dims, perplexity);
        System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + ArrayString.printDoubleArray(Y));
        Plot2DPanel plot = new Plot2DPanel();
        
        ColoredScatterPlot setosaPlot = new ColoredScatterPlot("setosa", Y, labels);
        plot.plotCanvas.setNotable(true);
        plot.plotCanvas.setNoteCoords(true);
        plot.plotCanvas.addPlot(setosaPlot);
                
        FrameView plotframe = new FrameView(plot);
        plotframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        plotframe.setVisible(true);
    }
    
    public static void tsne_mnist_icons(int nistSize) {
        System.out.println("Running example on " + nistSize + " MNIST digits...");
        double [][] X = nistReadStringDouble(ASCIIFile.read(new File("src/main/resources/datasets/mnist" + nistSize + "_X.txt")));
    	String [] imgfiles = new String[nistSize];
    	for (int i = 0; i < imgfiles.length; i++) {
			imgfiles[i] = "imgs/img" + i + ".png";
		}
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
        double [][] Y = TSne.tsne(X, 2, initial_dims, perplexity, 2000, false);
        System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + ArrayString.printDoubleArray(Y));
        Plot2DPanel plot = new Plot2DPanel();
        
        IconScatterPlot setosaPlot = new IconScatterPlot("setosa", Y, imgfiles);
        plot.plotCanvas.setNotable(true);
        plot.plotCanvas.setNoteCoords(true);
        plot.plotCanvas.addPlot(setosaPlot);
                
        FrameView plotframe = new FrameView(plot);
        plotframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        plotframe.setVisible(true);
    }
    
    public static void main(String [] args) {
        System.out.println("TSneDemo: Runs t-SNE on various dataset.");
        //pca_iris();
        //tsne_iris();
        //tsne_mnist(250);
        //tsne_mnist_icons(500);
        tsne_mnist(500);
        //tsne_mnist(1000);
        //tsne_mnist(2500);
    }

}
