package com.jujutsu.tsne.demos;

import java.io.File;
import java.util.ArrayList;
import java.util.Date;
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

import com.jujutsu.tsne.FastTSne;
import com.jujutsu.tsne.PrincipalComponentAnalysis;
import com.jujutsu.tsne.SimpleTSne;
import com.jujutsu.tsne.TSne;
import com.jujutsu.utils.MatrixUtils;

public class TSneDemo {
	
	static double perplexity = 20.0;
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
    	double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/resources/datasets/iris_X.txt"), ",");
    	System.out.println("Input is = " + X.length + " x " + X[0].length + " => \n" + ArrayString.printDoubleArray(X));
        PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
    	double [][] Y = pca.pca(X,2);
    	System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + ArrayString.printDoubleArray(Y));
    }

	public static void pca_mnist(int nistSize) {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/resources/datasets/mnist" + nistSize + "_X.txt"));
    	String [] labels = MatrixUtils.simpleReadLines(new File("src/main/resources/datasets/mnist2500_labels.txt"));
    	for (int i = 0; i < labels.length; i++) {
			labels[i] = labels[i].trim().substring(0, 1);
		}
    	System.out.println("Input is = " + X.length + " x " + X[0].length + " => \n" + ArrayString.printDoubleArray(X));
        PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
    	double [][] Y = pca.pca(X,2);
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

	
    public static void tsne_iris() {
    	double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/resources/datasets/iris_X.txt"), ",");
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
        TSne tsne = new SimpleTSne();
		double [][] Y = tsne.tsne(X, 2, initial_dims, perplexity);
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
        run_tsne_mnist(nistSize,new SimpleTSne());
    }
    
    public static void fast_tsne_mnist(int nistSize) {
        run_tsne_mnist(nistSize,new FastTSne());
    }
    
    public static void fast_tsne_mnist(String filename, String labelfilename) {
    	TSne tsne = new FastTSne();
    	System.out.println("Running TSne on " + filename);
        double [][] X = MatrixUtils.simpleRead2DMatrix(new File(filename));
    	String [] labels = MatrixUtils.simpleReadLines(new File(labelfilename));
    	for (int i = 0; i < labels.length; i++) {
			labels[i] = labels[i].trim().substring(0, 1);
		}
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
        System.out.println("Starting TSNE: " + new Date());
        double [][] Y = tsne.tsne(X, 2, initial_dims, perplexity, 10);
        System.out.println("Finished TSNE: " + new Date());
        //System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + ArrayString.printDoubleArray(Y));
        System.out.println("Result is = " + Y.length + " x " + Y[0].length);
        ASCIIFile.write(new File("Java-tsne-result.txt"), ArrayString.printDoubleArray(Y));
        Plot2DPanel plot = new Plot2DPanel();
        
        ColoredScatterPlot setosaPlot = new ColoredScatterPlot("setosa", Y, labels);
        plot.plotCanvas.setNotable(true);
        plot.plotCanvas.setNoteCoords(true);
        plot.plotCanvas.addPlot(setosaPlot);
                
        FrameView plotframe = new FrameView(plot);
        plotframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        plotframe.setVisible(true);
    }

    
    public static void run_tsne_mnist(int nistSize, TSne tsne) {
        System.out.println("Running FAST on " + nistSize + " MNIST digits...");
        double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/resources/datasets/mnist" + nistSize + "_X.txt"));
    	String [] labels = MatrixUtils.simpleReadLines(new File("src/main/resources/datasets/mnist2500_labels.txt"));
    	for (int i = 0; i < labels.length; i++) {
			labels[i] = labels[i].trim().substring(0, 1);
		}
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
        System.out.println("Starting TSNE: " + new Date());
        double [][] Y = tsne.tsne(X, 2, initial_dims, perplexity);
        System.out.println("Finished TSNE: " + new Date());
        //System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + ArrayString.printDoubleArray(Y));
        System.out.println("Result is = " + Y.length + " x " + Y[0].length);
        ASCIIFile.write(new File("Java-tsne-result.txt"), ArrayString.printDoubleArray(Y));
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
        double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/resources/datasets/mnist" + nistSize + "_X.txt"));
    	String [] imgfiles = new String[nistSize];
    	for (int i = 0; i < imgfiles.length; i++) {
			imgfiles[i] = "src/main/resources/nistimgs/img" + i + ".png";
		}
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
        TSne tsne = new SimpleTSne();
        double [][] Y = tsne.tsne(X, 2, initial_dims, perplexity, 1000, true);
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
        if(args.length==0) {
        	System.out.println("usage: For the data format that the TSneDemo accepts, look at the file 'src/main/resources/datasets/minst2500_X.txt' file and accompaning label file 'src/main/resources/datasets/mnist2500_labels.txt'.");
        	System.out.println("       The label file mus have as meny rows as the input matrix");
        	System.out.println("usage: Example using the data and label file in: tsne-demos/src/main/resources/datasets/");
        	System.out.println("usage: java -cp target/tsne-demos-0.0.1-SNAPSHOT.jar com.jujutsu.tsne.demos.TSneDemo minst2500_X.txt mnist2500_labels.txt");
        	System.exit(0);
        }
        //pca_iris();
        //pca_mnist(1000);
        //tsne_iris();
        //tsne_mnist(250);
        //tsne_mnist_icons(500);
        //tsne_mnist(500);
        //tsne_mnist(1000);
        //tsne_mnist(1000);
        //fast_tsne_mnist(2500);
        fast_tsne_mnist(args[0], args[1]);
    }

}
