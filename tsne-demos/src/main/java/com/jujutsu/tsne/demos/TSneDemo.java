package com.jujutsu.tsne.demos;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;

import javax.swing.JFrame;

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
import com.jujutsu.tsne.TSneConfiguration;
import com.jujutsu.utils.MatrixOps;
import com.jujutsu.utils.MatrixUtils;
import com.jujutsu.utils.TSneUtils;

public class TSneDemo {
	
	static double perplexity = 20.0;
	private static int initial_dims = 50;

	public static void saveFile(File file, String text) {
		saveFile(file,text,false);
	}
	
	public static void saveFile(File file, String text, boolean append) {
        try (FileWriter fw = new FileWriter(file, append);
            BufferedWriter bw = new BufferedWriter(fw)) {
            bw.write(text);
            bw.close();
        } catch (IOException e) {
            throw new IllegalArgumentException(e);
        }
    }
	
	public static void pca_iris() {
    	double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/resources/datasets/iris_X.txt"), ",");
    	System.out.println("Input is = " + X.length + " x " + X[0].length + " => \n" + MatrixOps.doubleArrayToPrintString(X));
        PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
    	double [][] Y = pca.pca(X,2);
    	System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + MatrixOps.doubleArrayToPrintString(Y));
    	plotIris(Y);
    }

	public static void pca_mnist(int nistSize) {
		double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/resources/datasets/mnist" + nistSize + "_X.txt"));
    	String [] labels = MatrixUtils.simpleReadLines(new File("src/main/resources/datasets/mnist2500_labels.txt"));
    	for (int i = 0; i < labels.length; i++) {
			labels[i] = labels[i].trim().substring(0, 1);
		}
        PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
    	double [][] Y = pca.pca(X,2);
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
        System.out.println(MatrixOps.doubleArrayToPrintString(X, ", ", 50,10));
        TSne tsne = new SimpleTSne();
		TSneConfiguration config = TSneUtils.buildConfig(X, 2, initial_dims, perplexity, 1000);
		double [][] Y = tsne.tsne(config);
        System.out.println(MatrixOps.doubleArrayToPrintString(Y, ", ", 50,10));
        plotIris(Y);
    }

	static void plotIris(double[][] Y) {
		double [][]        setosa = new double[initial_dims][];
        String []     setosaNames = new String[initial_dims];
        double [][]    versicolor = new double[initial_dims][];
        String [] versicolorNames = new String[initial_dims];
        double [][]     virginica = new double[initial_dims][];
        String []  virginicaNames = new String[initial_dims];
        
        int cnt = 0;
        for (int i = 0; i < initial_dims; i++, cnt++) {
        	setosa[i] = Y[cnt];
            	setosaNames[i] = "setosa";
        }
        for (int i = 0; i < initial_dims; i++, cnt++) {
        	versicolor[i] = Y[cnt];
        	versicolorNames[i] = "versicolor";
        }
        for (int i = 0; i < initial_dims; i++, cnt++) {
        	virginica[i] = Y[cnt];
        	virginicaNames[i] = "virginica";
        }
        
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
    	System.out.println("Running SimpleTSne on " + nistSize + " MNIST digits...");
        run_tsne_mnist(nistSize,new SimpleTSne());
    }
    
    public static void fast_tsne_mnist(int nistSize) {
    	System.out.println("Running FastTSne on " + nistSize + " MNIST digits...");
        run_tsne_mnist(nistSize,new FastTSne());
    }
    
    public static void fast_tsne(String filename, String labelfilename) {
    	TSne tsne = new FastTSne();
    	int iters = 1000;
    	System.out.println("Running " + iters + " iterations of TSne on " + filename);
        double [][] X = MatrixUtils.simpleRead2DMatrix(new File(filename), " ");
    	String [] labels = MatrixUtils.simpleReadLines(new File(labelfilename));
    	for (int i = 0; i < labels.length; i++) {
			labels[i] = labels[i].trim().substring(0, 1);
		}
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
        System.out.println("Starting TSNE: " + new Date());
		TSneConfiguration config = TSneUtils.buildConfig(X, 2, initial_dims, perplexity, iters);
		double [][] Y = tsne.tsne(config);
        System.out.println("Finished TSNE: " + new Date());
        //System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + MatrixOps.doubleArrayToString(Y));
        System.out.println("Result is = " + Y.length + " x " + Y[0].length);
        saveFile(new File("Java-tsne-result.txt"), MatrixOps.doubleArrayToString(Y));
        Plot2DPanel plot = new Plot2DPanel();
        
        ColoredScatterPlot setosaPlot = new ColoredScatterPlot("setosa", Y, labels);
        plot.plotCanvas.setNotable(true);
        plot.plotCanvas.setNoteCoords(true);
        plot.plotCanvas.addPlot(setosaPlot);
                
        FrameView plotframe = new FrameView(plot);
        plotframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        plotframe.setVisible(true);
    }
    
    public static void fast_tsne_no_labels(String filename) {
    	TSne tsne = new FastTSne();
    	int iters = 2000;
    	System.out.println("Running " + iters + " iterations of TSne on " + filename);
        double [][] X = MatrixUtils.simpleRead2DMatrix(new File(filename), ",");
        //System.out.println("X:" + MatrixOps.doubleArrayToString(X));
        X = MatrixOps.log(X, true);
        //System.out.println("X:" + MatrixOps.doubleArrayToString(X));
        X = MatrixOps.centerAndScale(X);
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
        System.out.println("Starting TSNE: " + new Date());
		TSneConfiguration config = TSneUtils.buildConfig(X, 2, initial_dims, perplexity, iters);
		double [][] Y = tsne.tsne(config);
        System.out.println("Finished TSNE: " + new Date());
        //System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + MatrixOps.doubleArrayToString(Y));
        System.out.println("Result is = " + Y.length + " x " + Y[0].length);
        saveFile(new File("Java-tsne-result.txt"), MatrixOps.doubleArrayToString(Y));
        Plot2DPanel plot = new Plot2DPanel();
        
        ScatterPlot setosaPlot = new ScatterPlot("setosa", Color.BLACK, Y);
        plot.plotCanvas.setNotable(true);
        plot.plotCanvas.setNoteCoords(true);
        plot.plotCanvas.addPlot(setosaPlot);
                
        FrameView plotframe = new FrameView(plot);
        plotframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        plotframe.setVisible(true);
    }
    
    public static void pca_no_labels(String filename) {
    	System.out.println("Running PCA on " + filename);
        double [][] X = MatrixUtils.simpleRead2DMatrix(new File(filename), ",");
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
        System.out.println("Starting PCA: " + new Date());
        PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
        double [][] Y = pca.pca(X, 2);
        System.out.println("Finished PCA: " + new Date());
        //System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + MatrixOps.doubleArrayToString(Y));
        System.out.println("Result is = " + Y.length + " x " + Y[0].length);
        saveFile(new File("Java-tsne-result.txt"), MatrixOps.doubleArrayToString(Y));
        Plot2DPanel plot = new Plot2DPanel();
        
        ScatterPlot setosaPlot = new ScatterPlot("setosa", Color.BLACK, Y);
        plot.plotCanvas.setNotable(true);
        plot.plotCanvas.setNoteCoords(true);
        plot.plotCanvas.addPlot(setosaPlot);
                
        FrameView plotframe = new FrameView(plot);
        plotframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        plotframe.setVisible(true);
    }

    
    public static void run_tsne_mnist(int nistSize, TSne tsne) {
        double [][] X = MatrixUtils.simpleRead2DMatrix(new File("src/main/resources/datasets/mnist" + nistSize + "_X.txt"));
    	String [] labels = MatrixUtils.simpleReadLines(new File("src/main/resources/datasets/mnist2500_labels.txt"));
    	for (int i = 0; i < labels.length; i++) {
			labels[i] = labels[i].trim().substring(0, 1);
		}
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
        System.out.println("Starting TSNE: " + new Date());
		TSneConfiguration config = TSneUtils.buildConfig(X, 2, initial_dims, perplexity, 1000);
		double [][] Y = tsne.tsne(config);
        System.out.println("Finished TSNE: " + new Date());
        //System.out.println("Result is = " + Y.length + " x " + Y[0].length + " => \n" + MatrixOps.doubleArrayToString(Y));
        System.out.println("Result is = " + Y.length + " x " + Y[0].length);
        saveFile(new File("Java-tsne-result.txt"), MatrixOps.doubleArrayToString(Y));
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
		TSneConfiguration config = TSneUtils.buildConfig(X, 2, initial_dims, perplexity, 1000);
		double [][] Y = tsne.tsne(config);
        System.out.println(MatrixOps.doubleArrayToPrintString(Y));
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
        if(args.length<1||args.length>2) {
        	System.out.println("usage: For the data format TSneDemo accepts, look at the file 'src/main/resources/datasets/mnist2500_X.txt' file and accompaning label file 'src/main/resources/datasets/mnist2500_labels.txt'.");
        	System.out.println("       The label file must have as meny rows as the input matrix");
        	System.out.println("usage: Example using the data and label file in: tsne-demos/src/main/resources/datasets/");
        	System.out.println("usage: java -cp target/tsne-demos-X.X.X.jar com.jujutsu.tsne.demos.TSneDemo src/main/resources/datasets/mnist2500_X.txt src/main/resources/datasets/mnist2500_labels.txt");
        	System.out.println("usage: Example using only data file in: tsne-demos/src/main/resources/datasets/");
        	System.out.println("usage: java -cp target/tsne-demos-X.X.X.jar com.jujutsu.tsne.demos.TSneDemo src/main/resources/datasets/mnist2500_X.txt");
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
        //pca_no_labels(args[0]);
        //pca_no_labels(args[0]);
        if(args.length==1)
        	fast_tsne_no_labels(args[0]);
        else
        	fast_tsne(args[0], args[1]);
    }

}
