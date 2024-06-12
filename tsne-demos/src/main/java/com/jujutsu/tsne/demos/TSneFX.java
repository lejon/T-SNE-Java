package com.jujutsu.tsne.demos;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

import com.jujutsu.tsne.FastTSne;
import com.jujutsu.tsne.TSne;
import com.jujutsu.utils.MatrixUtils;

import javafx.application.Application;
import javafx.collections.ObservableList;
import javafx.geometry.Bounds;
import javafx.geometry.Point2D;
import javafx.scene.Scene;
import javafx.scene.chart.Axis;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart;
import javafx.scene.chart.XYChart.Series;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import javafx.stage.Stage;

public class TSneFX extends Application {
	
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
	
	static Color[] palette = {
			Color.BLUE,
			Color.RED,
			Color.GREEN,
			Color.GRAY,
			Color.BLUEVIOLET,
			Color.ORANGE,
			Color.AQUAMARINE,
			Color.INDIGO,
			Color.DARKGOLDENROD
	};
	
	double zoomx, zoomy;
	double x0, x1, y0, y1;
    NumberAxis xaxis = new NumberAxis(), yaxis = new NumberAxis();
    ScatterChart<Number, Number> chart;
    
    public void fast_tsne(Stage stage) {
    	TSne tsne = new FastTSne();
    	int iters = 2000;
    	System.out.println("Running " + iters + " iterations of TSne on " + filename);
    	System.out.println("Reading data");
        double [][] X = MatrixUtils.simpleRead2DMatrix(new File(filename), ",");
        String[] labels;
        Map<Integer, Integer> colors = new HashMap<>();
        if (labelfilename != null) {
        	System.out.println("Reading labels");
        	labels = MatrixUtils.simpleReadLines(new File(labelfilename));
        	for (int i = 0; i < labels.length; i++) {
        		labels[i] = labels[i].trim();
        		if (labels[i].contains(",")) {
        			int comma = labels[i].indexOf(',');
        			colors.put(i,  Integer.parseInt(labels[i].substring(comma+1)));
        			labels[i] = labels[i].substring(0, comma);
        		}
        	}
        } else {
        	labels = new String[X.length];
        	for (int i = 0; i < X.length; i++) {
        		labels[i] = Integer.toString(i);
        	}
        }
        System.out.println("Shape is: " + X.length + " x " + X[0].length);
        System.out.println("Starting TSNE: " + new Date());
        double [][] Y = tsne.tsne(X, 2, initial_dims, perplexity, iters);
        System.out.println("Finished TSNE: " + new Date());
        System.out.println("Result is = " + Y.length + " x " + Y[0].length);
        
        xaxis.setAutoRanging(true);
        yaxis.setAutoRanging(true);
        chart = new ScatterChart<>(xaxis, yaxis);
        
        Series<Number, Number> series = new XYChart.Series<>();
        ObservableList<XYChart.Data<Number, Number>> data = series.getData();
        for (int i = 0; i < Y.length; i++) {
        	XYChart.Data<Number, Number> datum = new XYChart.Data<>(Y[i][0], Y[i][1]);
        	datum.setExtraValue(labels[i]);
        	Text t = new Text(labels[i]);
        	t.setFont(new Font("Helvetica", 200/Math.sqrt(Y.length)));
        	if (colors.get(i) != null) {
        		int ci = Math.abs(colors.get(i) % palette.length);
        		t.setFill(palette[ci]);
        	}
        	
        	datum.setNode(t);
        	data.add(datum);
        }
        final Scene sc = new Scene(chart);
        chart.getData().add(series);
        chart.setOnZoomStarted(ze -> {
        
        	Point2D p = xaxis.sceneToLocal(ze.getSceneX(), ze.getSceneY());
        	
        	x0 = xaxis.getLowerBound();
        	x1 = xaxis.getUpperBound();
        	y0 = yaxis.getLowerBound();
        	y1 = yaxis.getUpperBound();
        	
        	Bounds b = chart.getLayoutBounds();
        	
        	zoomx = x0 + p.getX() / b.getWidth() * (x1 - x0);
        	zoomy = y0 - p.getY() / b.getHeight() * (y1 - y0);
 
        	System.out.println("" + x0 + " " + x1 + " " + y0 + " " + y1 + "p=" + p);
        
        	System.out.println("zoom started at " + zoomx + " " + zoomy + " point " + p);
        	System.out.println("chart size = " + chart.getLayoutBounds());
        });
        chart.setOnZoom(ze -> {
        	double z = 1/ze.getZoomFactor();
        	System.out.println("zoom = " + z);
        	System.out.println(" range was " + x0 + "-" + x1 + ", " + y0 + "-" + y1);
            double xm = zoomx;
            double ym = zoomy;
            double x0n = xm * (1-z) + x0 * z, x1n = xm * (1 - z) + x1 * z;
            double y0n = ym * (1-z) + y0 * z, y1n = ym * (1 - z) + y1 * z;
            double ticks = 1.0;
            while ((x1n-x0n)/ticks > 100) ticks *= 10.0;
            if ((x1n-x0n)/ticks > 20) ticks *= 5.0;
            else if ((x1n - x0n)/ticks > 5) ticks *= 2.0;
//            x0n = ticks * Math.ceil(x0n/ticks);
//            x1n = ticks * Math.floor(x1n/ticks);
//            y0n = ticks * Math.ceil(y0n/ticks);
//            y1n = ticks * Math.floor(y1n/ticks);
            x0 = x0n; x1 = x1n; y0 = y0n; y1 = y1n;
            System.out.println("updating range to " + x0 + "-" + x1 + ", " + y0 + "-" + y1);
            xaxis = new NumberAxis(x0, x1, ticks);
            yaxis = new NumberAxis(y0, y1, ticks);
            ScatterChart<Number, Number> chart2 = new ScatterChart<>(xaxis, yaxis);
            chart2.getData().add(series);
            sc.setRoot(chart2);
            chart2.setOnZoom(chart.getOnZoom());
            chart2.setOnZoomStarted(chart.getOnZoomStarted());
            chart2.setOnZoomFinished(chart.getOnZoomFinished());
            chart = chart2;
        });
        chart.setOnZoomFinished(ze -> { System.out.println("Zoom finished"); });

        stage.setScene(sc);
        stage.show();
    }
    
    
    String filename;
    String labelfilename; // null if no labels
    
    public static void main(String [] args) {
        System.out.println("TSneFX: Runs t-SNE on various datasets.");
        if (args.length<1 || args.length>2) {
        	System.out.println("usage: For the data format TSneFX accepts, look at the file 'src/main/resources/datasets/mnist2500_X.txt' file and accompaning label file 'src/main/resources/datasets/mnist2500_labels.txt'.");
        	System.out.println("       The label file must have as meny rows as the input matrix");
        	System.out.println("usage: Example using the data and label file in: tsne-demos/src/main/resources/datasets/");
        	System.out.println("usage: java -cp target/tsne-demos-X.X.X.jar com.jujutsu.tsne.demos.TSneFX minst2500_X.txt mnist2500_labels.txt");
        	System.out.println("usage: Example using only data file in: tsne-demos/src/main/resources/datasets/");
        	System.out.println("usage: java -cp target/tsne-demos-X.X.X.jar com.jujutsu.tsne.demos.TSneFX minst2500_X.txt");
        	System.exit(0);
        }
        Application.launch(args);
    }
    
	@Override
	public void start(Stage stage) throws Exception {
		Application.Parameters params = getParameters();
		filename = params.getUnnamed().get(0);
		if (params.getUnnamed().size() > 1) {
			labelfilename = params.getUnnamed().get(1);
		}
		fast_tsne(stage);
	}

}
