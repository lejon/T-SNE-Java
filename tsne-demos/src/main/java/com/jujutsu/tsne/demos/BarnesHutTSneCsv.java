package com.jujutsu.tsne.demos;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;

import joinery.DataFrame;
import joinery.DataFrame.NumberDefault;
import joinery.DataFrame.Predicate;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;
import org.math.plot.FrameView;
import org.math.plot.Plot2DPanel;
import org.math.plot.PlotPanel;
import org.math.plot.plots.ColoredScatterPlot;
import org.math.plot.plots.ScatterPlot;

import com.jujutsu.tsne.barneshut.BHTSne;
import com.jujutsu.tsne.barneshut.BarnesHutTSne;
import com.jujutsu.tsne.barneshut.ParallelBHTsne;
import com.jujutsu.utils.MatrixOps;
import com.jujutsu.utils.MatrixUtils;

public class BarnesHutTSneCsv {
	static int     initial_dims    = -1;
	static int     output_dims     = 2;
	static double  perplexity      = 20.0;
	static double  theta           = 0.5;
	static boolean hasLabels       = true;
	static boolean scale_log       = false;
	static boolean normalize       = false;
	static boolean addNoise        = false;
	static boolean subSample       = false;
	static boolean hasHeader       = true;
	static boolean doPlot          = true;
	static boolean doSave          = false;
	static boolean transpose_after = false;
	static String  output_fn       = null;
	static String  label_fn        = null;
	static String  naString        = null;
	static int     iterations      = 2000;
	static int     label_col_no    = 0;
	static String  label_col_name  = null;
	static boolean parallel        = false;

	public static DataFrame<Object> processCommandline(String [] args) throws IOException {
		CommandLineParser parser = new PosixParser();

		Options options = new Options();
		options.addOption( "perp", "perplexity",    true, 
				"set the perplexity of the t-SNE algorithm (default " + perplexity +")" );
		options.addOption( "th", "theta",           true, 
				"sets the theta value of the Barnes Hut algorithm (default " + theta +")" );
		options.addOption( "idims", "initial_dims", true, 
				"scale the dataset to initial dims with PCA before running t-SNE (default " + initial_dims + "). Negative number indicates no scaling"  );
		options.addOption( "nolbls", "no_labels",   false, 
				"The dataset does not contain any labels (if not set labels are assumed to be in the first column)" );
		options.addOption( "lblf", "label_file",    true, 
				"Separate input file with dataset labels, one label per row, must contain at least as many rows as in the dataset. Extra labels will be thrown away" );
		options.addOption( "na", "na_string",   true, 
				"Dataset can contain N/A, the given string is parsed as N/A in the dataset" );
		options.addOption( "nohdr", "no_headers",   false, 
				"If set, won't try to read a first row of column headers / names" );
		options.addOption( "log", "scale_log",      false, 
				"Scale the dataset by first taking the log of each datapoint (keeping zeros) " );
		options.addOption( "norm", "normalize",     false, 
				"Normalize the data by subtracting the mean and dividing by the stdev (this is done after eventual log) " );
		options.addOption( "iter", "iterations",    true, 
				"How many iterations to run, default is " + iterations);
		options.addOption( "noplt", "no_plot",           false, 
				"Don't plot the resulting dataset " );
		options.addOption( "lblcolno", "label_column_no",    true, 
				"If labels are not in first column, this option gives the index of the label column" );
		options.addOption( "lblcolnme", "label_column_name",    true, 
				"If labels are not in first column, this option gives the name of the label column. Requires headers in the dataset" );
		options.addOption( "shw", "show",           false, 
				"Show displays the tabular data of a data frame in a gui window " );
		options.addOption( "dn",  "drop_name",      true, 
				"drop column names. Takes a list of names (Example: \"Customer Name,Comment,Id\") representing the cloumn names to drop. This is done AFTER any drop_column!" );
		options.addOption( "dc",  "drop_column",    true, 
				"drop column no's. Takes a list of integers (Example: \"1,2,8,11\") representing the cloumns to drop" );
		options.addOption( "sep", "separator",      true, 
				"column separator ',' , ';' , '\\t' (',' per default). '\\t' denotes tab" );
		options.addOption( "dbl", "double_default", false, 
				"use Double as number format (Long is default but even with Long default, numbers with decimals will still be parsed as Double)" );
		options.addOption( "trsp", "transpose", false, 
				"transpose the dataset first" );
		options.addOption( "trspa", "transpose_after", false, 
				"transpose the dataset after t-SNE is done" );
		options.addOption( "out", "output_file",    true, 
				"Save the result to the given filename" );
		options.addOption( "no", "noise",    false, 
				"add a small amount of noise to each column. This can be useful with highly structured datasets which can otherwise cause problems" );
		options.addOption( "ss", "subsample",    false, 
				"the current implementation does not handle very large datasets due to memory and time constraints. Adding this flag will uniformly subsample the dataset" );
		options.addOption( "pa", "parallel",     false, 
				"Run parts of algorithm in parallel. Using this option will hog your CPUs!" );
		options.addOption( "odim", "output_dims",    true, 
				"Alternatives are '2D' or '3D' default is (" + output_dims + "D " );


		CommandLine parsedCommandLine = null;
		HelpFormatter formatter = new HelpFormatter();

		// Try to parse the command line
		try {
			parsedCommandLine = parser.parse( options, args );
		} catch (org.apache.commons.cli.ParseException e) {
			System.out.println("TSneCsv: Could not parse command line due to:  " + e.getMessage());
			System.out.println("Args where:");
			for (int i = 0; i < args.length; i++) {
				System.out.print(args[i] + ", ");
			}
			formatter.printHelp("TSneCsv [options] <csv file>", options );
			System.exit(-1);
		}

		// Make sure we got a CSV file
		if(parsedCommandLine.getArgs().length==0) {
			System.out.println("No CSV file given...");
			formatter.printHelp("TSneCsv [options] <csv file>", options );
			System.exit(255);
		}

		final List<DataFrame<Object>> frames = new ArrayList<>();

		String sep = parsedCommandLine.hasOption( "separator" ) ? (String) parsedCommandLine.getOptionValue("separator").trim()	: ",";

		if(!(sep.equals(",") || sep.equals(";") || sep.equals("\\t"))) {
			System.out.println("Only the separators ',' , ';' or '\\t' is currently supported...");
			formatter.printHelp( DataFrame.class.getCanonicalName() + " [options] <csv file>", options );
			System.exit(255);
		}

		DataFrame<Object> df = null;

		// Read DataFrame with any given separator
		String csvFilename = parsedCommandLine.getArgs()[0];
		
		if (parsedCommandLine.hasOption( "iterations" )) {
			iterations = Integer.parseInt(parsedCommandLine.getOptionValue("iterations").trim());
		}
		
		if (parsedCommandLine.hasOption( "label_column_no" )) {
			label_col_no = Integer.parseInt(parsedCommandLine.getOptionValue("label_column_no").trim());
			label_col_no--; // Label col is 1 indexed
		}

		if (parsedCommandLine.hasOption( "label_column_name" )) {
			label_col_name = parsedCommandLine.getOptionValue( "label_column_name" );
		}
		
		if (parsedCommandLine.hasOption( "no_headers" )) {
			hasHeader = false;
		}
		
		if (parsedCommandLine.hasOption( "na_string" )) {
			naString = parsedCommandLine.getOptionValue("na_string").trim();
		}

		System.out.println("TSneCsv: Running " + iterations + " iterations of t-SNE on " + csvFilename);
		System.out.println("NA string is: " + naString) ;

		if(parsedCommandLine.hasOption( "double_default" ) ) {
			df = DataFrame.readCsv(csvFilename,sep,NumberDefault.DOUBLE_DEFAULT, naString, hasHeader);
			frames.add(df);
		} else {
			df = DataFrame.readCsv(csvFilename,sep,NumberDefault.LONG_DEFAULT, naString, hasHeader);
			frames.add(df);
		}
		
		System.out.println("Loaded CSV with: " + df.length()+ " rows and " + df.size() +" columns.");

		if (parsedCommandLine.hasOption( "transpose" )) {
			df = df.transpose();
		}
		if (parsedCommandLine.hasOption( "perplexity" )) {
			perplexity = Double.parseDouble(parsedCommandLine.getOptionValue("perplexity").trim());
		}
		if (parsedCommandLine.hasOption( "theta" )) {
			theta = Double.parseDouble(parsedCommandLine.getOptionValue("theta").trim());
		}
		if (parsedCommandLine.hasOption( "initial_dims" )) {
			initial_dims = Integer.parseInt(parsedCommandLine.getOptionValue("initial_dims").trim());
		}
		if (parsedCommandLine.hasOption( "no_labels" )) {
			hasLabels = false;
			//if(initial_dims>df.size()) 
			//	throw new IllegalArgumentException("Requested more initial dims " + initial_dims + "than available in the dataset: " + df.size());
		}
		if (parsedCommandLine.hasOption( "label_file" )) {
			label_fn = parsedCommandLine.getOptionValue("label_file").trim();
		}
		if (parsedCommandLine.hasOption( "scale_log" )) {
			System.out.println("Log transforming dataset...");
			scale_log = true;
		}
		if (parsedCommandLine.hasOption( "normalize" )) {
			System.out.println("Normalizing dataset...");
			normalize = true;
		}
		if (parsedCommandLine.hasOption( "noise" )) {
			System.out.println("Adding noise...");
			addNoise = true;
		}

		if (parsedCommandLine.hasOption( "subsample" )) {
			System.out.println("Subsampling dataset...");
			subSample = true;
		}
		
		if (parsedCommandLine.hasOption( "parallel" )) {
			System.out.println("Using parallel Barnes Hut t-SNE...");
			parallel = true;
		}
		
		if (parsedCommandLine.hasOption( "output_dims" )) {
			String dimSpec = parsedCommandLine.getOptionValue("output_dims").trim();
			if(dimSpec.equalsIgnoreCase("2d")) {
				output_dims = 2;
			} else if(dimSpec.equalsIgnoreCase("3d")) {
				output_dims = 3;
			} else {
				throw new IllegalArgumentException("Only legal output_dims options are '2D' or '3D'.");
			}
		}

		// Now process any drop column arguments, it makes sense to
		// drop the integer indexed first since named ones are 
		// independent of index and thus we can do those after
		// and get the drops correct in case of both index and
		// named drops
		if (parsedCommandLine.hasOption( "drop_column" )) {
			String [] colIdxs = parsedCommandLine.getOptionValue("drop_column").split(",");
			List<Integer> dropIdxs = new ArrayList<>();
			for (int i = 0; i < colIdxs.length; i++) {
				String colVal = colIdxs[i].trim();
				if(colVal.length()>0) dropIdxs.add(Integer.parseInt(colVal));
			}
			Integer [] idxs = dropIdxs.toArray(new Integer[0]);
			if(dropIdxs.size()>0) df = df.drop(idxs);
		}

		if (parsedCommandLine.hasOption( "drop_name" )) {
			String [] colNames = parsedCommandLine.getOptionValue("drop_name").split(",");
			for (String colName : colNames) {
				String trimmedName = colName.trim();
				if(trimmedName.length()>0) df = df.drop(trimmedName);				
			}
		}

		if (parsedCommandLine.hasOption( "no_plot" )) {
			doPlot = false;
		}

		if (parsedCommandLine.hasOption( "output_file" )) {
			output_fn = parsedCommandLine.getOptionValue("output_file");
			doSave = true;
		}
		
		final double noRows = df.length();
		if(subSample) {
			df = df.select(new Predicate<Object>() {
				double prob = Math.min(1.0, 2500.0 / noRows);
				Random rnd = new Random();
				@Override
				public Boolean apply(List<Object> values) {
					return rnd.nextDouble() < prob;
				}});
		}

		if (parsedCommandLine.hasOption( "show" )) {
			if (frames.size() == 1) {
				frames.get(0).show();
				return df;
			}
		}

		return df;
	}

	static void tsne_csv(String [] args) throws IOException {
		DataFrame<Object> df = processCommandline(args);

		System.out.println("Dataset types:" + df.types());
		System.out.println(df.head(10));

		String [] labels = null;
		// IF the labels are in the DataFrame extract them as numbers to 'ys' 
		// and strings to 'plotLabels' and remove that column from the DataFrame
		// Also build a mapping from numeric <-> string labels
		if(hasLabels) {
			labels = new String[df.length()];
			int idx = 0;
			List<Object> col = null;
			if(label_col_name!=null){
				System.out.println("Using labels from colum name: " + label_col_name);
				col = df.col(label_col_name);
			} else {
				System.out.println("Using labels from colum index: " + (label_col_no+1) );
				col = df.col(label_col_no);
			}
			for (Object lbl : col) {
				labels[idx++] = lbl.toString();
			}

			// Now drop the label column
			if(label_col_name!=null){
				df = df.drop(label_col_name);
			} else {
				df = df.drop(label_col_no);
			}
		}

		// Load labels from external file if they are not in the data set
		if(!hasLabels && label_fn!=null) {
			System.out.println("Loading labels from:" + label_fn);
			String [] readLabels = MatrixUtils.simpleReadLines(new File(label_fn));

			// Got too many labels, remove the extra ones
			if(readLabels.length > df.length()) {
				labels = Arrays.copyOf(readLabels, df.length());
			} else {
				labels = readLabels;
			}
		}
				
		int data_len = transpose_after ? df.size() : df.length();
		if(labels != null && labels.length < data_len) {
			throw new IllegalArgumentException("The number of labels (" + labels.length + ") is not the same (or more) as the number of rows in the dataset (" + data_len + ").");
		}

		final DataFrame<Double> ddf = df.cast(Double.class);
		
		//double[][] matrix = ddf.fillna(0.0).toArray(double[][].class);
		double[][] matrix = ddf.toModelMatrix(0.0);
		if(scale_log) matrix = MatrixOps.log(matrix, true);
		if(normalize) matrix = MatrixOps.centerAndScale(matrix);
		if(addNoise)  matrix = MatrixOps.addNoise(matrix);
		System.out.println(MatrixOps.doubleArrayToPrintString(matrix,5,5,20));

		BarnesHutTSne tsne;
		long t1 = System.currentTimeMillis();
		if(parallel) {			
			tsne = new ParallelBHTsne();
		} else {
			tsne = new BHTSne();
		}
		double [][] Y = tsne.tsne(matrix, output_dims, initial_dims, perplexity, iterations, true, theta);
		if(transpose_after) Y = MatrixOps.transposeSerial(matrix);
		long t2 = System.currentTimeMillis();
		System.out.println("TSne took: " + ((double) (t2-t1) / 1000.0) + " seconds");
		
		if(doSave) {
			DataFrame<Object> outdf;
			if(labels!=null) {
				String [] colnames = new String[output_dims+1];
					colnames[0] = "label";
					colnames[1] = "X";
					colnames[2] = "Y";
				if(output_dims==3) {
					colnames[3] = "Z";
				}
				outdf = new DataFrame<Object>(colnames);
				int rowIdx = 0;
				for(double [] row : Y) {
					if(output_dims==2) {
						outdf.append(Arrays.asList(labels[rowIdx++], row[0],row[1]));
					} else {
						outdf.append(Arrays.asList(labels[rowIdx++], row[0],row[1],row[2]));
					}
				}
				System.out.println(outdf);
				outdf.writeCsv(output_fn);				
			} else {
				String [] colnames = new String[output_dims];
				colnames[0] = "X";
				colnames[1] = "Y";
				if(output_dims==3) {
					colnames[2] = "Z";
				}
				outdf = new DataFrame<Object>(colnames);
				for(double [] row : Y) {
					if(output_dims==2) {
						outdf.append(Arrays.asList(row[0],row[1]));
					} else {
						outdf.append(Arrays.asList(row[0],row[1],row[2]));
					}
				}
				System.out.println(outdf);
				outdf.writeCsv(output_fn);
			}
		}
		
		if(doPlot && output_dims==2)
			plot2D(labels, Y);
	}

	static void plot2D(String[] labels, double[][] Y) {
		Plot2DPanel plot = new Plot2DPanel();
		if(labels != null) {
			ColoredScatterPlot setosaPlot = new ColoredScatterPlot("TSne Result", Y, labels);
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
	
	/*static void plot2DJChart(String[] labels, double[][] Y) {
		JFreeChartColorScatterPlot plt = new JFreeChartColorScatterPlot("TSne REsult", Y, labels);
		plt.pack();
		Rectangle s = plt.getGraphicsConfiguration().getBounds();
		Dimension f = plt.getSize();
		int w = Math.max(s.width - f.width, 0);
		int h = Math.max(s.height - f.height, 0);
		int x = (int) (0.5 * w) + s.x;
		int y = (int) (0.5 * h) + s.y;
		plt.setBounds(x, y, f.width, f.height);
		plt.setVisible(true);
	}*/

	public static void printMtx(double [][]mtx) {
		for (int i = 0; i < mtx.length; i++) {
			for (int j = 0; j < mtx[0].length; j++) {
				System.out.print(mtx[i][j] + ", ");
			}
			System.out.println();
		}
	}

	public static void printMtx(Object []mtx) {
		for (int i = 0; i < mtx.length; i++) {
			System.out.print(mtx[i].getClass() + "=>" + mtx[i] + ", ");
		}
		System.out.println();
	}

	public static void main(String [] args) throws IOException {
		tsne_csv(args);
	}

}
