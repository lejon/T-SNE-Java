package org.math.plot.plots;

import java.awt.Color;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import javax.imageio.ImageIO;

import org.math.plot.render.AbstractDrawer;

public class IconScatterPlot extends ScatterPlot {

	Random rand = new Random();
	String [] imagefiles; 
	Map<String,Color> colors = new HashMap<String,Color>();
	int colorCnt = 0;
	Image [] imgs;
	float alpha = 1.0f;
	double width = 0.4;
	double height = 0.4;
	
	public IconScatterPlot(String n, double[][] _XY, String [] imagefiles) {
		super(n, AbstractDrawer.DEFAULT_COLOR, AbstractDrawer.ROUND_DOT, AbstractDrawer.DEFAULT_DOT_RADIUS, _XY);
		
		this.imagefiles = imagefiles;
		imgs = new Image[imagefiles.length];
		for (int i = 0; i < imagefiles.length; i++) {
			BufferedImage img = null;
			try {
			    img = ImageIO.read(new File(imagefiles[i]));
			    imgs[i] = img;
			    if(img.getHeight()<=0 || img.getWidth()<=0) {
			    	throw new IllegalArgumentException("Image: " + imagefiles[i] 
			    			+ " is broken, width or height is <= 0!");
			    }
			} catch (IOException e) {
				System.err.println("File is: " + imagefiles[i]);
				e.printStackTrace();
			}
		}
	}
	
	public void plot(AbstractDrawer draw, Color col) {
        if (!visible) {
            return;
        }
        for (int i = 0; i < XY.length; i++) {
        	height = (double)draw.canvas.getHeight() / 500.0;
        	width  = (double)draw.canvas.getWidth() / 500.0;
        	double [] swc = {XY[i][0]-(width/2),XY[i][1]-(height/2)};
        	double [] sec = {XY[i][0]+(width/2),XY[i][1]-(height/2)};
        	double [] nwc = {XY[i][0]-(width/2),XY[i][1]+(height/2)};
            draw.drawImage(imgs[i], alpha, swc, sec, nwc);
        }
    }
	
//	public static void main(String[] args) {
//		Plot2DPanel p2 = new Plot2DPanel();
//		for (int i = 0; i < 1; i++) {
//			double[][] XYZ = new double[10][2];
//			for (int j = 0; j < XYZ.length; j++) {
//				XYZ[j][0] =/*1 + */Math.random();
//				XYZ[j][1] = /*100 * */Math.random();
//			}
//			p2.addScatterPlot("toto" + i, XYZ);
//		}
//		//public RasterImage(File _source, float _alpha,double[] _xyzSW, double[] _xyzSE,double[] _xyzNW) {
//		p2.addPlotable(new RasterImage(new File("imgs/img2.jpg"), 1.0f,new double[] {0.0,0.0},new double[] {0.05,0.0},new double[] {0.0,0.05}));
//
//		p2.setLegendOrientation(PlotPanel.SOUTH);
//		new FrameView(p2).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		
//		
//		
//		Plot3DPanel p = new Plot3DPanel();
//		for (int i = 0; i < 1; i++) {
//			double[][] XYZ = new double[10][3];
//			for (int j = 0; j < XYZ.length; j++) {
//				XYZ[j][0] = /*1 +*/ Math.random();
//				XYZ[j][1] = /*100 **/ Math.random();
//				XYZ[j][2] = /*0.0001 **/ Math.random();
//			}
//			p.addScatterPlot("toto" + i, XYZ);
//		}
//
//		p.addPlotable(new RasterImage(new File("imgs/img0.jpg"), 0.5f, new double[] {0.0,0.0,0.0},new double[] {1,0,0.0},new double[] {0.0,0,1}));
//		p.addPlotable(new RasterImage(new File("imgs/img1.jpg"), 0.5f, new double[] {0.0,0.0,0.0},new double[] {0,1,0.0},new double[] {0,0.0,1}));
//		p.addPlotable(new RasterImage(new File("imgs/img2.jpg"), 0.5f, new double[] {0.0,0.0,0.0},new double[] {1,0,0.0},new double[] {0,1.0,0}));
//		// TODO this following case is not totally supported...
//		//p.addPlotable(new PlotImage(new File("test.jpg"),0.5f,  new double[] {1,0,0},new double[] {1,1,0},new double[] {0,0,1}));
//
//		
//		p.setLegendOrientation(PlotPanel.SOUTH);
//		p.setPreferredSize(new Dimension(600,600));
//		new FrameView(p).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);		
//	}
}
