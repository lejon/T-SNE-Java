package org.math.plot.plots;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.math.plot.PlotPanel;
import org.math.plot.render.AbstractDrawer;
import org.math.plot.plots.ScatterPlot;

public class ColoredScatterPlot extends ScatterPlot {

	// New random with a seed so we get the same order of colors each time
	Random rand = new Random(4711);
	String [] classes; 
	Map<String,Color> colors = new HashMap<String,Color>();
	int colorCnt = 0;
	
	public ColoredScatterPlot(String n, double[][] _XY, String [] classes) {
		super(n, AbstractDrawer.DEFAULT_COLOR, AbstractDrawer.ROUND_DOT, AbstractDrawer.DEFAULT_DOT_RADIUS, _XY);
		this.classes = classes;
	}
	
	public void plot(AbstractDrawer draw, Color col) {
        if (!visible) {
            return;
        }
        for (int i = 0; i < XY.length; i++) {
        	Color c = colors.get(classes[i]);
        	if( c == null ) {
        		c = oldNewColor();
        		//c = newContrastColor();
        		//c = newColor();
        		colors.put(classes[i],c);
        		colorCnt++;
        		System.out.println(colorCnt + ": Class: " + classes[i] + " => "+ c);
        	}
        	draw.setColor(c);
            draw.drawDot(XY[i]);
        }
    }
	
	Color newColor() {
		// Java 'Color' class takes 3 floats, from 0 to 1.
		float r = rand.nextFloat();
		float g = rand.nextFloat();
		float b = rand.nextFloat();
		//Then to finally create the colour, pass the primary colours into the constructor:

		return new Color(r, g, b);
	}
	
	Color oldNewColor() {
		return PlotPanel.COLORLIST[colorCnt % PlotPanel.COLORLIST.length];
	}
	
	Color newContrastColor() {
		return kellysMaxContrastSet().get(colorCnt++ % kellysMaxContrastSet().size());
	}

	static List<Color> colorBrewer12s1()
	{
		ArrayList<Color> cols = new ArrayList<Color>();
		cols.add(UIntToColor(0xFFa6cee3));
		cols.add(UIntToColor(0xFF1f78b4));
		cols.add(UIntToColor(0xFFb2df8a));
		cols.add(UIntToColor(0xFF33a02c));
		cols.add(UIntToColor(0xFFfb9a99));
		cols.add(UIntToColor(0xFFe31a1c));
		cols.add(UIntToColor(0xFFfdbf6f));
		cols.add(UIntToColor(0xFFff7f00));
		cols.add(UIntToColor(0xFFcab2d6));
		cols.add(UIntToColor(0xFF6a3d9a));
		cols.add(UIntToColor(0xFFffff99));
		cols.add(UIntToColor(0xFFb15928));
		
		return cols;
	}
	
	static List<Color> colorBrewer12s2()
	{
		ArrayList<Color> cols = new ArrayList<Color>();
		cols.add(UIntToColor(0xFF8dd3c7));
		cols.add(UIntToColor(0xFFffffb3));
		cols.add(UIntToColor(0xFFbebada));
		cols.add(UIntToColor(0xFFfb8072));
		cols.add(UIntToColor(0xFF80b1d3));
		cols.add(UIntToColor(0xFFfdb462));
		cols.add(UIntToColor(0xFFb3de69));
		cols.add(UIntToColor(0xFFfccde5));
		cols.add(UIntToColor(0xFFd9d9d9));
		cols.add(UIntToColor(0xFFbc80bd));
		cols.add(UIntToColor(0xFFccebc5));
		cols.add(UIntToColor(0xFFffed6f));
		return cols;
	}
	
	static List<Color> kellysMaxContrastSet()
	{
		ArrayList<Color> cols = new ArrayList<Color>();
		cols.add(UIntToColor(0xFFFFB300)); //Vivid Yellow
		cols.add(UIntToColor(0xFF803E75)); //Strong Purple
		cols.add(UIntToColor(0xFFFF6800)); //Vivid Orange
		cols.add(UIntToColor(0xFFA6BDD7)); //Very Light Blue
		cols.add(UIntToColor(0xFFC10020)); //Vivid Red
		cols.add(UIntToColor(0xFFCEA262)); //Grayish Yellow
		cols.add(UIntToColor(0xFF817066)); //Medium Gray

		//The following will not be good for people with defective color vision
		cols.add(UIntToColor(0xFF007D34)); //Vivid Green
		cols.add(UIntToColor(0xFFF6768E)); //Strong Purplish Pink
		cols.add(UIntToColor(0xFF00538A)); //Strong Blue
		cols.add(UIntToColor(0xFFFF7A5C)); //Strong Yellowish Pink
		cols.add(UIntToColor(0xFF53377A)); //Strong Violet
		cols.add(UIntToColor(0xFFFF8E00)); //Vivid Orange Yellow
		cols.add(UIntToColor(0xFFB32851)); //Strong Purplish Red
		cols.add(UIntToColor(0xFFF4C800)); //Vivid Greenish Yellow
		cols.add(UIntToColor(0xFF7F180D)); //Strong Reddish Brown
		cols.add(UIntToColor(0xFF93AA00)); //Vivid Yellowish Green
		cols.add(UIntToColor(0xFF593315)); //Deep Yellowish Brown
		cols.add(UIntToColor(0xFFF13A13)); //Vivid Reddish Orange
		cols.add(UIntToColor(0xFF232C16)); //Dark Olive Green
		return cols;
	}

	
	static List<Color> boyntonOptimized()
	{
		ArrayList<Color> cols = new ArrayList<Color>();
		cols.add(new Color(0, 0, 255));      //Blue
		cols.add(new Color(255, 0, 0));      //Red
		cols.add(new Color(0, 255, 0));      //Green
		cols.add(new Color(255, 255, 0));    //Yellow
		cols.add(new Color(255, 0, 255));    //Magenta
		cols.add(new Color(255, 128, 128));  //Pink
		cols.add(new Color(128, 128, 128));  //Gray
		cols.add(new Color(128, 0, 0));      //Brown
		cols.add(new Color(255, 128, 0));    //Orange
		return cols;
	}

	public static Color UIntToColor(long color)
	{
	    int a = (int) (0x000000FF & (color >> 24));
	    int r = (int) (0x000000FF & (color >> 16));
	    int g = (int) (0x000000FF & (color >> 8));
	    int b = (int) (0x000000FF & (color >> 0));
	    return new Color(r, g, b, a);
	}
	
}
