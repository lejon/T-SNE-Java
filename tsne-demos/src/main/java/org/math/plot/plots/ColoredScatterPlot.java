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

	Random rand = new Random();
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
        		c = newColor();
        		colors.put(classes[i],c);
        		colorCnt++;
        	}
        	draw.setColor(c);
            //draw.drawDot(XY[i]);
            draw.drawText(classes[i], XY[i]);
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
	    int a = (byte)(color >> 24);
	    int r = (byte)(color >> 16);
	    int g = (byte)(color >> 8);
	    int b = (byte)(color >> 0);
	    return new Color(r, g, b,a);
	}
	
}
