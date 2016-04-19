package org.math.plot.plots;

import java.awt.Color;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.math.plot.render.AbstractDrawer;

public class ColoredScatterPlot extends ScatterPlot {

	// New random with a seed so we get the same order of colors each time
	Random rand = new Random(4711);
	String [] classes; 
	Map<String,Color> colors = new HashMap<String,Color>();
	int colorCnt = 0;
	Color [] palette;
	
	public ColoredScatterPlot(String n, double[][] _XY, String [] classes) {
		super(n, AbstractDrawer.DEFAULT_COLOR, AbstractDrawer.ROUND_DOT, AbstractDrawer.DEFAULT_DOT_RADIUS, _XY);
		this.classes = classes;
		Set<String> unique = new HashSet<>(Arrays.asList(classes));
		boolean colorBlindSave = false;
        ColorBrewer[] sequentialPalettes = ColorBrewer.getDivergingColorPalettes(colorBlindSave);  

        ColorBrewer myBrewer = sequentialPalettes[0];

        System.out.println( "Name of this color brewer: " + myBrewer);

        Color[] myGradient = myBrewer.getColorPalette(unique.size());

        palette = myGradient;
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
        		System.out.println(colorCnt + ": Class: " + classes[i] + " => "+ c);
        	}
        	draw.setColor(c);
            draw.drawDot(XY[i]);
        }
    }
	
	Color newColor() {
		return palette[colorCnt++];
	}
	
}
