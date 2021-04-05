package com.marginallyclever.makelangelo.nodes.fractals;

import com.marginallyclever.core.Translator;
import com.marginallyclever.core.node.NodeConnectorBoundedInt;
import com.marginallyclever.core.node.NodeConnectorInteger;
import com.marginallyclever.core.turtle.Turtle;
import com.marginallyclever.makelangelo.nodes.TurtleGenerator;

/**
 * @see <a href='https://en.wikipedia.org/wiki/Sierpi%C5%84ski_arrowhead_curve'>Wikipedia</a>
 * @author Dan Royer 
 * @since 2016-12-12
 *
 */
public class Generator_SierpinskiTriangle extends TurtleGenerator {
	// controls complexity of curve
	private NodeConnectorInteger inputOrder = new NodeConnectorBoundedInt(Translator.get("Fractal.inputOrder"),15,1,4);
	
	private double xMax, xMin, yMax, yMin;
	private double maxSize;
		
	public Generator_SierpinskiTriangle() {
		super();
		inputs.add(inputOrder);
		inputOrder.setDescription(Translator.get("Fractal.inputOrder.tooltip"));
	}
	
	@Override
	public String getName() {
		return Translator.get("Generator_SierpinskiTriangle.name");
	}
	
	@Override
	public boolean iterate() {
		Turtle turtle = new Turtle();

		double pw = inputWidth.getValue();
		double ph = inputHeight.getValue();
		yMin = -ph/2;
		yMax = ph/2;
		xMin = -pw/2;
		xMax = pw/2;

		turtle = new Turtle();
		
		double xx = xMax - xMin;
		double yy = yMax - yMin;
		maxSize = Math.tan(Math.toRadians(30))*(xx < yy ? xx : yy)*2;
		double jj = Math.asin(Math.toRadians(30))*(xx < yy ? xx : yy);

		// move to starting position
		if(xMax>yMax) {
			turtle.moveTo(-jj,yMin);
		} else {
			turtle.moveTo(xMax,-jj);
			turtle.turn(90);
		}
		turtle.penDown();
		// do the curve
		int order = inputOrder.getValue();
		
		if( (order&1) == 0 ) {
			drawCurve(turtle,order, maxSize,-60);
		} else {
			turtle.turn(60);
			drawCurve(turtle,order, maxSize,-60);
		}

		outputTurtle.setValue(turtle);
		
	    return false;
	}


	private void drawCurve(Turtle turtle,int n, double distance,double angle) {
		if (n == 0) {
			turtle.forward(distance);
			return;
		}
		
		drawCurve(turtle,n-1,distance/2.0f,-angle);
		turtle.turn(angle);
		drawCurve(turtle,n-1,distance/2.0f,angle);
		turtle.turn(angle);
		drawCurve(turtle,n-1,distance/2.0f,-angle);
	}
}
