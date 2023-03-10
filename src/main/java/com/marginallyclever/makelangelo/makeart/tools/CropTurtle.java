package com.marginallyclever.makelangelo.makeart.tools;

import com.marginallyclever.convenience.Clipper2D;
import com.marginallyclever.convenience.MathHelper;
import com.marginallyclever.convenience.Point2D;
import com.marginallyclever.makelangelo.turtle.MovementType;
import com.marginallyclever.makelangelo.turtle.Turtle;
import com.marginallyclever.makelangelo.turtle.TurtleMove;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;

@Deprecated
public class CropTurtle {
	private static final Logger logger = LoggerFactory.getLogger(CropTurtle.class);
	
	public static void run(Turtle turtle,Rectangle2D.Double rectangle) {
		logger.debug("crop start @ {}", turtle.history.size());

		List<TurtleMove> newHistory = new ArrayList<>();
		// limits we will need for rectangle
		Point2D rMax = new Point2D(rectangle.getMaxX(),rectangle.getMaxY());
		Point2D rMin = new Point2D(rectangle.getMinX(),rectangle.getMinY());
		// working space for clipping
		Point2D P0 = new Point2D(); 
		Point2D P1 = new Point2D(); 
		
		TurtleMove prev=null;
		
		for (TurtleMove m : turtle.history ) {
			switch (m.type) {
				case DRAW_LINE, TRAVEL -> {
					if (prev != null) {
						P0.set(prev.x, prev.y);
						P1.set(m.x, m.y);
						boolean result = Clipper2D.clipLineToRectangle(P0, P1, rMax, rMin);
						// !result means full crop, do nothing.
						if (result) {
							// partial crop.  Which end(s)?
							boolean startCropped = MathHelper.lengthSquared(P0.x - prev.x, P0.y - prev.y) > 1e-8;
							boolean endCropped = MathHelper.lengthSquared(P1.x - m.x, P1.y - m.y) > 1e-8;

							if (startCropped && endCropped) {
								// crosses rectangle, both ends out.
								newHistory.add(new TurtleMove(P0.x, P0.y, MovementType.TRAVEL));
								newHistory.add(m);
								TurtleMove m2 = new TurtleMove(P1.x, P1.y, m.type);
								newHistory.add(m2);
							} else if (!startCropped && !endCropped) {
								newHistory.add(m);
							} else if (endCropped) {
								// end cropped, leaving the rectangle
								TurtleMove m2 = new TurtleMove(P1.x, P1.y, m.type);
								newHistory.add(m2);
							} else {
								// start cropped, coming back into rectangle
								newHistory.add(new TurtleMove(P0.x, P0.y, MovementType.TRAVEL));
								newHistory.add(m);
							}
						}
					}
					prev = m;
				}
				default -> newHistory.add(m);
			}
		}

		turtle.history.clear();
		turtle.history.addAll(newHistory);
		
		// There may be some dumb travel moves left. (several travels in a row.)
	
		logger.debug("crop end @ {}", turtle.history.size());

	}
}
