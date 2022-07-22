package com.marginallyclever.makelangelo.makeart.io.vector;

import com.marginallyclever.convenience.Bezier;
import com.marginallyclever.convenience.ColorRGB;
import com.marginallyclever.convenience.Point2D;
import com.marginallyclever.makelangelo.turtle.Turtle;
import io.sf.carte.echosvg.anim.dom.*;
import io.sf.carte.echosvg.bridge.*;
import io.sf.carte.echosvg.dom.svg.SVGItem;
import io.sf.carte.echosvg.dom.util.SAXDocumentFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.w3c.dom.svg.*;

import javax.swing.filechooser.FileNameExtensionFilter;
import javax.vecmath.Matrix3d;
import javax.vecmath.Vector3d;
import java.awt.geom.Rectangle2D;
import java.io.InputStream;
import java.util.List;

/**
 * @author Dan Royer
 * See <a href="https://www.w3.org/TR/SVG/paths.html">w3.org</a>
 */
public class LoadSVG implements TurtleLoader {
	private static final Logger logger = LoggerFactory.getLogger(LoadSVG.class);

	private static final String LABEL_STROKE="stroke:";

	private static final FileNameExtensionFilter filter = new FileNameExtensionFilter("Scaleable Vector Graphics 1.1", "svg");
	private Turtle myTurtle;

	private boolean isNewPath;  // for cubic paths
	private final Vector3d pathFirstPoint = new Vector3d();
	private final Vector3d pathPoint = new Vector3d();

	@Override
	public FileNameExtensionFilter getFileNameFilter() {
		return filter;
	}

	@Override
	public boolean canLoad(String filename) {
		String ext = filename.substring(filename.lastIndexOf('.'));
		return ext.equalsIgnoreCase(".svg");
	}

	@Override
	public Turtle load(InputStream in) throws Exception {
		if (in == null) {
			throw new NullPointerException("Input stream is null");
		}

		logger.debug("Loading...");

		Document document = newDocumentFromInputStream(in);
		initSVGDOM(document);

		myTurtle = new Turtle();
		myTurtle.setColor(new ColorRGB(0,0,0));
		parseAll(document);

		Rectangle2D.Double r = myTurtle.getBounds();
		myTurtle.translate(-r.width/2,-r.height/2);
		myTurtle.scale(1, -1);

		return myTurtle;
	}

	private void parseAll(Document document) throws Exception {
		SVGOMSVGElement documentElement = (SVGOMSVGElement)document.getDocumentElement();

		logger.debug("...parse path");			parsePathElements(    documentElement.getElementsByTagName( "path"     ));
		logger.debug("...parse polylines");		parsePolylineElements(documentElement.getElementsByTagName( "polyline" ));
		logger.debug("...parse polygons");		parsePolylineElements(documentElement.getElementsByTagName( "polygon"  ));
		logger.debug("...parse lines");			parseLineElements(    documentElement.getElementsByTagName( "line"     ));
		logger.debug("...parse rects");			parseRectElements(    documentElement.getElementsByTagName( "rect"     ));
		logger.debug("...parse circles");		parseCircleElements(  documentElement.getElementsByTagName( "circle"   ));
		logger.debug("...parse ellipses");		parseEllipseElements( documentElement.getElementsByTagName( "ellipse"  ));
	}

	/**
	 * Parse through all the SVG polyline elements and raster them to gcode.
	 * @param pathNodes the source of the elements
	 */
	private void parsePolylineElements(NodeList pathNodes) {
	    int pathNodeCount = pathNodes.getLength();
		logger.debug("{} elements", pathNodeCount);
	    for( int iPathNode = 0; iPathNode < pathNodeCount; iPathNode++ ) {
	    	SVGPointShapeElement element = (SVGPointShapeElement)pathNodes.item( iPathNode );
			if(isElementStrokeNone(element)) 
				continue;

			Matrix3d m = getMatrixFromElement(element);

			SVGPointList pointList = element.getAnimatedPoints();
			int numPoints = pointList.getNumberOfItems();
			//logger.debug("New Node has "+pathObjects+" elements.");

			SVGPoint item = pointList.getItem(0);
			Vector3d v2 = transform(item.getX(),item.getY(),m);
			myTurtle.jumpTo(v2.x,v2.y);

			for( int i=1; i<numPoints; ++i ) {
				item = pointList.getItem(i);
				v2 = transform(item.getX(),item.getY(),m);
				myTurtle.moveTo(v2.x,v2.y);
			}
		}
	}

	private void parseLineElements(NodeList node) {
		Vector3d v2;
	    int pathNodeCount = node.getLength();
		logger.debug("{} elements", pathNodeCount);
	    for( int iPathNode = 0; iPathNode < pathNodeCount; iPathNode++ ) {
			Element element = (Element)node.item( iPathNode );
			if(isElementStrokeNone(element))
				continue;

			Matrix3d m = getMatrixFromElement(element);

			double x1=0,y1=0;
			double x2=0,y2=0;

			if(element.hasAttribute("x1")) x1 = Double.parseDouble(element.getAttribute("x1"));
			if(element.hasAttribute("y1")) y1 = Double.parseDouble(element.getAttribute("y1"));
			if(element.hasAttribute("x2")) x2 = Double.parseDouble(element.getAttribute("x2"));
			if(element.hasAttribute("y2")) y2 = Double.parseDouble(element.getAttribute("y2"));
			v2 = transform(x1,y1,m);
			myTurtle.jumpTo(v2.x,v2.y);
			v2 = transform(x2,y2,m);
			myTurtle.moveTo(v2.x,v2.y);
		}
	}

	private boolean isElementStrokeNone(Element element) {
		if(element.hasAttribute("style")) {
			String style = element.getAttribute("style").toLowerCase().replace("\s","");
			if(style.contains(LABEL_STROKE)) {
				int k = style.indexOf(LABEL_STROKE);
				String strokeStyleName = style.substring(k+LABEL_STROKE.length());
				// it is!  bail.
				return strokeStyleName.contentEquals("none") || strokeStyleName.contentEquals("white");
			} else {
				// default SVG stroke is "none", which isn't even transparent - it's nothing!
				return false;
			}
		}
		return false;
	}

	/**
	 * Draw rectangles that may have rounded corners.
	 * given corners
	 *    x0 x1 x2 x3
	 * y0    a  b
	 * y1 c  i  j  d
	 * y2 e  m  k  f
	 * y3    g  h
	 * draw a-b-d-f-h-g-e-c-a.
	 *
	 * See <a href="https://developer.mozilla.org/en-US/docs/Web/SVG/Element/rect">SVG specification</a>
	 */
	private void parseRectElements(NodeList node) {
	    int pathNodeCount = node.getLength();
		logger.debug("{} elements", pathNodeCount);
	    for( int iPathNode = 0; iPathNode < pathNodeCount; iPathNode++ ) {
			Element element = (Element)node.item( iPathNode );
			if(isElementStrokeNone(element))
				continue;

			Matrix3d m = getMatrixFromElement(element);

			double x=0,y=0;
			double rx=0,ry=0;

			if(element.hasAttribute("x")) x = Double.parseDouble(element.getAttribute("x"));
			if(element.hasAttribute("y")) y = Double.parseDouble(element.getAttribute("y"));
			if(element.hasAttribute("rx")) {
				rx = Double.parseDouble(element.getAttribute("rx"));
				if(element.hasAttribute("ry")) {
					ry = Double.parseDouble(element.getAttribute("ry"));
				} else {
					// ry defaults to rx if specified
					ry = rx;
				}
			} else if(element.hasAttribute("ry")) {
				// rx defaults to ry if specified
				rx = ry = Double.parseDouble(element.getAttribute("ry"));

			}
			double w = Double.parseDouble(element.getAttribute("width"));
			double h = Double.parseDouble(element.getAttribute("height"));

			//double x0=x;
			double x1=x+rx;
			double x2=x+w-rx;
			//double x3=x+w;
			double y0=y;
			double y1=y+ry;
			double y2=y+h-ry;
			//double y3=y+h;

			Vector3d v2 = transform(x1,y0,m);
			myTurtle.jumpTo(v2.x,v2.y);
			arcTurtle(myTurtle, x2,y1, rx,ry, Math.PI * -0.5,Math.PI *  0.0,m);
			arcTurtle(myTurtle, x2,y2, rx,ry, Math.PI *  0.0,Math.PI *  0.5,m);
			arcTurtle(myTurtle, x1,y2, rx,ry, Math.PI * -1.5,Math.PI * -1.0,m);
			arcTurtle(myTurtle, x1,y1, rx,ry, Math.PI * -1.0,Math.PI * -0.5,m);
		}
	}

	/**
	 *
	 * @param cx center position
	 * @param cy center position
	 * @param rx radius on X
	 * @param ry radius on Y
	 * @param p0 radian start angle.
	 * @param p1 radian end angle.
	 */
	private void arcTurtle(Turtle turtle,double cx,double cy,double rx,double ry,double p0,double p1,Matrix3d m) {
		Vector3d v2;
		double steps=1;
		if(rx>0 && ry>0) {
			double r = Math.max(rx, ry);
			double circ = Math.PI*r*2.0;  // radius to circumference
			steps = Math.ceil(circ/4.0);  // 1/4 circumference
			steps = Math.max(steps,1);
		}
		steps = steps/4;
		for(double p = 0;p<=steps;++p) {
			double pFraction = ((p1-p0)*(p/steps) + p0);
			double c = Math.cos(pFraction) * rx;
			double s = Math.sin(pFraction) * ry;
			v2 = transform(cx+c,cy+s,m);
			turtle.moveTo(v2.x,v2.y);
		}
	}

	private void parseCircleElements(NodeList node) {
		Vector3d v2;

	    int pathNodeCount = node.getLength();
		logger.debug("{} elements", pathNodeCount);
	    for( int iPathNode = 0; iPathNode < pathNodeCount; iPathNode++ ) {
			Element element = (Element)node.item( iPathNode );
			if(isElementStrokeNone(element))
				continue;

			Matrix3d m = getMatrixFromElement(element);

			double cx=0,cy=0,r=0;
			if(element.hasAttribute("cx")) cx = Double.parseDouble(element.getAttribute("cx"));
			if(element.hasAttribute("cy")) cy = Double.parseDouble(element.getAttribute("cy"));
			if(element.hasAttribute("r" )) r  = Double.parseDouble(element.getAttribute("r"));
			v2 = transform(cx+r,cy,m);
			myTurtle.jumpTo(v2.x,v2.y);

			double circ = Math.PI * 2.0 * r;
			circ = Math.ceil(Math.min(Math.max(3,circ),360));

			logger.debug("circ={}", circ);
			printEllipse(m, cx, cy, r, r, circ);
		}
	}

	private void parseEllipseElements(NodeList node) {
		Vector3d v2;
	    int pathNodeCount = node.getLength();
		logger.debug("{} elements", pathNodeCount);
	    for( int iPathNode = 0; iPathNode < pathNodeCount; iPathNode++ ) {
			Element element = (Element)node.item( iPathNode );
			if(isElementStrokeNone(element))
				continue;

			Matrix3d m = getMatrixFromElement(element);

			double cx=0,cy=0,rx=0,ry=0;
			if(element.hasAttribute("cx")) cx = Double.parseDouble(element.getAttribute("cx"));
			if(element.hasAttribute("cy")) cy = Double.parseDouble(element.getAttribute("cy"));
			if(element.hasAttribute("rx")) rx = Double.parseDouble(element.getAttribute("rx"));
			if(element.hasAttribute("ry")) ry = Double.parseDouble(element.getAttribute("ry"));
			v2 = transform(cx+rx,cy,m);
			myTurtle.jumpTo(v2.x,v2.y);

			double perimeterOfAnEllipseApprox = Math.PI * 2.0 * Math.sqrt((ry*ry + rx*rx)/2.0);
			double steps = Math.max(3,perimeterOfAnEllipseApprox);
			steps = Math.min(60,steps);
			printEllipse(m, cx, cy, rx, ry, steps);
		}
	}

	private void printEllipse(Matrix3d m, double cx, double cy, double rx, double ry, double steps) {
		Vector3d v2;
		for(double i = 1; i<steps; ++i) {
			double v = (Math.PI*2.0) * (i/steps);
			double s=ry*Math.sin(v);
			double c=rx*Math.cos(v);
			v2 = transform(cx+c,cy+s,m);
			myTurtle.moveTo(v2.x,v2.y);
		}
		v2 = transform(cx+rx,cy,m);
		myTurtle.moveTo(v2.x,v2.y);
	}

	/**
	 * Parse through all the SVG path elements and raster them to {@link Turtle}.
	 * @param paths the source of the elements
	 */
	private void parsePathElements(NodeList paths) throws Exception {
		int pathCount = paths.getLength();
		logger.debug("{} elements", pathCount);
		for( int iPath = 0; iPath < pathCount; iPath++ ) {
			if(paths.item( iPath ) instanceof SVGOMPolylineElement) {
				logger.debug("Node is a polyline.");
				parsePolylineElements(paths);
				continue;
			}
			SVGOMPathElement element = ((SVGOMPathElement)paths.item( iPath ));
			if(isElementStrokeNone(element))
				continue;

			Matrix3d m = getMatrixFromElement(element);

			SVGPathSegList pathList = element.getPathSegList();
			int itemCount = pathList.getNumberOfItems();
			logger.debug("Node has {} elements.", itemCount);
			isNewPath=true;

			for(int i=0; i<itemCount; i++) {
				SVGPathSeg item = pathList.getItem(i);
				logger.debug(((SVGItem)item).getValueAsString());
				switch( item.getPathSegType() ) {
					case SVGPathSeg.PATHSEG_MOVETO_ABS 			-> doMoveToAbs(item,m);  	// M
					case SVGPathSeg.PATHSEG_MOVETO_REL 			-> doMoveToRel(item,m);     // m
					case SVGPathSeg.PATHSEG_LINETO_ABS 			-> doLineToAbs(item,m);  	// L
					case SVGPathSeg.PATHSEG_LINETO_REL 			-> doLineToRel(item,m);  	// l
					case SVGPathSeg.PATHSEG_CURVETO_CUBIC_ABS 	-> doCubicCurveAbs(item,m);	// C
					case SVGPathSeg.PATHSEG_CURVETO_CUBIC_REL 	-> doCubicCurveRel(item,m);	// c
					case SVGPathSeg.PATHSEG_CLOSEPATH 			-> doClosePath(); 			// Z z
					case SVGPathSeg.PATHSEG_ARC_ABS				-> doArcAbs(item,m);  // A
					//case SVGPathSeg.PATHSEG_ARC_REL				-> doArcRel(item,m);  // a
					case SVGPathSeg.PATHSEG_LINETO_VERTICAL_ABS -> doLineToVerticalAbs(item,m);  // V
					case SVGPathSeg.PATHSEG_LINETO_VERTICAL_REL -> doLineToVerticalRel(item,m);  // v
					case SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_ABS -> doLineToHorizontalAbs(item,m);  // H
					case SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_REL -> doLineToHorizontalRel(item,m);  // h

					default -> throw new Exception("Found unknown SVGPathSeg type "+item.getPathSegType());
				}
			}
		}
	}

	private double [] polarToCartesian(double cx,double cy,double radius,double degrees) {
		var radians = Math.toRadians(degrees);
		var x = cx + radius * Math.cos(radians);
		var y = cy + radius * Math.sin(radians);
		return new double[]{x,y};
	}

	/**
	 * Calculate the angle between vector a vector b
	 *
	 * @param u
	 * @param v
	 * @return the angle in radian
	 */
	public static double getAngle(Vector3d u, Vector3d v) {
		double cos = u.dot(v) / (u.length() * v.length());
		double result = Math.acos(cos);
		double sign = Math.signum(u.x*v.y - u.y*v.x);
		return sign * Math.abs(result);
	}

	private void doArcAbs(SVGPathSeg item, Matrix3d m) {
		var path = (SVGPathSegArcAbs) item;
		boolean sweep = path.getSweepFlag();
		boolean large = path.getLargeArcFlag();
		var r1 = path.getR1();
		var r2 = path.getR2();
		var angleDegrees = path.getAngle();
		var end = transform(path.getX(), path.getY(), m);

		// with start (aka pathPoint), end, rx, ry find the two centers.
		getArcCenter(pathPoint, end, r1, r2, angleDegrees, large, sweep);
	}

	private void getArcCenter(Vector3d start, Vector3d end, double rx, double ry, double angleDegrees, boolean fA, boolean fS) {
		// translate to mid point of start-end
		var prime = new Vector3d(
				(start.x-end.x)/2.0,
				(start.y-end.y)/2.0,
				0);
		// rotate
		Matrix3d m = new Matrix3d();
		m.rotZ(Math.toRadians(-angleDegrees));
		m.transform(prime);
		var x1p = prime.x;
		var y1p = prime.y;

		// Ensure radii are large enough
		// Based on http://www.w3.org/TR/SVG/implnote.html#ArcOutOfRangeParameters
		// Step (a): Ensure radii are non-zero
		// Step (b): Ensure radii are positive
		rx = Math.abs(rx);
		ry = Math.abs(ry);
		// Step (c): Ensure radii are large enough
		var lambda = ( (x1p * x1p) / (rx * rx) ) + ( (y1p * y1p) / (ry * ry) );
		if(lambda > 1) {
			rx = Math.sqrt(lambda) * rx;
			ry = Math.sqrt(lambda) * ry;
		}


		// Step 2: Compute (cx′, cy′)
		var sign = (fA == fS)? -1 : 1;

		// Bit of a hack, as presumably rounding errors were making his negative inside the square root!
		double v = ((rx * rx * ry * ry) - (rx * rx * y1p * y1p) - (ry * ry * x1p * x1p)) / ((rx * rx * y1p * y1p) + (ry * ry * x1p * x1p));
		double co = (v < 1e-7) ? 0 : sign * Math.sqrt(v);
		Vector3d Cprime = new Vector3d(rx*y1p/ry,-ry*x1p/rx,0);
		Cprime.scale(v * co);

		// Step 3: Compute (cx, cy) from (cx′, cy′)
		prime.set(Cprime);
		m.invert();
		m.transform(prime);

		var middle2 = new Vector3d(
				(start.x+end.x)/2.0,
				(start.y+end.y)/2.0,
				0);
		middle2.add(prime);

		// Step 4: compute start angle and sweep angle
		Vector3d s1 = new Vector3d(
				(x1p-Cprime.x)/rx,
				(y1p-Cprime.y)/ry,
				0);
		double startAngleRadians = getAngle(new Vector3d(1,0,0),s1);

		Vector3d s2 = new Vector3d(
				(-x1p-Cprime.x)/rx,
				(-y1p-Cprime.y)/ry,
				0);
		double sweepRadians = getAngle(s1,s2) % (Math.PI*2.0);
		if(!fS) {
			if(sweepRadians > 0) sweepRadians -= (Math.PI*2.0);
		} else {
			if(sweepRadians < 0) sweepRadians += (Math.PI*2.0);
		}


		int steps = (int)Math.abs(sweepRadians)*10;
		double angleRadians = Math.toRadians(angleDegrees);
		Vector3d rx2 = new Vector3d( Math.cos(angleRadians), Math.sin(angleRadians),0 );
		Vector3d ry2 = new Vector3d( -rx2.y, rx2.x, 0 );
		rx2.scale(rx);
		ry2.scale(ry);

		Vector3d rx3 = new Vector3d();
		Vector3d ry3 = new Vector3d();
		Vector3d c2 = new Vector3d();

		logger.debug("angleDegrees={}",angleDegrees);
		logger.debug("steps={}",steps);

		for(int theta = 0; theta < steps; ++theta ) {
			double vv = sweepRadians * ((double)theta / (double)steps);
			double r = startAngleRadians + (fS ? vv : -vv);
			rx3.set(rx2);
			rx3.scale(Math.cos(r));

			ry3.set(ry2);
			ry3.scale(Math.sin(r));

			c2.set(middle2);
			c2.add(rx3);
			c2.add(ry3);

			myTurtle.moveTo(c2.x,c2.y);
		}

		pathPoint.set(end);
		myTurtle.moveTo(pathPoint.x,pathPoint.y);
	}

	private void doLineToVerticalAbs(SVGPathSeg item, Matrix3d m) {
		var path = (SVGPathSegLinetoVerticalAbs)item;
		Vector3d p = transform(0,path.getY(),m);
		pathPoint.y = p.y;
		myTurtle.moveTo(pathPoint.x,pathPoint.y);
	}

	private void doLineToVerticalRel(SVGPathSeg item, Matrix3d m) {
		var path = (SVGPathSegLinetoVerticalAbs)item;
		Vector3d p = transform(0,path.getY(),m);
		pathPoint.y += p.y;
		myTurtle.moveTo(pathPoint.x,pathPoint.y);
	}

	private void doLineToHorizontalAbs(SVGPathSeg item, Matrix3d m) {
		var path = (SVGPathSegLinetoHorizontalAbs)item;
		Vector3d p = transform(path.getX(),0,m);
		pathPoint.x = p.x;
		myTurtle.moveTo(pathPoint.x,pathPoint.y);
	}

	private void doLineToHorizontalRel(SVGPathSeg item, Matrix3d m) {
		var path = (SVGPathSegLinetoHorizontalAbs)item;
		Vector3d p = transform(path.getX(),0,m);
		pathPoint.x += p.x;
		myTurtle.moveTo(pathPoint.x,pathPoint.y);
	}

	private void doMoveToAbs(SVGPathSeg item, Matrix3d m) {
		var path = (SVGPathSegMovetoAbs)item;
		Vector3d p = transform(path.getX(),path.getY(),m);
		pathPoint.set(p);
		myTurtle.jumpTo(pathPoint.x,pathPoint.y);

		rememberIfFirstMove(pathPoint);
	}

	private void doMoveToRel(SVGPathSeg item, Matrix3d m) {
		var path = (SVGPathSegMovetoRel)item;
		Vector3d p = transform(path.getX(),path.getY(),m);
		pathPoint.add(p);
		myTurtle.jumpTo(pathPoint.x,pathPoint.y);

		rememberIfFirstMove(pathPoint);
	}

	private void rememberIfFirstMove(Vector3d p) {
		if(isNewPath) pathFirstPoint.set(p);
		isNewPath=false;
	}

	private void doLineToRel(SVGPathSeg item, Matrix3d m) {
		var path = (SVGPathSegLinetoRel)item;
		Vector3d p = transform(path.getX(),path.getY(),m);
		pathPoint.add(p);
		myTurtle.moveTo(pathPoint.x,pathPoint.y);
	}

	private void doLineToAbs(SVGPathSeg item, Matrix3d m) {
		var path = (SVGPathSegLinetoAbs)item;
		Vector3d p = transform(path.getX(),path.getY(),m);
		pathPoint.set(p);
		myTurtle.moveTo(pathPoint.x,pathPoint.y);
	}

	private void doCubicCurveAbs(SVGPathSeg item, Matrix3d m) {
		var path = (SVGPathSegCurvetoCubicAbs)item;
		// x0,y0 is the first point
		Vector3d p0 = pathPoint;
		// x1,y1 is the first control point
		Vector3d p1 = transform(path.getX1(),path.getY1(),m);
		// x2,y2 is the second control point
		Vector3d p2 = transform(path.getX2(),path.getY2(),m);
		// x3,y3 is the end point
		Vector3d p3 = transform(path.getX(),path.getY(),m);

		Bezier b = new Bezier(
				p0.x,p0.y,
				p1.x,p1.y,
				p2.x,p2.y,
				p3.x,p3.y);
		List<Point2D> points = b.generateCurvePoints(0.1);
		for(Point2D p : points) myTurtle.moveTo(p.x,p.y);
		pathPoint.set(p3);
	}

	private void doCubicCurveRel(SVGPathSeg item, Matrix3d m) {
		var path = (SVGPathSegCurvetoCubicRel)item;
		// x0,y0 is the first point
		Vector3d p0 = pathPoint;
		// x1,y1 is the first control point
		Vector3d p1 = transform(path.getX1(),path.getY1(),m);
		// x2,y2 is the second control point
		Vector3d p2 = transform(path.getX2(),path.getY2(),m);
		// x3,y3 is the end point
		Vector3d p3 = transform(path.getX(),path.getY(),m);

		p1.add(p0);
		p2.add(p0);
		p3.add(p0);

		Bezier b = new Bezier(
				p0.x,p0.y,
				p1.x,p1.y,
				p2.x,p2.y,
				p3.x,p3.y);
		List<Point2D> points = b.generateCurvePoints(0.1);
		for(Point2D p : points) myTurtle.moveTo(p.x,p.y);
		pathPoint.set(p3);
	}

	private void doClosePath() {
		myTurtle.moveTo(pathFirstPoint.x,pathFirstPoint.y);
		isNewPath=true;
	}

	private Vector3d transform(double x, double y, Matrix3d m) {
		Vector3d p = new Vector3d(x,y,1);
		m.transform(p);
		return p;
	}

	private Matrix3d getMatrixFromElement(Element element) {
		if(!(element instanceof SVGGraphicsElement)) {
			Matrix3d m = new Matrix3d();
			m.setIdentity();
			return m;
		}

		Matrix3d m = new Matrix3d();

		try {
			SVGGraphicsElement svgge = (SVGGraphicsElement)element;

			SVGMatrix svgMatrix = svgge.getCTM();
			// [ a c e ]
			// [ b d f ]
			// [ 0 0 1 ]
			m.m00 = svgMatrix.getA();	m.m10 = svgMatrix.getB();	m.m20 = 0;
			m.m01 = svgMatrix.getC();	m.m11 = svgMatrix.getD();	m.m21 = 0;
			m.m02 = svgMatrix.getE();	m.m12 = svgMatrix.getF();	m.m22 = 1;
		}
		catch(Exception e) {
			m.setIdentity();
		}
		return m;
	}

	/**
	 * Enhance the SVG DOM for the given document to provide CSS- and
	 * SVG-specific DOM interfaces.
	 * @param document The document to enhance.
	 * @link <a href="https://cwiki.apache.org/confluence/display/XMLGRAPHICSBATIK/BootSvgAndCssDom">apache.org</a>
	 */
	private void initSVGDOM(Document document) {
		UserAgent userAgent = new UserAgentAdapter();
		DocumentLoader loader = new DocumentLoader(userAgent);
		BridgeContext bridgeContext = new BridgeContext(userAgent, loader);
		bridgeContext.setDynamicState(BridgeContext.STATIC);

		// Enable CSS- and SVG-specific enhancements.
		(new GVTBuilder()).build(bridgeContext, document);
	}

	private static SVGDocument newDocumentFromInputStream(InputStream in) throws Exception {
		SAXDocumentFactory factory = new SAXSVGDocumentFactory();
		return (SVGDocument)factory.createDocument(null,in);
	}
}
