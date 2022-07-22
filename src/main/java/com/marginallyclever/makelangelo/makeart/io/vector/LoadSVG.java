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

	private boolean isNewSubPath;  // for cubic paths
	private final Vector3d pathFirstPoint = new Vector3d();
	private final Vector3d pathPoint = new Vector3d();

	private final Vector3d smoothControlPointA = new Vector3d();
	private final Vector3d smoothControlPointB = new Vector3d();
	private boolean previousPathSegWasBezier;

	private Matrix3d myMatrix;

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
		pathPoint.set(myTurtle.getX(),myTurtle.getY(),0);
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

			myMatrix = getMatrixFromElement(element);

			SVGPointList pointList = element.getAnimatedPoints();
			int numPoints = pointList.getNumberOfItems();
			//logger.debug("New Node has "+pathObjects+" elements.");

			SVGPoint item = pointList.getItem(0);
			jumpTo(item.getX(),item.getY());

			for( int i=1; i<numPoints; ++i ) {
				item = pointList.getItem(i);
				moveTo(item.getX(),item.getY());
			}
		}
	}

	private void parseLineElements(NodeList node) {
	    int pathNodeCount = node.getLength();
		logger.debug("{} elements", pathNodeCount);
	    for( int iPathNode = 0; iPathNode < pathNodeCount; iPathNode++ ) {
			Element element = (Element)node.item( iPathNode );
			if(isElementStrokeNone(element))
				continue;

			myMatrix = getMatrixFromElement(element);

			double x1=0,y1=0;
			double x2=0,y2=0;

			if(element.hasAttribute("x1")) x1 = Double.parseDouble(element.getAttribute("x1"));
			if(element.hasAttribute("y1")) y1 = Double.parseDouble(element.getAttribute("y1"));
			if(element.hasAttribute("x2")) x2 = Double.parseDouble(element.getAttribute("x2"));
			if(element.hasAttribute("y2")) y2 = Double.parseDouble(element.getAttribute("y2"));
			jumpTo(x1,y1);
			moveTo(x2,y2);
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

			myMatrix = getMatrixFromElement(element);

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

			jumpTo(x1,y0);
			arcTurtle(x2,y1, rx,ry, Math.PI * -0.5,Math.PI *  0.0);
			arcTurtle(x2,y2, rx,ry, Math.PI *  0.0,Math.PI *  0.5);
			arcTurtle(x1,y2, rx,ry, Math.PI * -1.5,Math.PI * -1.0);
			arcTurtle(x1,y1, rx,ry, Math.PI * -1.0,Math.PI * -0.5);
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
	private void arcTurtle(double cx,double cy,double rx,double ry,double p0,double p1) {
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
			moveTo(cx+c,cy+s);
		}
	}

	private void parseCircleElements(NodeList node) {
	    int pathNodeCount = node.getLength();
		logger.debug("{} elements", pathNodeCount);
	    for( int iPathNode = 0; iPathNode < pathNodeCount; iPathNode++ ) {
			Element element = (Element)node.item( iPathNode );
			if(isElementStrokeNone(element))
				continue;

			myMatrix = getMatrixFromElement(element);

			double cx=0,cy=0,r=0;
			if(element.hasAttribute("cx")) cx = Double.parseDouble(element.getAttribute("cx"));
			if(element.hasAttribute("cy")) cy = Double.parseDouble(element.getAttribute("cy"));
			if(element.hasAttribute("r" )) r  = Double.parseDouble(element.getAttribute("r"));
			jumpTo(cx+r,cy);

			double circ = Math.PI * 2.0 * r;
			circ = Math.ceil(Math.min(Math.max(3,circ),360));

			logger.debug("circ={}", circ);
			printEllipse(cx, cy, r, r, circ);
		}
	}

	private void parseEllipseElements(NodeList node) {
	    int pathNodeCount = node.getLength();
		logger.debug("{} elements", pathNodeCount);
	    for( int iPathNode = 0; iPathNode < pathNodeCount; iPathNode++ ) {
			Element element = (Element)node.item( iPathNode );
			if(isElementStrokeNone(element))
				continue;

			myMatrix = getMatrixFromElement(element);

			double cx=0,cy=0,rx=0,ry=0;
			if(element.hasAttribute("cx")) cx = Double.parseDouble(element.getAttribute("cx"));
			if(element.hasAttribute("cy")) cy = Double.parseDouble(element.getAttribute("cy"));
			if(element.hasAttribute("rx")) rx = Double.parseDouble(element.getAttribute("rx"));
			if(element.hasAttribute("ry")) ry = Double.parseDouble(element.getAttribute("ry"));
			jumpTo(cx+rx,cy);

			double perimeterOfAnEllipseApprox = Math.PI * 2.0 * Math.sqrt((ry*ry + rx*rx)/2.0);
			double steps = Math.max(3,perimeterOfAnEllipseApprox);
			steps = Math.min(60,steps);
			printEllipse(cx, cy, rx, ry, steps);
		}
	}

	private void printEllipse(double cx, double cy, double rx, double ry, double steps) {
		for(double i = 1; i<steps; ++i) {
			double v = (Math.PI*2.0) * (i/steps);
			double s=ry*Math.sin(v);
			double c=rx*Math.cos(v);
			moveTo(cx+c,cy+s);
		}
		moveTo(cx+rx,cy);
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
			var element = ((SVGOMPathElement)paths.item( iPath ));
			if(isElementStrokeNone(element))
				continue;

			myMatrix = getMatrixFromElement(element);

			SVGPathSegList pathList = element.getNormalizedPathSegList();
			int itemCount = pathList.getNumberOfItems();
			logger.debug("Node has {} elements.", itemCount);
			isNewSubPath =true;
			previousPathSegWasBezier=false;

			for(int i=0; i<itemCount; i++) {
				SVGPathSeg item = pathList.getItem(i);
				logger.debug(((SVGItem)item).getValueAsString());
				switch( item.getPathSegType() ) {
					case SVGPathSeg.PATHSEG_MOVETO_ABS 			-> doMoveToAbs(item);  	// M
					case SVGPathSeg.PATHSEG_MOVETO_REL 			-> doMoveToRel(item);     // m
					case SVGPathSeg.PATHSEG_LINETO_ABS 			-> doLineToAbs(item);  	// L
					case SVGPathSeg.PATHSEG_LINETO_REL 			-> doLineToRel(item);  	// l
					case SVGPathSeg.PATHSEG_CURVETO_CUBIC_ABS 	-> doCubicCurveAbs(item);	// C
					case SVGPathSeg.PATHSEG_CURVETO_CUBIC_REL 	-> doCubicCurveRel(item);	// c
					case SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_ABS -> doCubicCurveSmoothAbs(item);  // S
					case SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_REL -> doCubicCurveSmoothRel(item);  // s
					case SVGPathSeg.PATHSEG_CLOSEPATH 			-> doClosePath(); 			// Z z
					case SVGPathSeg.PATHSEG_ARC_ABS				-> doArcAbs(item);  // A
					case SVGPathSeg.PATHSEG_ARC_REL				-> doArcRel(item);  // a
					case SVGPathSeg.PATHSEG_LINETO_VERTICAL_ABS -> doLineToVerticalAbs(item);  // V
					case SVGPathSeg.PATHSEG_LINETO_VERTICAL_REL -> doLineToVerticalRel(item);  // v
					case SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_ABS -> doLineToHorizontalAbs(item);  // H
					case SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_REL -> doLineToHorizontalRel(item);  // h

					default -> throw new Exception("Found unknown SVGPathSeg type "+item.getPathSegType());
				}
			}
		}
	}

	/**
	 * Calculate the angle between vector u vector v
	 * @param u vector 1
	 * @param v vector 2
	 * @return the angle in radian
	 */
	public static double getAngle(Vector3d u, Vector3d v) {
		double cos = u.dot(v) / (u.length() * v.length());
		double result = Math.abs(Math.acos(cos));
		double sign = Math.signum(u.x*v.y - u.y*v.x);
		return (sign<0) ? -result : result;
	}

	private void doArcAbs(SVGPathSeg item) {
		var path = (SVGPathSegArcAbs) item;
		boolean fS = path.getSweepFlag();
		boolean fA = path.getLargeArcFlag();
		var r1 = path.getR1();
		var r2 = path.getR2();
		var angleDegrees = path.getAngle();
		var end = new Vector3d(path.getX(), path.getY(), 0);

		// with start (aka pathPoint), end, rx, ry find the two centers.
		drawArc(pathPoint, end, r1, r2, angleDegrees, fA, fS);
	}

	private void doArcRel(SVGPathSeg item) {
		var path = (SVGPathSegArcAbs) item;
		boolean fS = path.getSweepFlag();
		boolean fA = path.getLargeArcFlag();
		var r1 = path.getR1();
		var r2 = path.getR2();
		var angleDegrees = path.getAngle();
		var end = new Vector3d(path.getX(), path.getY(), 0);
		end.add(pathPoint);

		// with start (aka pathPoint), end, rx, ry find the two centers.
		drawArc(pathPoint, end, r1, r2, angleDegrees, fA, fS);
	}

	private void drawArc(Vector3d start, Vector3d end, double r1, double r2, double angleDegrees, boolean fA, boolean fS) {
		// translate to mid point of start-end
		var prime = new Vector3d(
				(start.x-end.x)/2.0,
				(start.y-end.y)/2.0,
				0);
		// rotate
		Matrix3d m = new Matrix3d();
		double angleRadians = Math.toRadians(angleDegrees);
		m.rotZ(-angleRadians);
		m.transform(prime);
		var x1p = prime.x;
		var y1p = prime.y;

		// Ensure radii are large enough
		// Based on http://www.w3.org/TR/SVG/implnote.html#ArcOutOfRangeParameters
		// Step (a): Ensure radii are non-zero
		// Step (b): Ensure radii are positive
		r1 = Math.abs(r1);
		r2 = Math.abs(r2);
		// Step (c): Ensure radii are large enough
		var lambda = ( (x1p * x1p) / (r1 * r1) ) + ( (y1p * y1p) / (r2 * r2) );
		if(lambda > 1) {
			r1 = Math.sqrt(lambda) * r1;
			r2 = Math.sqrt(lambda) * r2;
		}

		// Step 2: Compute (cx′, cy′)
		var sign = (fA == fS)? -1 : 1;

		double v = ((r1 * r1 * r2 * r2) - (r1 * r1 * y1p * y1p) - (r2 * r2 * x1p * x1p)) / ((r1 * r1 * y1p * y1p) + (r2 * r2 * x1p * x1p));
		Vector3d CPrime = new Vector3d(r1*y1p/r2,-r2*x1p/r1,0);
		CPrime.scale(sign * Math.sqrt(v));

		// Step 3: Compute (cx, cy) from (cx′, cy′)
		Vector3d center = new Vector3d(CPrime);
		m.transpose();
		m.transform(center);
		center.add(new Vector3d(
				(start.x+end.x)/2.0,
				(start.y+end.y)/2.0,
				0));

		// Step 4: compute start angle and sweep angle
		Vector3d s1 = new Vector3d(
				( x1p - CPrime.x)/r1,
				( y1p - CPrime.y)/r2,
				0);
		double startAngleRadians = getAngle(new Vector3d(1,0,0),s1);
		logger.debug("startAngleRadians={}",startAngleRadians);

		Vector3d s2 = new Vector3d(
				(-x1p - CPrime.x)/r1,
				(-y1p - CPrime.y)/r2,
				0);
		double sweep = getAngle(s1,s2);
		logger.debug("sweep={}",sweep);

		final double TWO_PI=Math.PI*2.0;

		double sweepRadians = sweep % TWO_PI;
		if(!fS) {
			if(sweepRadians > 0) sweepRadians -= TWO_PI;
		} else {
			if(sweepRadians < 0) sweepRadians += TWO_PI;
		}

		int steps = (int)Math.abs(sweepRadians)*5;
		logger.debug("steps={}",steps);

		logger.debug("angleDegrees={}",angleDegrees);
		double c=Math.cos(angleRadians);
		double s=Math.sin(angleRadians);
		Vector3d rx2 = new Vector3d( c, s,0 );
		Vector3d ry2 = new Vector3d( -s, c, 0 );

		Vector3d rx3 = new Vector3d();
		Vector3d ry3 = new Vector3d();
		Vector3d c2 = new Vector3d();
		double scale = fS? 1:-1;

		for(int theta = 0; theta < steps; ++theta ) {
			double r = startAngleRadians + scale * sweepRadians * ((double)theta / (double)steps);
			rx3.set(rx2);
			rx3.scale(r1 * Math.cos(r));

			ry3.set(ry2);
			ry3.scale(r2 * Math.sin(r));

			c2.set(center);
			c2.add(rx3);
			c2.add(ry3);

			moveTo(c2.x,c2.y);
		}

		pathPoint.set(end);
		moveTo(pathPoint.x,pathPoint.y);
		previousPathSegWasBezier=false;
		isNewSubPath =false;
	}

	private void doLineToVerticalAbs(SVGPathSeg item) {
		var path = (SVGPathSegLinetoVerticalAbs)item;
		pathPoint.y = path.getY();
		moveTo(pathPoint.x,pathPoint.y);
		previousPathSegWasBezier=false;
		isNewSubPath =false;
	}

	private void doLineToVerticalRel(SVGPathSeg item) {
		var path = (SVGPathSegLinetoVerticalAbs)item;
		pathPoint.y += path.getY();
		moveTo(pathPoint.x,pathPoint.y);
		previousPathSegWasBezier=false;
		isNewSubPath =false;
	}

	private void doLineToHorizontalAbs(SVGPathSeg item) {
		var path = (SVGPathSegLinetoHorizontalAbs)item;
		pathPoint.x = path.getX();
		moveTo(pathPoint.x,pathPoint.y);
		previousPathSegWasBezier=false;
		isNewSubPath =false;
	}

	private void doLineToHorizontalRel(SVGPathSeg item) {
		var path = (SVGPathSegLinetoHorizontalAbs)item;
		pathPoint.x += path.getX();
		moveTo(pathPoint.x,pathPoint.y);
		previousPathSegWasBezier=false;
		isNewSubPath =false;
	}

	private void doMoveToAbs(SVGPathSeg item) {
		var path = (SVGPathSegMovetoAbs)item;
		Vector3d p = new Vector3d(path.getX(),path.getY(),0);
		pathPoint.set(p);
		jumpTo(pathPoint.x,pathPoint.y);

		rememberNewSubPathStart(pathPoint);
		previousPathSegWasBezier=false;
	}

	private void doMoveToRel(SVGPathSeg item) {
		var path = (SVGPathSegMovetoRel)item;
		Vector3d p = new Vector3d(path.getX(),path.getY(),0);
		if(isNewSubPath) pathPoint.set(p);
		else pathPoint.add(p);
		jumpTo(pathPoint.x,pathPoint.y);

		rememberNewSubPathStart(pathPoint);
		previousPathSegWasBezier=false;
	}

	private void rememberNewSubPathStart(Vector3d p) {
		pathFirstPoint.set(p);
		isNewSubPath=false;
	}

	private void doLineToRel(SVGPathSeg item) {
		var path = (SVGPathSegLinetoRel)item;
		Vector3d p = new Vector3d(path.getX(),path.getY(),0);
		pathPoint.add(p);
		moveTo(pathPoint.x,pathPoint.y);
		previousPathSegWasBezier=false;
		isNewSubPath=false;
	}

	private void doLineToAbs(SVGPathSeg item) {
		var path = (SVGPathSegLinetoAbs)item;
		Vector3d p = new Vector3d(path.getX(),path.getY(),0);
		pathPoint.set(p);
		moveTo(pathPoint.x,pathPoint.y);
		previousPathSegWasBezier=false;
		isNewSubPath=false;
	}

	private void doCubicCurveSmoothAbs(SVGPathSeg item) {
		var path = (SVGPathSegCurvetoCubicSmoothAbs)item;

		Vector3d p0 = pathPoint;
		Vector3d p1;
		if(!previousPathSegWasBezier) {
			p1 = pathPoint;
		} else {
			// p1 is a reflection of the previous control point.
			// aka (A-B)+A aka 2A-B
			p1 = new Vector3d(smoothControlPointA);
			p1.scale(2);
			p1.sub(smoothControlPointB);
		}
		// x2,y2 is the second control point
		Vector3d p2 = new Vector3d(path.getX2(), path.getY2(), 0);
		// x3,y3 is the end point
		Vector3d p3 = new Vector3d(path.getX(), path.getY(), 0);

		drawBezier(p0, p1, p2, p3);
	}

	private void doCubicCurveSmoothRel(SVGPathSeg item) {
		var path = (SVGPathSegCurvetoCubicSmoothRel)item;

		Vector3d p0 = pathPoint;
		Vector3d p1;
		if(!previousPathSegWasBezier) {
			p1 = pathPoint;
		} else {
			// p1 is a reflection of the previous control point.
			// aka (A-B)+A aka 2A-B
			p1 = new Vector3d(smoothControlPointA);
			p1.add(smoothControlPointA);
			p1.sub(smoothControlPointB);
		}

		// x2,y2 is the second control point
		Vector3d p2 = new Vector3d(path.getX2(), path.getY2(), 0);
		// x3,y3 is the end point
		Vector3d p3 = new Vector3d(path.getX(), path.getY(),0);

		p2.add(p0);
		p3.add(p0);

		drawBezier(p0, p1, p2, p3);
	}

	private void doCubicCurveAbs(SVGPathSeg item) {
		var path = (SVGPathSegCurvetoCubicAbs) item;
		// x0,y0 is the first point
		Vector3d p0 = pathPoint;
		// x1,y1 is the first control point
		Vector3d p1 = new Vector3d(path.getX1(), path.getY1(), 0);
		// x2,y2 is the second control point
		Vector3d p2 = new Vector3d(path.getX2(), path.getY2(), 0);
		// x3,y3 is the end point
		Vector3d p3 = new Vector3d(path.getX(), path.getY(), 0);

		drawBezier(p0, p1, p2, p3);
	}

	private void doCubicCurveRel(SVGPathSeg item) {
		var path = (SVGPathSegCurvetoCubicRel)item;
		// x0,y0 is the first point
		Vector3d p0 = pathPoint;
		// x1,y1 is the first control point
		Vector3d p1 = new Vector3d(path.getX1(),path.getY1(),0);
		// x2,y2 is the second control point
		Vector3d p2 = new Vector3d(path.getX2(),path.getY2(),0);
		// x3,y3 is the end point
		Vector3d p3 = new Vector3d(path.getX(),path.getY(),0);

		p1.add(p0);
		p2.add(p0);
		p3.add(p0);

		drawBezier(p0, p1, p2, p3);
	}

	private void drawBezier(Vector3d p0,Vector3d p1,Vector3d p2,Vector3d p3) {
		Bezier b = new Bezier(
				p0.x,p0.y,
				p1.x,p1.y,
				p2.x,p2.y,
				p3.x,p3.y);
		List<Point2D> points = b.generateCurvePoints(0.1);
		for(Point2D p : points) moveTo(p.x,p.y);
		pathPoint.set(p3);
		smoothControlPointA.set(p3);
		smoothControlPointB.set(p2);
		previousPathSegWasBezier=true;
		isNewSubPath=false;
	}

	private void doClosePath() {
		logger.debug("closing to {}",pathFirstPoint);
		moveTo(pathFirstPoint.x,pathFirstPoint.y);
		isNewSubPath=true;
	}

	private void jumpTo(double x,double y) {
		Vector3d p = transform(x,y,myMatrix);
		myTurtle.jumpTo(p.x,p.y);
	}

	private void moveTo(double x,double y) {
		Vector3d p = transform(x,y,myMatrix);
		myTurtle.moveTo(p.x,p.y);
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
