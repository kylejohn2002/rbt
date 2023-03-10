package com.marginallyclever.makelangelo.makeart.imagefilter;

import com.marginallyclever.makelangelo.makeart.TransformedImage;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.awt.image.BufferedImage;

/**
 * Converts an image to N shades of grey.
 * @author Dan Royer
 */
public class Filter_Greyscale extends ImageFilter {
	private static final Logger logger = LoggerFactory.getLogger(Filter_Greyscale.class);
	private double levels = 2;
	private int mode = 1;

	public Filter_Greyscale(int _levels) {
		levels = (double) _levels;
	}

	public TransformedImage filter(TransformedImage img) {
		return switch (mode) {
			case 0 -> filterLevels(img);
			case 1 -> filterTone(img);
			case 2 -> filterSimple(img);
			default -> null;
		};
	}

	protected TransformedImage filterLevels(TransformedImage img) {
		int h = img.getSourceImage().getHeight();
		int w = img.getSourceImage().getWidth();
		int x, y, i;

		double max_intensity = -1000;
		double min_intensity = 1000;

		BufferedImage bi = img.getSourceImage();
		for (y = 0; y < h; ++y) {
			for (x = 0; x < w; ++x) {
				i = decode32bit(bi.getRGB(x, y));
				if (max_intensity < i) max_intensity = i;
				if (min_intensity > i) min_intensity = i;
			}
		}
		double intensity_range = max_intensity - min_intensity;

		double ilevels = 1;
		if (levels != 0)
			ilevels = 1.0 / levels;

		// logger.debug("min_intensity="+min_intensity);
		// logger.debug("max_intensity="+max_intensity);
		// logger.debug("levels="+levels);
		// logger.debug("inverse="+ilevels);

		double pixel;

		TransformedImage after = new TransformedImage(img);
		BufferedImage afterBI = after.getSourceImage();

		for (y = 0; y < h; ++y) {
			for (x = 0; x < w; ++x) {
				pixel = decode32bit(bi.getRGB(x, y));

				double a = (pixel - min_intensity) / intensity_range;
				double c = a * levels * ilevels;
				int b = (int) Math.max(Math.min(c * 255.0, 255), 0);
				// if(b==255) logger.debug(x+"\t"+y+"\t"+i+"\t"+b);
				afterBI.setRGB(x, y, ImageFilter.encode32bit(b));
			}
		}

		return after;
	}

	// accepts and returns a number between 0 and 255, inclusive.
	private double toneControl(double v) {
		v /= 255.0;
		v = 0.017 * Math.exp(3.29 * v) + 0.005 * Math.exp(7.27 * v);
		return Math.min(1, Math.max(0, v)) * 255.0;
	}

	public TransformedImage filterTone(TransformedImage img) {
		int h = img.getSourceImage().getHeight();
		int w = img.getSourceImage().getWidth();

		BufferedImage bi = img.getSourceImage();
		TransformedImage after = new TransformedImage(img);
		BufferedImage afterBI = after.getSourceImage();

		int x, y;
		for (y = 0; y < h; ++y) {
			for (x = 0; x < w; ++x) {
				double pixel = decode32bit(bi.getRGB(x, y));
				double v2 = toneControl(pixel);
				int rgb = (int) Math.min(255, Math.max(0, v2));
				afterBI.setRGB(x, y, ImageFilter.encode32bit(rgb));
			}
		}
		return after;
	}

	public TransformedImage filterSimple(TransformedImage img) {
		int h = img.getSourceImage().getHeight();
		int w = img.getSourceImage().getWidth();

		BufferedImage bi = img.getSourceImage();
		TransformedImage after = new TransformedImage(img);
		BufferedImage afterBI = after.getSourceImage();

		int x, y;
		for (y = 0; y < h; ++y) {
			for (x = 0; x < w; ++x) {
				double pixel = decode32bit(bi.getRGB(x, y));
				int rgb = (int) Math.min(255, Math.max(0, pixel));
				afterBI.setRGB(x, y, ImageFilter.encode32bit(rgb));
			}
		}
		return after;
	}
}