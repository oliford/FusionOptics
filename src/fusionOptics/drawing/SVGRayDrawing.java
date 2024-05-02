package fusionOptics.drawing;

import net.jafama.FastMath;

import java.util.ArrayList;
import java.util.List;

import fusionOptics.types.Element;
import fusionOptics.types.Optic;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;

import uk.co.oliford.jolu.OneLiners;
import uk.co.oliford.svg.SVGSplitView3D;

/** Simple ray drawing to SVG, without polarisation information. */
public class SVGRayDrawing extends LineDrawer {
	private SVGSplitView3D svg3D;
	private int currentColour;
	private boolean drawingOptic = false;
	
	/**
	 * @param fileNamePrefix
	 * @param bbox { x0, y0, z0, x1, y1, z1 };
	 * @param drawWholePaths If true, each ray path is drawn from start to end. If false, branches are drawn from branch point only 
	 * @param rotMat rotation matrix to apply to each coordinate before drawing
	 */
	public SVGRayDrawing(String fileNamePrefix, double[] bbox, boolean drawWholePaths, double rotMat[][]) {
		super(drawWholePaths, FastMath.sqrt(
				FastMath.pow2(bbox[3] - bbox[0]) +
				FastMath.pow2(bbox[4] - bbox[1]) +
				FastMath.pow2(bbox[5] - bbox[2])) / 200
			);
		
		OneLiners.makePath(fileNamePrefix);
		svg3D = new SVGSplitView3D(fileNamePrefix, bbox, rotMat);
		
	}
	
	/** @param fileNamePrefix
	 * @param bbox { x0, y0, z0, x1, y1, z1 };
	 * @param drawWholePaths If true, each ray path is drawn from start to end. If false, branches are drawn from branch point only 
	 */
	public SVGRayDrawing(String fileNamePrefix, double[] bbox, boolean drawWholePaths) {
		this(fileNamePrefix, bbox, drawWholePaths, null);
	}
		
	/** Sets the colour table and line width. Must be called before first draw command 
	 * @param colourTable		Colour tables [index][r/g/b]
	 * @param lineWidth			Width of lines to draw rays and optics.
	 */
	public void generateLineStyles(double colourTable[][], double lineWidth) {
		generateLineStyles(colourTable, lineWidth, lineWidth);
	}
	
	/** Sets the colour table and line widths. Must be called before first draw command 
	 * 
	 * @param colourTable		Colour tables [index][r/g/b]
	 * @param rayWidth			Width of lines to draw rays
	 * @param opticLineWidth	Width of lines to draw optics
	 */
	public void generateLineStyles(double colourTable[][], double rayWidth, double opticLineWidth) {
	    for(int i=0; i<colourTable.length; i++){
			String cStr = "#" +
					((colourTable[i][0]*255) < 16 ? "0" : "") + Integer.toHexString((int)(colourTable[i][0]*255)) +
					((colourTable[i][1]*255) < 16 ? "0" : "") + Integer.toHexString((int)(colourTable[i][1]*255)) + 
					((colourTable[i][2]*255) < 16 ? "0" : "") + Integer.toHexString((int)(colourTable[i][2]*255)) ;
			svg3D.addLineStyle("c"+i, "none", rayWidth, cStr);
		}			
		svg3D.addLineStyle("optic", "none", opticLineWidth, "green");
	}
	
	/** Draws a ray with the given color. 
	 * 
	 * @param ray
	 * @param colourNumber	Index into the color table given in generateLineStyles()
	 */
	public void drawRay(RaySegment ray, int colourNumber) {
		this.currentColour = colourNumber;
		super.drawRay(ray);
	}

	@Override
	/** Draw a line from a List of points double[x/y/z] */
	public void drawSinglePath(List<double[]> path) {
		svg3D.addLine(path, drawingOptic ? "optic" : "c" + currentColour);
	}
	
	@Override
	/** Draw a line from an array of points double[index][x/y/z] */
	public void drawSinglePath(double[][] path) {
		svg3D.addLine(path, drawingOptic ? "optic" : "c" + currentColour);
	}
	
	/** Draw a line from an array of points double[index][x/y/z] */
	public void drawSinglePath(double[][] path, String style) {
		svg3D.addLine(path, style);
	}
	
	@Override
	/** Draw an optic element */
	public void drawElement(Element elem) {
		drawingOptic = true;
		super.drawElement(elem);
		drawingOptic = false;
	}

	/** Finalise and close the SVG */
	public void destroy() {
		//finish the SVG off
		svg3D.destroy();
	}

	public SVGSplitView3D getSVG3D() { return svg3D; }
	@Override
	public void startGroup(String name) { svg3D.startGroup(name); }
	@Override
	public void endGroup() { svg3D.endGroup(); }

	/** Sets the colour used to draw rays.
	 * @param colour	Integer index into the colour table set in SVGRayDrawing.generateLineStyles()
	 */
	public void setRayColour(int colour){ this.currentColour = colour; }

	public void setPrecision(int precision){ svg3D.setPrecision(precision); }
	
}
