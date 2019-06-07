package fusionOptics.drawing;

import net.jafama.FastMath;

import java.util.ArrayList;
import java.util.List;

import fusionOptics.types.Element;
import fusionOptics.types.Optic;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;

import algorithmrepository.Algorithms;

import oneLiners.OneLiners;

import svg.SVGLineWriter;
import svg.SVGSplitView3D;

/** Special version of the SVG drawing, for drawing cylindrical (e.g. poloidal for Tokamaks) projections */
public class SVGCylindricalProjection extends LineDrawer {
	private SVGLineWriter svg;
	double smallLineLength;
	private int currentColour;
	
	/**
	 * @param fileNamePrefix
	 * @param bbox { x0, y0, z0, x1, y1, z1 };
	 * @param drawWholePaths If true, each ray path is drawn from start to end. If false, branches are drawn from branch point only 
	 * @param rotMat rotation matrix to apply to each coordinate before drawing
	 */
	public SVGCylindricalProjection(String fileName, double maxR, double Z0, double Z1, boolean drawWholePaths) {
		super(drawWholePaths, FastMath.sqrt(
				FastMath.pow2(maxR) +
				FastMath.pow2(Z1 - Z0)) / 200
			);
		
		svg = new SVGLineWriter(fileName, new double[]{ -maxR, Z0, maxR, Z1});
		
	}
	
	public void generateLineStyles(double colourTable[][], double lineWidth) {
		generateLineStyles(colourTable, lineWidth, lineWidth);
	}
	
	/** Sets the colour table, must be called before drawing */
	public void generateLineStyles(double colourTable[][], double lineWidth, double opticLineWidth) {
		for(int i=0; i<colourTable.length; i++){
			String cStr = "#" +
					((colourTable[i][0]*255) < 16 ? "0" : "") + Integer.toHexString((int)(colourTable[i][0]*255)) +
					((colourTable[i][1]*255) < 16 ? "0" : "") + Integer.toHexString((int)(colourTable[i][1]*255)) + 
					((colourTable[i][2]*255) < 16 ? "0" : "") + Integer.toHexString((int)(colourTable[i][2]*255)) ;
			svg.addLineStyle("c"+i, "none", lineWidth, cStr);
			svg.addLineStyle("optic", "none", opticLineWidth, "green");
		}			
	}
	
	public void drawRay(RaySegment ray, int colourNumber) {
		this.currentColour = colourNumber;
		super.drawRay(ray);
	}

	@Override
	public void drawSinglePath(List<double[]> path) {
		addLine(path, "c" + currentColour);
	}
	

	@Override
	public void drawSinglePath(double[][] path) {
		addLine(path, "c" + currentColour);
	}
	
	public void addLines(List<double[][]> lines, String style) {
        for(double line[][] : lines){
            addLine(line[0], line[1], line[2], style);
        }
    }
	
	 public void addLine(double posXYZ[][], String style) {
	        addLine(posXYZ[0], posXYZ[1], posXYZ[2], style);
	   }
	
	 public void addLine(List<double[]> posXYZ, String style) {
	        double x[] = new double[posXYZ.size()];
	        double y[] = new double[posXYZ.size()];
	        double z[] = new double[posXYZ.size()];
	        
	        for(int i=0; i < posXYZ.size(); i++){
	            double p[] = posXYZ.get(i);
	            x[i] = p[0];
	            y[i] = p[1];
	            z[i] = p[2];
	        }
	        
	        addLine(x, y, z, style);
	    }
	    
	 
	 public void addLine(double x[], double y[], double z[], String style){
		 ArrayList<double[]> posRZ = new ArrayList<double[]>(2*x.length);
		 
		 double Rl = FastMath.sqrt(x[0]*x[0] + y[0]*y[0]);
		 posRZ.add(new double[]{ Rl, z[0] });
		 
		 for(int i=1; i < x.length; i++){
			 double R = FastMath.sqrt(x[i]*x[i] + y[i]*y[i]);
			 double minR = FastMath.max(FastMath.min(R, Rl), 0.01);
			 double minL = minR * 0.1; 
				 
			 double l = FastMath.sqrt(FastMath.pow2(y[i] - y[i-1]) + FastMath.pow2(x[i] - x[i-1]));
			 int n = (int)((l / minL) + 0.5);
			 if(n <= 1){
				 posRZ.add(new double[]{ R, z[i] });
				 
			 }else{
				 double dl = l / (n - 1);
				 for(int j=1; j < n; j++){
					 double xSD = x[i-1] + j * dl * (x[i] - x[i-1]);
					 double ySD = y[i-1] + j * dl * (y[i] - y[i-1]);					 
					 double zSD = z[i-1] + j * dl * (z[i] - z[i-1]);					 
					 double RSD = FastMath.sqrt(xSD*xSD + ySD*ySD);						
					 posRZ.add(new double[]{ RSD, zSD});					 
				 }
			 }
			 Rl = R;
		 }
		 
		 double R[] = new double[posRZ.size()];
		 double Z[] = new double[R.length];
		 for(int i=0; i < R.length; i++){
			 double p[] = posRZ.get(i);
			 R[i] = p[0];
			 Z[i] = p[1];
		 }
		 svg.addLine(R, Z, style);
	 }
	
	public void destroy() {
		//finish the SVG off
		svg.destroy();
	}

	public SVGLineWriter getSVG() { return svg; }

	@Override
	public void startGroup(String name) { svg.startGroup(name); }
	@Override
	public void endGroup() { svg.endGroup(); }

}
