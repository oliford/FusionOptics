package fusionOptics.tracer;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.lenses.Nikon50mmF11;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import otherSupport.ColorMaps;
import oneLiners.OneLiners;
import svg.SVGSplitView3D;

/** Shortest possible code to produce a nice imaging SVG
 * @author oliford  */
public class SuccinctImagingSVGExample {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/succinctImaging";	
	final static int nRays = 500;	
	final static double z[] = OneLiners.linSpace(-0.2, 0.2, 6);
	
	public static void main(String[] args) {
		Nikon50mmF11 lens = new Nikon50mmF11(new double[]{ 1, 0, 0 });		
		Square imgPlane = new Square("imgPlane", new double[]{ 1.0526, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.040, 0.040, Absorber.ideal());
		Optic all = new Optic("all", new Element[]{ lens, imgPlane });
		
		double col[][] = ColorMaps.jet(z.length);
		SVGRayDrawing svgOut = new SVGRayDrawing(outPath + "/imgTest", new double[]{ 0, -1, -1, 2, 1, 1 }, true );
		svgOut.generateLineStyles(col, 0.0002);
		
		for(int iZ=0; iZ < z.length; iZ++) {
			for(int i=0; i < nRays; i++) {					
				RaySegment ray = new RaySegment(new double[]{ 0, 0, z[iZ] }, lens);					
				Tracer.trace(all, ray, 30, 0.01, true);			
				svgOut.drawRay(ray, iZ);
				Pol.recoverAll();
			}				
		}
	
		svgOut.drawElement(all);
		svgOut.destroy();
	}		
}
