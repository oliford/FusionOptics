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

/** Shortest possible code to produce a nice imaging VRML output
 * @author oliford  */
public class SuccinctImagingVRMLExample {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/succinctImaging";	
	final static int nRays = 100;	// number of rays for each point on grid
	final static double y[] = OneLiners.linSpace(-0.2, 0.2, 6); //object grid Y positions
	final static double z[] = OneLiners.linSpace(-0.2, 0.2, 6); //object grid Z positions
	
	public static void main(String[] args) {	
		
		Nikon50mmF11 lens = new Nikon50mmF11(new double[]{ 1, 0, 0 });		
		Square imgPlane = new Square("imgPlane", new double[]{ 1.0526, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.040, 0.040, Absorber.ideal());
		Optic all = new Optic("all", new Element[]{ lens, imgPlane });
		
		double col[][] = ColorMaps.alternating2D2x2(y.length, z.length); //alternating colours over point grid		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/imgTest.vrml", 0.0001);
		vrmlOut.setDrawPolarisationFrames(false);
		
		for(int iZ=0; iZ < z.length; iZ++) {
			for(int iY=0; iY < y.length; iY++) {
				for(int i=0; i < nRays; i++) {					
					RaySegment ray = new RaySegment(new double[]{ 0, y[iY], z[iZ] }, lens);					
					Tracer.trace(all, ray, 30, 0.01, true);					
					vrmlOut.drawRay(ray, col[iY*z.length+iZ]);
					Pol.recoverAll();
				}				
			}
		}
		vrmlOut.drawOptic(all);
		vrmlOut.destroy();
	}		
}
