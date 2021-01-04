package fusionOptics.examples.tracer;

import net.jafama.FastMath;
import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.lenses.EdmundOptics50mmAspheric;
import fusionOptics.lenses.Nikon50mmF11;
import fusionOptics.materials.BK7;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.SimplePlanarAsphericLens;
import fusionOptics.optics.SimplePlanarConvexLens;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import otherSupport.ColorMaps;
import otherSupport.RandomManager;
import oneLiners.OneLiners;
import svg.SVGSplitView3D;

/** Imaging example of aspheric lens.
 * 
 * Uses coefficients from Edmund Optics to check that our coefficients  
 * 
 * 
 * @author oliford  */
public class AsphericImagingExample {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/asphericImaging";	
	final static int nRays = 50;	// number of rays for each point on grid
	final static double angY[] = new double[]{ 0 }; //OneLiners.linSpace(-5.0 * Math.PI / 180, 5.0 * Math.PI / 180, 5);
	final static double angZ[] = OneLiners.linSpace(-5.0 * Math.PI / 180, 5.0 * Math.PI / 180, 5);
		
	
	final static double wavelen = EdmundOptics50mmAspheric.designWavelen;
		
	// rho_max = 53.78288 ???	
	public static void main(String[] args) {
		
		EdmundOptics50mmAspheric lens = new EdmundOptics50mmAspheric(new double[]{0,0,0}, new double[]{1,0,0});
		
		System.out.println("n("+ wavelen/1e-9 +"nm) = " + EdmundOptics50mmAspheric.mat.getRefractiveIndex(0, wavelen, 300) + "\tArticle says 1.51680"); //
		
		double imgPlanePos = lens.getPlanarSurface().getCentre()[0] + EdmundOptics50mmAspheric.backFocalDistance;
				
		Square imgPlane = new Square("imgPlane", new double[]{ imgPlanePos, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.040, 0.040, Absorber.ideal());
		Optic all = new Optic("all", new Element[]{ lens, imgPlane });
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/imgTest.vrml", 0.002);
		vrmlOut.setDrawIntersectionNormals(true);
		vrmlOut.setDrawPolarisationFrames(false);
		
		double col[][] = ColorMaps.alternating2D2x2(angY.length, angZ.length); //alternating colours over point grid		
		SVGRayDrawing svgOut = new SVGRayDrawing(outPath + "/imgTest", new double[]{ 0, -1, -1, 2, 1, 1 }, true );
		svgOut.generateLineStyles(col, 0.00005);
		svgOut.setSmallLineLength(0.002);
		svgOut.setDrawIntersectionNormals(true);
		
		for(int iY=0; iY < angY.length; iY++) {
			for(int iZ=0; iZ < angZ.length; iZ++) {
				for(int i=0; i < nRays; i++) {					
					//RaySegment ray = new RaySegment(new double[]{ 0, y[iY] + 0.5, z[iZ] }, lens);
					RaySegment ray = new RaySegment();
					
					double pointAtLens[] = new double[]{ 
							1.0, 
							RandomManager.instance().nextUniform(-EdmundOptics50mmAspheric.diameter/2, EdmundOptics50mmAspheric.diameter/2), 
							RandomManager.instance().nextUniform(-EdmundOptics50mmAspheric.diameter/2, EdmundOptics50mmAspheric.diameter/2) 
						};
					ray.dir = Util.reNorm(new double[]{ FastMath.cos(angZ[iZ]) *  FastMath.cos(angY[iY]), FastMath.sin(angY[iY]), FastMath.sin(angZ[iZ]) });
					ray.startPos = Util.minus(pointAtLens, ray.dir);
					
					ray.wavelength = wavelen;
					ray.up = new double[]{ 0, 1, 0 };
					ray.E0 = new double[][]{{1,0,0,0}}; //lin.pol. in 'up' dir
					
					Tracer.trace(all, ray, 30, 0.01, true);
					
					vrmlOut.drawRay(ray, col[iY*angZ.length+iZ]);
					svgOut.drawRay(ray, iZ);
					Pol.recoverAll();
				}
			}
		}
		
		vrmlOut.drawOptic(all);
		vrmlOut.destroy();
		
		svgOut.drawElement(all);
		svgOut.destroy();
	}		
}
