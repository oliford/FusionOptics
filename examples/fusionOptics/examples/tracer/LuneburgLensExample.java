package fusionOptics.examples.tracer;

import net.jafama.FastMath;
import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.lenses.Nikon50mmF11;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.surfaces.Sphere;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import otherSupport.ColorMaps;
import otherSupport.RandomManager;
import oneLiners.OneLiners;
import svg.SVGSplitView3D;

/** The Lüneburg lens - 
 * A sphere of glass with graded refractive index ranging as sqrt(2 - r²/R²)
 * Any parallel input light should focus exactly to a single point on the opposite surface of the sphere. 
 * 
 * It's implemented here as a series of spherical surfaces with discrete indices.
 * Unfortunately, this doesn't approximate very well and takes a lot of surfaces (~1000s) to focus well.
 * It improves with increasing number of surfaces, so presumably it limits to perfect for infinite shells. 
 * 
 * @author oliford  
 */
public class LuneburgLensExample {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/luneBurg";	
	final static int nRays = 50;	
	final static int nSpheres = 2000;
	final static double angle[] = OneLiners.linSpace(-70*Math.PI/180, 70 * Math.PI/180, 5);
	
	public static void main(String[] args) {
		double maxR = 0.050;
		double radii[] = OneLiners.linSpace(maxR, 0.001, nSpheres);
		double n[] = new double[nSpheres];
		Sphere spheres[] = new Sphere[radii.length];
		
		Optic luneburgLens = new Optic("luneburgLens");
		Medium lastMedium = null; 
		double dR = radii[1]-radii[0];
		for(int i=0; i < radii.length; i++){
			double avgR = radii[i] + dR/ 2;
			n[i] = FastMath.pow(2 - avgR*avgR/maxR/maxR, 0.5);
			
			Medium nextMedium = new Medium(new IsotropicFixedIndexGlass(n[i]));
			spheres[i] = new Sphere("sphere" + i, new double[]{ 1, 0, 0}, radii[i], nextMedium, lastMedium, IsoIsoInterface.ideal());
			spheres[i].nRingsPerHemisphere = 30;
			spheres[i].nSectors = 15;
			spheres[i].nPointsInSector = 5;
			luneburgLens.addElement(spheres[i]);
			lastMedium = nextMedium;
			
			System.out.println("Sphere " + i + ", <R>=" + avgR + ", n=" + n[i]);
		}
			
		Square imgPlane = new Square("imgPlane", new double[]{ 1.100, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.100, 0.100, Absorber.ideal());
		Optic all = new Optic("all", new Element[]{ luneburgLens, imgPlane });
		
		double col[][] = ColorMaps.jet(angle.length);
		SVGRayDrawing svgOut = new SVGRayDrawing(outPath + "/imgTest", new double[]{ 0, -1, -1, 2, 1, 1 }, true );
		svgOut.generateLineStyles(col, 0.0001);
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/imgTest.vrml", 0.0001);
		vrmlOut.setDrawPolarisationFrames(false);
		
		
		for(int iA=0; iA < angle.length; iA++) { //for each input angle
			for(int i=0; i < nRays; i++) { //for each ray
				double dir[] = Util.reNorm(new double[]{ FastMath.cos(angle[iA]), 0, FastMath.sin(angle[iA]) });
				
				double start[] = Tracer.generateRandomInfinityRayStartPos(dir, luneburgLens, 0.500); //
				//double start[] = new double[]{ 0.5, 0, RandomManager.instance().nextUniform(-maxR, maxR) };
				
				RaySegment ray = new RaySegment(start, dir);
				Tracer.trace(all, ray, 30000, 0.01, true);
				
				svgOut.drawRay(ray, iA);			
				vrmlOut.drawRay(ray, col[iA]);
				Pol.recoverAll();
				
				System.out.println("Angle " +iA + " / " + angle.length + ", ray " + i + " / " + nRays);
			}
		}
	
		//svgOut.drawElement(all);
		svgOut.drawElement(imgPlane);
		svgOut.drawElement(spheres[0]);
		svgOut.destroy();
		//vrmlOut.drawOptic(all);
		vrmlOut.drawElement(spheres[0]);
		vrmlOut.destroy();	
	}		
}
