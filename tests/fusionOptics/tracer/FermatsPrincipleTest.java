package fusionOptics.tracer;

import java.util.List;

import net.jafama.FastMath;
import uk.co.oliford.jolu.ColorMaps;
import uk.co.oliford.jolu.OneLiners;
import junit.framework.TestCase;
import fusionOptics.MinervaOpticsSettings;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.materials.BK7;
import fusionOptics.optics.SimplePlanarConvexLens;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;

/** Test Fermat's principle - The total optical path length (wave count + phase) should be the same for
 * all rays imaged onto the same point.
 * 
 * @author oliford  */
public class FermatsPrincipleTest extends TestCase {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/fermatTest";
	final static int nRays = 10;	
	final static double a[] = OneLiners.linSpace(-5.0* Math.PI/180, 5.0* Math.PI/180, 6) ; //+/- 20mm
	final static double wavelen = 650e-9;
	final static double focalLength = 0.300; //80mm 
	
	public void testNikonLens() {
		
		//Nikon50mmF11 lens = new Nikon50mmF11(new double[]{ 1, 0, 0 });
		Medium lensMedium = new Medium(new BK7());
		SimplePlanarConvexLens lens = SimplePlanarConvexLens.fromFocalLengthAndCentreThickness("lens", 
												new double[]{ 0.050, 0, 0 }, //centre at x=50mm
												new double[]{ 1, 0, 0 }, //normal
												0.020, //radius = 20mm
												focalLength, //focalLength = 50mm
												0.010, //centreThickness = 10mm
												lensMedium, IsoIsoInterface.ideal(), wavelen);
		
		Iris iris = new Iris("iris", new double[]{ 0.050, 0, 0 }, new double[]{1,0,0}, 0.100, 0.0198, Absorber.ideal());
				
		Square imgPlane = new Square("imgPlane", new double[]{ 0.050 + focalLength - 0.010, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.100, 0.100, Absorber.ideal());
		Optic all = new Optic("all", new Element[]{ iris, lens, imgPlane });
		
		//AutoFocus autoFocus = new AutoFocus(all, imgPlane, wavelen);
		//autoFocus.setMoveableElement(imgPlane, 0.100);
		//autoFocus.setParallelRays(lens.getFrontSurface(), new double[]{1,0,0}, 5*Math.PI/180, 5);
		//autoFocus.go(100);
		
		double col[][] = ColorMaps.jet(a.length);
		SVGRayDrawing svgOut = new SVGRayDrawing(outPath + "/fermat", new double[]{ 0, -1, -1, 2, 1, 1 }, true );
		svgOut.generateLineStyles(col, 0.0002);
		
		for(int iA=0; iA < a.length; iA++) {
			double dir[] = { FastMath.cos(a[iA]), 0, FastMath.sin(a[iA]) };
			
			
			for(int i=0; i < nRays; i++) {					
				//RaySegment ray = new RaySegment(new double[]{ 0, 0, z[iZ] }, lens);
									 
				double startPos[] = Tracer.generateRandomInfinityRayStartPos(dir, lens.getFrontSurface(), 0.050);
				RaySegment ray = new RaySegment(startPos, dir);
				Tracer.trace(all, ray, 30, 0.01, true);			
				svgOut.drawRay(ray, iA);
				
				List<Intersection> hits = ray.getIntersections(imgPlane);
				double rayOPL = Double.NaN;
				double calcOPL = Double.NaN;
				if(hits.size() > 0){
					Intersection hit = hits.get(0);
					rayOPL = FastMath.atan2(hit.incidentRay.E1[0][1], hit.incidentRay.E1[0][0]) / 2*Math.PI;
					rayOPL += hit.incidentRay.nWaves;
				
				
					calcOPL = 0;
					while(hit.incidentRay != null && hit.incidentRay.startHit != null){
						double n = 1.0;
						if(!Double.isNaN(hit.incidentRay.raySpecificRefractiveIndex) && hit.incidentRay.raySpecificRefractiveIndex != 0){
							n = hit.incidentRay.raySpecificRefractiveIndex;
						}else if(hit.incidentRay.medium != null){
							n = hit.incidentRay.medium.getRefractiveIndex(0, hit.incidentRay.wavelength);
						}
						calcOPL += hit.incidentRay.length * n;						
						
						hit = hit.incidentRay.startHit;
					}
				}
						
				System.out.println(iA + "\t" + i + "\t" + rayOPL + "\t" + (calcOPL / wavelen));
				
				Pol.recoverAll();
				
			}				
		}
	
		svgOut.drawElement(all);
		svgOut.destroy();
		
		fail("No, just no :(");
	}		
}
