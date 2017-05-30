package fusionOptics.tracer;


import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.AsciiOutForWendel;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.interfaces.Reflector;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.Box;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.surfaces.Triangle;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;
import junit.framework.TestCase;

/** Test the surface proximity and multi-hit warning system 
 *  */
public class ProximityWarningTest  extends TestCase {
	
	/** Output path for drawings and data */
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/proximityTest";
	
	private final static int nRays = 1000;
	
	private final static double wavelen = 500e-9;	// Wavelength / m 
	
	private Square sq1, sq2, sq3;
	private Optic all;
	
	@Override
	public void setUp(){
		/** Number of rays from each start position */
		
		
		//A square at x = 1.00000
		sq1 = new Square("sq1", new double[]{ 1.0, 0.0, 0.0 }, new double[]{ 1.0, 0.0, 0.0 }, new double[]{ 0.0, 0.0, 1.0 }, 0.2, 0.2, null, null, NullInterface.ideal());
		
		//A second square, which we'll move near the first one
		sq2 = new Square("sq2", new double[]{ 2.0, 0.0, 0.0 }, new double[]{ 1.0, 0.0, 0.0 }, new double[]{ 0.0, 0.0, 1.0 }, 0.2, 0.2, null, null, NullInterface.ideal());
		
		//End target, mostly just for aiming through the first two
		sq3 = new Square("sq3", new double[]{ 3.0, 0.0, 0.0 }, new double[]{ 1.0, 0.0, 0.0 }, new double[]{ 0.0, 0.0, 1.0 }, 0.2, 0.2, null, null, NullInterface.ideal());
				
		// the 'all' Optic just contains all the other elements (optics and surfaces) 
		all = new Optic("all", new Element[]{ sq1, sq2, sq3 });
	}
	
	public void testFarApart(){
		//first, big 1m gap between sq1 and sq2 
		Tracer.proximityWarningCount = 0;
		sq1.setCentre(new double[]{ 1.0, 0.0, 0.0 });
		sq2.setCentre(new double[]{ 2.0, 0.0, 0.0 });
		fireSomeRays();
		assertEquals(0, Tracer.proximityWarningCount);
	}
	
	public void testAlmostTogether(){
		//almost together, but sufficiently far apart that it should still work
		// (This will fail if you set the tolerance too small)
		Tracer.proximityWarningCount = 0;
		sq1.setCentre(new double[]{ 1.0, 0.0, 0.0 });
		sq2.setCentre(new double[]{ 1.0 + 2.0 * Tracer.reHitTolerance, 0.0, 0.0 });
		fireSomeRays();
		assertEquals(0, Tracer.proximityWarningCount);
	}
	
	public void testTogether(){
		//actually on top of each other - this should ALWAYS fail and the ray will bounce back and forth	
		Tracer.proximityWarningCount = 0;
		sq1.setCentre(new double[]{ 1.0, 0.0, 0.0 });
		sq2.setCentre(new double[]{ 1.0, 0.0, 0.0 });
		fireSomeRays();
		assertTrue(Tracer.proximityWarningCount > nRays); // should warn /at least/ once for every ray 
	}
	
	public void testSideBySide(){
		//square are in the same plane but more than tolerance apart in the centre	
		Tracer.proximityWarningCount = 0;
		sq1.setCentre(new double[]{ 1.0, -0.1 - 2*Tracer.reHitTolerance, 0.0 });
		sq2.setCentre(new double[]{ 1.0, 0.1 + 2*Tracer.reHitTolerance, 0.0 });
		fireSomeRays();
		assertEquals(0, Tracer.proximityWarningCount);
	}
	
	public void fireSomeRays(){
				
		for(int j=0;j < nRays;j++){
			
			RaySegment ray = new RaySegment(); 
					
			//start the ray off heading randomly toward the first lens primary surface
			ray.startPos = new double[]{ 0, 0, 0 };				
			ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, sq3);
				
			//ray starts with infinte elength and no end-hit. These will get filled by the tracer.
			ray.length = Double.POSITIVE_INFINITY;
			ray.startHit = null;
			ray.E0 = new double[][]{ {1,0,0,0} }; // linear polarisation not important for this example.
			ray.wavelength = wavelen;
			ray.up = new double[]{ 0, 0, 1 };
				
			// Trace it through everything until it hits an absorber or after 30 hits 
			Tracer.trace(all, ray, 30, 1e-4, true);			
		}
	}
}
