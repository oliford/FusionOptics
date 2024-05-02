package fusionOptics.surfaces;

import fusionOptics.Util;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Reflector;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Intersection;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import junit.framework.TestCase;
import uk.co.oliford.jolu.ColorMaps;

/** Fires lots of random rays at a fairly bends aspheric surface
 * and tests that the returned intersection point is near the surface.
 * 
 * VRML and SVG outputs are produced.
 */
public class AsphericTest extends TestCase{
	
	public void testAspheric(){
		
		Aspheric aspheric = new Aspheric("aspheric", 
				new double[]{ 0, 0, 0 }, 
				Util.reNorm(new double[]{ 1.0, 0.2, 0.4 }),
				3.0, //rad of curv.
				1.0, //rim rad
				0.0, //conic constant (NOT SUPPORTED in surface finder)
				new double[]{ 0.3, -0.97, 0.6 }, //r^2, r^4, r^6 coeffs 
				//new double[]{ 0.0, 0.0, 0.0 }, //r^2, r^4, r^6 coeffs
				null, 
				null, 
				Reflector.ideal());
		/*Dish spheric = new Dish("spheric", 
				new double[]{ 1, 0, 0 }, 
				Util.reNorm(new double[]{ 1.0, 0.2, 0.4 }),
				3.0, //rad of curv.
				1.0, //rim rad
				null, 
				null, 
				Reflector.ideal());
		//*/
		
		int n = 2000;
		
		SVGRayDrawing svgOut = new SVGRayDrawing("/tmp/rayTracing/tests/aspheric.txt", Util.getBoundingBox(aspheric), true);
		VRMLDrawer vrmlOut = new VRMLDrawer("/tmp/rayTracing/tests/aspheric.vrml", 0.05);
		vrmlOut.setDrawIntersectionNormals(true);
		vrmlOut.setDrawPolarisationFrames(false);
		vrmlOut.setSkipRays(9);
		
		double col[][] = ColorMaps.jet(n);
		svgOut.generateLineStyles(col, 0.001);
		
		for(int i=3; i < n; i++){
			
			RaySegment ray = new RaySegment(new double[]{ -1, 0, 0 }, new double[]{ 1, 0, 0 });
			
			ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, aspheric);
			
			Tracer.trace(aspheric, ray, 10, 0.01, false);
			
			if(ray.endHit != null && Double.isNaN(Util.length(ray.endHit.normal)) ){
					ray.length = Double.POSITIVE_INFINITY;
					Tracer.trace(aspheric, ray, 10, 0.01, false);
				
					fail("Hit surface, but normal is NaN");
			}
			
			svgOut.drawRay(ray, i);
			vrmlOut.drawRay(ray, col[i]);
			
			Pol.recoverAll();
			
		}
		
		vrmlOut.drawElement(aspheric);
		svgOut.drawElement(aspheric);		
		//vrmlOut.drawElement(spheric);
		
		vrmlOut.destroy();		
		svgOut.destroy();
	
		
		System.out.println("Intersections dropped out to brute force:" + aspheric.surfaceBruteForces + " / " + n);
		System.out.println("Intersections failed:" + aspheric.surfaceFailures);
		if(aspheric.surfaceFailures > 0){				
			fail("Failed to find surface " + aspheric.surfaceFailures + " times");
		}
	}

}
