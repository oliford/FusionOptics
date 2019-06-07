package fusionOptics.surfaces;

import fusionOptics.Util;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.interfaces.Reflector;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import otherSupport.ColorMaps;
import junit.framework.TestCase;
import oneLiners.OneLiners;

/** Fires lots of random rays at a fairly bend paraboloid surface
 * and tests that the returned intersection point is near the surface.
 * 
 * VRML and SVG outputs are produced.
 */
public class ParaboloidTest extends TestCase{
	
	public void testAspheric(){
		try {
			double normal[] = Util.reNorm(new double[]{ 1.0, 0.3, 0.2 });
			double tip[] = new double[] { 0.4, 1.0, 3 };
			double up[] = Util.reNorm(Util.cross(Util.cross(normal, new double[] {0,0,1}), normal));
			double r0 = -1;
			double r1 = 1;
			double u0 = -1;
			double u1 = 1;
			double cr = 2.0;
			double cu = 2.0;
			
			//Paraboloid paraboloid = new Paraboloid("Paraboloid", tip, normal, up, cu, cr, new double[] { u0, r0, u1, r1}, Reflector.ideal());			
			//Sphere focusSphere = new Sphere("focusSphere", paraboloid.getFocus(), 0.050, Absorber.ideal());
			
			Sphere focusSphere = new Sphere("focusSphere", new double[] { 3.4, 1.0, 3.0 }, 0.100, Absorber.ideal());
			Paraboloid paraboloid = new Paraboloid("Paraboloid", new double[] { 0,0,0}, focusSphere.getCentre(), normal, 0.500, null, null, Reflector.ideal());
			
			OneLiners.dumpArray(focusSphere.getCentre());
			OneLiners.dumpArray(paraboloid.getFocus());
			
			Optic all = new Optic("all", new Element[] { paraboloid, focusSphere });
			
			int n = 50;
			
			SVGRayDrawing svgOut = new SVGRayDrawing("/tmp/rayTracing/tests/paraboloid.txt", Util.getBoundingBox(paraboloid), true);
			VRMLDrawer vrmlOut = new VRMLDrawer("/tmp/rayTracing/tests/paraboloid.vrml");
			vrmlOut.setDrawIntersectionNormals(true);
			vrmlOut.setDrawPolarisationFrames(false);
			vrmlOut.setSmallLineLength(0.01);
			vrmlOut.setSkipRays(0);
			
			
			double col[][] = ColorMaps.jet(n);
			svgOut.generateLineStyles(col, 0.001);
			
			//double start[] = new double[]{ 1000, 0, 0 };
			double start[] = Util.plus(paraboloid.getCentre(), Util.mul(paraboloid.getNormal(), 1000));
			double dir[] = Util.reNorm(Util.mul(normal, -1.0));
			
			for(int i=0; i < n; i++){
				
				 
				
				//double dir[];
				start = Tracer.generateRandomInfinityRayStartPos(dir.clone(), paraboloid, 10.0);
				//dir = Tracer.generateRandomRayTowardSurface(start, paraboloid);
				
				RaySegment ray = new RaySegment(start.clone(), dir.clone());
				
				Tracer.trace(all, ray, 10, 0.01, false);
								
				if(ray.endHit != null && Double.isNaN(Util.length(ray.endHit.normal)) ){
						ray.length = Double.POSITIVE_INFINITY;
						Tracer.trace(paraboloid, ray, 10, 0.01, false);
					
						fail("Hit surface, but normal is NaN");
				}
				
				svgOut.drawRay(ray, i);
				vrmlOut.drawRay(ray, col[i]);
				
				Pol.recoverAll();
				
			}
			
			vrmlOut.drawElement(all);
			svgOut.drawElement(all);	
			
			//Sphere boundingSphere = new Sphere("ParaboloidBoundary", paraboloid.getBoundarySphereCentre(), paraboloid.getBoundarySphereRadius(), NullInterface.ideal());
			//vrmlOut.drawElement(boundingSphere);
			//svgOut.drawElement(boundingSphere);	
			
			vrmlOut.destroy();		
			svgOut.destroy();

		}catch(RuntimeException e) {
			e.printStackTrace();
			throw e;
		}
	}

}
