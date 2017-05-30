package fusionOptics.tracer;

import fusionOptics.Util;
import fusionOptics.interfaces.Absorber;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import junit.framework.TestCase;


public class BoundingSphereTest extends TestCase {
	public static class BSTester extends Optic {
		public BSTester() {	super("bsTester"); }
		@Override
		public double getBoundarySphereRadius() { return 1; }

		@Override
		public double[] getBoundarySphereCentre() { return new double[]{ 2, 0, 0 }; } 
	}
	
	
	public void testBoundingSphereHitTest() {
		
		BSTester bsTester = new BSTester();
	
		RaySegment ray = new RaySegment();
		
		//enter and exit
		ray.startPos = new double[]{ 0, 0, 0 };
		ray.dir = new double[]{ 1, 0, 0 };
		ray.length = 5;
		assertEquals(bsTester.testBoundingSphere(ray), true);

		//completly outside
		ray.startPos = new double[]{ 0, 0, 8 };
		ray.dir = new double[]{ 1, 0, 0 };
		ray.length = 5;
		assertEquals(bsTester.testBoundingSphere(ray), false);
		
		//entry only
		ray.startPos = new double[]{ 0, 0.5, 0 };
		ray.dir = new double[]{ 1, 0, 0 };
		ray.length = 1.7;
		assertEquals(bsTester.testBoundingSphere(ray), true);
		
		//exit only
		ray.startPos = new double[]{ 2.3, 0.5, 0 };
		ray.dir = new double[]{ 1, 0, 0 };
		ray.length = 1;
		assertEquals(bsTester.testBoundingSphere(ray), true);
		
		//exit only, infinite length
		ray.startPos = new double[]{ 2.3, 0.5, 0 };
		ray.dir = new double[]{ 1, 0, 0 };
		ray.length = Double.POSITIVE_INFINITY;
		assertEquals(bsTester.testBoundingSphere(ray), true);
		
		//entirely inside
		ray.startPos = new double[]{ 1.5, -0.1, 0.1 };
		ray.dir = new double[]{ 0, 1, 0 };
		ray.length = 0.1;
		assertEquals(bsTester.testBoundingSphere(ray), true);
		
		//sphere off end of ray
		ray.startPos = new double[]{ 2, 4, 0 };
		ray.dir = new double[]{ 0, -1, 0 };
		ray.length = 0.1;
		assertEquals(bsTester.testBoundingSphere(ray), false);
		
		//sphere off start of ray
		ray.startPos = new double[]{ 2, 0, 4 };
		ray.dir = new double[]{ 0, 0, 1 };
		ray.length = 0.1;
		assertEquals(bsTester.testBoundingSphere(ray), false);
		
		//just in
		ray.startPos = new double[]{ 0, 0.5, 0 };
		ray.dir = new double[]{ 1, 0, 0 };
		ray.length = 2 - Math.sqrt(3.0/4.0) + 1e-6;
		assertEquals(bsTester.testBoundingSphere(ray), true);
		

		//almost in
		ray.startPos = new double[]{ 0, 0.5, 0 };
		ray.dir = new double[]{ 1, 0, 0 };
		ray.length = 2 - Math.sqrt(3.0/4.0) - 1e-6;
		assertEquals(bsTester.testBoundingSphere(ray), false);
		
	}
	
	
	public void testBoundingSphereSpeed(){
		
		Square sq1 = new Square("sq1", new double[]{1,0,0}, new double[]{-1,0,0}, new double[]{0,0,1}, 0.1, 0.1, Absorber.ideal());
		//Disc sq2 = new Disc("sq2", new double[]{1.1,0,0}, new double[]{-1,0,0}, 0.2, Absorber.ideal());
		Dish sq2 = new Dish("sq2", new double[]{1.1,0,0}, new double[]{-1,0,0}, 100, 0.2, Absorber.ideal());
		
		Optic o = new Optic("o", new Element[]{ sq1, sq2 });

		for(int j=0; j < 100; j++){
			o.enableBoundingCheck = false;
			
			long t0 = System.currentTimeMillis();
			for(int i=0; i < 100000; i++){
				RaySegment ray = new RaySegment();
	
				ray.startPos = new double[]{0,0,0};
				ray.dir = new double[]{1,0,0};
				ray.length = Double.POSITIVE_INFINITY;
				ray.up = Util.createPerp(ray.dir);
				ray.E0 = new double[][]{ { 1,0,0,0 } };
				
				Tracer.trace(o, ray, 2, 0, true);
				Pol.recoverAll();
			}
			System.out.print(System.currentTimeMillis() - t0 + "\t");
			
			o.enableBoundingCheck = true;
			
			t0 = System.currentTimeMillis();
			for(int i=0; i < 100000; i++){
				RaySegment ray = new RaySegment();
	
				ray.startPos = new double[]{0,0,0};
				ray.dir = new double[]{1,0,0};
				ray.length = Double.POSITIVE_INFINITY;
				ray.up = Util.createPerp(ray.dir);
				ray.E0 = new double[][]{ { 1,0,0,0 } };
				
				Tracer.trace(o, ray, 2, 0, true);
				Pol.recoverAll();
			}
			System.out.println(System.currentTimeMillis() - t0);
		}
	}
}
