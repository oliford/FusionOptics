package fusionOptics.tracer;


import java.util.Vector;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.interfaces.Reflector;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.RaySegment;

import oneLiners.OneLiners;

import svg.SVGSplitView3D;

/** First rays at a screen - test of correct solid angle distribution */
public class FireRaysAtSquareTest {
	final static double reHitTolerence = 1e-6;
	//final static String outPath = System.getProperty("java.io.tmpdir") + "/rayTracing";
	//final static String outPath = "e:/temp/rayTracingOutput";
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing";
	
	/** fire multiple rays at a square  */
	public static void main(String[] args) {
		
		int nRays = 10000000;
		
		//SVGSplitView3D svg = new SVGSplitView3D(outPath + "/multipleRaysAtASquare", new double[]{ 0, -0.5, -0.5,  1, 0.5, 0.5 });
		
		//svg.addLineStyle("optic", "none", 0.01, "green");
		
		Square sq = new Square("sq", new double[]{ 1.0, 0.0, 0.0}, new double[]{ -1.0, 0.0, 0.0}, new double[]{ 0.0, 0.0, 1.0}, 2, 4, Reflector.ideal());
		double hits[][] = new double[nRays][3];
		long t0 = System.currentTimeMillis();
		for(int i=0;i<nRays;i++){
			RaySegment ray = new RaySegment();
			
			ray.startPos = new double[]{ 0,0,0 };		
			ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, sq);
			ray.length = Double.POSITIVE_INFINITY;
			ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, new double[]{0,0,1})), ray.dir);
			ray.E0 = new double[][]{ {1,0,0,0} };
			
			Tracer.trace(sq, ray, 2, 1e-10, false);
			
			if(ray.endHit != null) {
				//svg.addSphere(ray.endHit.pos, 0.01);
				hits[i] = ray.endHit.pos;
			}
		}

		System.out.println((System.currentTimeMillis() - t0) + " ms");
		//svg.addLines(sq.draw(), "optic");
		
		//svg.destroy();
		
		binaryMatrixFile.BinaryMatrixFile.mustWrite(outPath + "/hits.txt", hits, false);
	} //*/
	
}
