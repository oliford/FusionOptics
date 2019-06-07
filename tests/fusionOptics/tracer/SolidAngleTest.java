package fusionOptics.tracer;


import net.jafama.FastMath;

import java.util.Vector;

import junit.framework.TestCase;
import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.collection.IntensityInfo;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.Reflector;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Optic;
import fusionOptics.types.RaySegment;
import oneLiners.OneLiners;
import svg.SVGSplitView3D;

/** Test the target and effective solid angle calculations */
public class SolidAngleTest extends TestCase {
	
	final static double reHitTolerence = 1e-6;
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing";
	
	/** fire multiple rays at a square  */
	public void testSolidAngle(){
		
		double diameter = 0.040;
		double distance = 1.900;
		double allowedFractionalError = 0.03; //pretty poor, but I think it has to do with the theoretical assumptions
		
		int nRays = 100000;
				
		//Disc target = new Disc("target", new double[]{ distance, 0.0, 0.0}, new double[]{ -1.0, 0.0, 0.0}, diameter/2, Absorber.ideal());
		
		Square target = new Square("target", new double[]{ distance, 0.0, 0.0}, new double[]{ -1.0, 0.0, 0.0}, new double[]{ 0.0, 0.0, 1.0}, 1.2*diameter, 1.2*diameter, Absorber.ideal());		
		Iris aperture = new Iris("aperture", new double[]{ distance, 0, 0}, new double[]{ -1, 0, 0}, 0.8 * diameter, diameter / 2, Absorber.ideal());
		Optic sys = new Optic("sys", new Element[]{ target, aperture });
		
		IntensityInfo iInfo = new IntensityInfo(target);
		double targetSolidAngle = Double.NaN;
		
		for(int i=0;i<nRays;i++){
			RaySegment ray = new RaySegment();
			
			ray.startPos = new double[]{ 0,0,0 };		
			
			ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, sys, true);			
			targetSolidAngle = Util.length(ray.dir); 
			ray.dir = Util.reNorm(ray.dir);
			
			ray.length = Double.POSITIVE_INFINITY;
			ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, new double[]{0,0,1})), ray.dir);
			ray.E0 = new double[][]{ {1,0,0,0} };
			
			Tracer.trace(sys, ray, 2, 1e-10, false);
			
			ray.processIntersections(target, iInfo);
		}
		
		double theoretical = Math.PI * diameter*diameter / 4 / distance / distance;
		double effectiveSolidAngle = iInfo.getSourceSolidAng(target, targetSolidAngle, nRays);
		
		iInfo.dump(true);		
		
		System.out.println("target solid angle = " + (targetSolidAngle * 1e6) + " µSR");		
		System.out.println("effective solid angle = " + effectiveSolidAngle * 1e6 + " µSR");
		System.out.println("theoretical solid angle = " + theoretical * 1e6 + " µSR");
		System.out.println("error = " + FastMath.abs((theoretical - effectiveSolidAngle)/theoretical * 100) + " %");
		
		assertEquals(theoretical, effectiveSolidAngle, allowedFractionalError * theoretical);
		
	} 
	
}
