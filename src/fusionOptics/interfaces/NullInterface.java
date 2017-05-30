package fusionOptics.interfaces;

import jafama.FastMath;


import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.List;

import fusionOptics.Util;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;


/** Does nothing but entirely transmit the ray.
 * Useful for target, testing and debugging surfaces.
 */ 
public class NullInterface extends DualMediumInterface {
	/** Keep a single global ideal nothing. */
	private static NullInterface ideal = new NullInterface();
	
	/** Don't instantiate, use .ideal() */
	private NullInterface() { }
	
	/** Returns a perfect null interface (global single instance) */
	public static NullInterface ideal(){ return ideal; }

	@Override
	public void calcIntersection(Intersection hit, Medium incidentMedium, Medium transmissionMedium, double minIntensity){

			hit.transmittedOrdinary = new RaySegment();
			hit.transmittedOrdinary.startHit = hit;
			hit.transmittedOrdinary.startPos = hit.pos;
					
			hit.transmittedOrdinary.dir = hit.incidentRay.dir.clone();
			
			hit.transmittedOrdinary.up = hit.incidentRay.up.clone();

			//now we are free to just copy the polarisation info
			hit.transmittedOrdinary.E0 = Pol.copyAll(hit.incidentRay.E1);	
			
			hit.transmittedOrdinary.length = Double.POSITIVE_INFINITY;
			
			hit.transmittedOrdinary.raySpecificRefractiveIndex = hit.incidentRay.raySpecificRefractiveIndex;
			
			hit.transmittedOrdinary.medium = transmissionMedium;
			hit.transmittedOrdinary.wavelength = hit.incidentRay.wavelength;
			hit.transmittedOrdinary.endHit = null;
			
	}

	@Override
	public void checkCompatibility(Surface surface) {
		//everything is ok
	}
	
	//all nulls are alike
	@Override
	public int hashCode() { return 0; } //must be null-like hashcode
	@Override
	public boolean equals(Object obj) { return (obj == null) || obj instanceof Absorber; }	
}
