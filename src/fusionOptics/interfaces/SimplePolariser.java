package fusionOptics.interfaces;

import net.jafama.FastMath;

import java.util.Arrays;
import java.util.List;

import fusionOptics.Util;
import fusionOptics.materials.Vacuum;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;


/** Very simplified polariser. Can represent either a reflecting
 * or absorbing polariser.
 * 
 * @author oliford */
public class SimplePolariser implements Interface {
	
	/** The amount of the perpendicular polarisation
	 * that is reflected rather than absorbed */
	private double ampltiudeReflectionCoeff = 0;
	
	/** Direction of passing polarisation */
	private double polDir[]; 
	
	public SimplePolariser(double polDir[], double reflectedIntensityFraction) {
		this.polDir = polDir;
		this.ampltiudeReflectionCoeff = FastMath.sqrt(ampltiudeReflectionCoeff);
	}

	@Override
	public void calcIntersection(Intersection hit, double minIntensity) {
		
		double c = Util.dot(hit.normal, hit.incidentRay.dir);
		
		double reflDir[] = Util.reNorm(new double[]{
				hit.incidentRay.dir[0] - (2 * hit.normal[0] * c),
				hit.incidentRay.dir[1] - (2 * hit.normal[1] * c),
				hit.incidentRay.dir[2] - (2 * hit.normal[2] * c),
			});
		
		
		//rotate incident ray's polarisation to the passing direction
		hit.incidentRay.rotatePolRefFrame(polDir);
		
		
		hit.transmittedOrdinary = new RaySegment();
		hit.transmittedOrdinary.startHit = hit;
		hit.transmittedOrdinary.startPos = hit.pos;
		hit.transmittedOrdinary.dir = hit.incidentRay.dir.clone();
		 
		
		hit.transmittedOrdinary.up = polDir;

		//now we are free to just copy the polarisation info
		
		hit.transmittedOrdinary.E0 = Pol.complexMulAll(hit.incidentRay.E1, 1, 1, 0, 0, 1);
			
		if(hit.transmittedOrdinary.startIntensity() < minIntensity)
			hit.transmittedOrdinary = null;
		else{
			hit.transmittedOrdinary.length = Double.POSITIVE_INFINITY;
			
			hit.transmittedOrdinary.medium = hit.incidentRay.medium; //still travelling in the same medium (whatever it was)
			hit.transmittedOrdinary.wavelength = hit.incidentRay.wavelength;
			hit.transmittedOrdinary.endHit = null;
		}
		
		if(FastMath.abs(ampltiudeReflectionCoeff) != 0) {
			hit.reflectedOrdinary = new RaySegment();
			hit.reflectedOrdinary.startHit = hit;
			hit.reflectedOrdinary.startPos = hit.pos;
			hit.reflectedOrdinary.dir = reflDir;
			
			hit.reflectedOrdinary.up = polDir.clone();
			
			/* We need to flip the polarisation in the exiting ray description because that ray's sense of left/right is 
			 * the opposite to that of the incident ray. See Reflector.pureReflection() */ 
			hit.reflectedOrdinary.E0 = Pol.complexMulAll(hit.incidentRay.E1, 0, 0, -1, -1, ampltiudeReflectionCoeff);
			
			if(hit.reflectedOrdinary.startIntensity() < minIntensity)
				hit.reflectedOrdinary = null;
			else{
				hit.reflectedOrdinary.length = Double.POSITIVE_INFINITY;
					
				hit.reflectedOrdinary.medium = hit.incidentRay.medium; //still travelling in the same medium (whatever it was)
				hit.reflectedOrdinary.wavelength = hit.incidentRay.wavelength;
				hit.reflectedOrdinary.endHit = null;
			}
		}
	}

	@Override
	public void checkCompatibility(Surface surface) {
		if( (surface.getFrontMedium() != null && 
				surface.getFrontMedium().getMaterial() != null && 
				!(surface.getFrontMedium().getMaterial() instanceof Vacuum)) ||
			(surface.getBackMedium() != null && 
					surface.getFrontMedium().getMaterial() != null && 
					!(surface.getBackMedium().getMaterial() instanceof Vacuum)) )
			throw new IllegalArgumentException("SimplePolariser interface can only be used on Vacuum-Vacuum surfaces.");
			
	}

	@Override
	public int hashCode() {
		long t; int r = 1;
		r = 31 * r + Arrays.hashCode(polDir);		
		t = Double.doubleToLongBits(ampltiudeReflectionCoeff); r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		SimplePolariser other = (SimplePolariser) obj;
		return Arrays.equals(polDir, other.polDir)
			&& Double.doubleToLongBits(ampltiudeReflectionCoeff) == Double.doubleToLongBits(other.ampltiudeReflectionCoeff);
			
	}
}
