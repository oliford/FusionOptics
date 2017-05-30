package fusionOptics.interfaces;

import jafama.FastMath;

import java.util.List;

import fusionOptics.Util;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;


/** Very simple ideal reflector, with no media.
 *  Can also absorb some at the interface.
 */ 
public class Reflector implements Interface {
	public boolean ignoreInterfaceCompatibility = false;
	
	/** Keep a single global ideal reflector. */
	private static Reflector ideal = new Reflector(0.0);

	private double nonAbsorbedAmplitudeCoeff;
	
	/** Inverts the amplitude (phase + 180°) of the p component upon reflection.
	 * This is normal for metalic mirrors and other simple reflections at
	 * near normal incidence. It means that the E field direction of the incoming
	 * and outgoing rays are parallel at the surface.
	 * 
	 * This generally isn't true at glancing incidence, or (apparently, but I'm not sure) 
	 * for dielectric mirrors.
	 * */
	private boolean invertPAmplitude = true;
	
	public Reflector(double absorbedIntensityFraction) {
		this.nonAbsorbedAmplitudeCoeff = FastMath.sqrt(1.0 - absorbedIntensityFraction);
		
	}
	public Reflector(double absorbedIntensityFraction, boolean invertPAmplitude){		
		this.nonAbsorbedAmplitudeCoeff = FastMath.sqrt(1.0 - absorbedIntensityFraction);
		this.invertPAmplitude = invertPAmplitude;
	}
	
	/** Returns a perfect reflection interface (global single instance) */
	public static Reflector ideal(){ return ideal; }
	
	
	@Override
	public void calcIntersection(Intersection hit, double minIntensity) {
		pureReflection(hit, minIntensity, nonAbsorbedAmplitudeCoeff, invertPAmplitude);
	}
	
	/** Exposed static handler, also used elsewhere for e.g. total internal reflection */ 
	public static void pureReflection(Intersection hit, double minIntensity, double nonAbsorbedAmplitudeCoeff, boolean invertPAmplitude) {
			hit.reflectedOrdinary = new RaySegment();
			hit.reflectedOrdinary.startHit = hit;
			hit.reflectedOrdinary.startPos = hit.pos;
			
			double c = Util.dot(hit.normal, hit.incidentRay.dir);
			
			hit.reflectedOrdinary.dir = Util.reNorm(new double[]{
					hit.incidentRay.dir[0] - (2 * hit.normal[0] * c),
					hit.incidentRay.dir[1] - (2 * hit.normal[1] * c),
					hit.incidentRay.dir[2] - (2 * hit.normal[2] * c),
				});
			
			//calculate the direction of the 's' polaraistion (the component in the plane)
			double sPolDir[] = Util.cross(hit.incidentRay.dir, hit.normal);
			if(Util.dot(sPolDir, sPolDir) == 0) 
				sPolDir = hit.incidentRay.up.clone(); //s and p are the same, keep the original
			else
				sPolDir = Util.reNorm(sPolDir); //otherwise, normalise it, so it's just a direction
			
			//rotate incident ray's polarisation definition to the s/p frame
			hit.incidentRay.rotatePolRefFrame(sPolDir);
			
			hit.reflectedOrdinary.up = sPolDir;

			
			/*For normal incident radiation, the input polarisation has to be exactly parallel to the reflected polarisation, by symmetry.
			 * To ensure this, we need to flip the polarisation in the exiting ray description because that ray's sense of left/right is 
			 * the opposite to that of the incident ray. 
			 * 
			 * People call this a '1/2 phase shift' - which is wrong as far as I'm concerned. It is JUST a change of coordinate description
			 * of the same polarisation direction. 
			 */
			if(invertPAmplitude){
				hit.reflectedOrdinary.E0 = Pol.complexMulAll(hit.incidentRay.E1, 1, 1, -1, -1, nonAbsorbedAmplitudeCoeff);
			}else{
				hit.reflectedOrdinary.E0 = Pol.complexMulAll(hit.incidentRay.E1, 1, 1, 1,1, nonAbsorbedAmplitudeCoeff);
			}
			
			hit.reflectedOrdinary.length = Double.POSITIVE_INFINITY;
			
			hit.reflectedOrdinary.medium = hit.incidentRay.medium; //still travelling in the same medium (whatever it was)
			hit.reflectedOrdinary.wavelength = hit.incidentRay.wavelength;
			hit.reflectedOrdinary.endHit = null;
			
	}

	@Override
	public void checkCompatibility(Surface surface) {
		if(!ignoreInterfaceCompatibility && (
				(surface.getFrontMedium() != null && surface.getFrontMedium().getMaterial().getNAxes() > 1) || 
				(surface.getBackMedium() != null && surface.getBackMedium().getMaterial().getNAxes() > 1)) )
			throw new IllegalArgumentException("Reflector can only handle isotropic (or null) media");
	}
	
	@Override
	public int hashCode() {
		long t; int r = 1;
		r = 31 * r + (invertPAmplitude ? 1231 : 1237);		
		t = Double.doubleToLongBits(nonAbsorbedAmplitudeCoeff); r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		Reflector other = (Reflector) obj;
		return invertPAmplitude == other.invertPAmplitude
			&& Double.doubleToLongBits(nonAbsorbedAmplitudeCoeff) == Double.doubleToLongBits(other.nonAbsorbedAmplitudeCoeff);
			
	}	
	
}
