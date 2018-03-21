package fusionOptics.interfaces;

import jafama.FastMath;

import java.util.List;

import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.types.Element;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;


/** Very simple ideal diffraction grating, with no media.
 *  
 */ 
public class DiffractionGrating implements Interface {
	public boolean ignoreInterfaceCompatibility = false;
	
	private int order;	
	private double gratingDensity;
	private double nonAbsorbedAmplitudeCoeff;
	private double[] ruleDirection;

	/** Inverts the amplitude (phase + 180°) of the p component upon reflection.
	 * This is normal for metalic mirrors and other simple reflections at
	 * near normal incidence. It means that the E field direction of the incoming
	 * and outgoing rays are parallel at the surface.
	 * 
	 * This generally isn't true at glancing incidence, or (apparently, but I'm not sure) 
	 * for dielectric mirrors.
	 * */
	private boolean invertPAmplitude = true;
	
	public DiffractionGrating(int order, double gratingDensity, double[] ruleDirection) {
		this.order = order;
		this.gratingDensity = gratingDensity;
		this.ruleDirection = ruleDirection;
		this.nonAbsorbedAmplitudeCoeff = 1.0;		
	}
	
	public DiffractionGrating(int order, double gratingDensity, double[] ruleDirection, double absorbedIntensityFraction) {
		this.order = order;
		this.gratingDensity = gratingDensity;
		this.ruleDirection = ruleDirection;
		this.nonAbsorbedAmplitudeCoeff = FastMath.sqrt(1.0 - absorbedIntensityFraction);		
	}
		
	@Override
	public void calcIntersection(Intersection hit, double minIntensity) {
		
		hit.reflectedOrdinary = new RaySegment();
			hit.reflectedOrdinary.startHit = hit;
			hit.reflectedOrdinary.startPos = hit.pos;
			
			//define directions: r=rule dir, d=diffraction dir in plane, a=incident dir
			double a[] = Util.reNorm(hit.incidentRay.dir);
			double d[] = Util.reNorm(Util.cross(hit.normal, ruleDirection));
			double r[] = Util.reNorm(Util.cross(d, hit.normal));
			
			
			double sinThetaI = Util.dot(a, d);
			double aDotR = Util.dot(a, r); 
			double aDotN = Util.dot(a, hit.normal); 
			double sinThetaM = sinThetaI - order * hit.incidentRay.wavelength * gratingDensity;
			
			double bDotN = FastMath.sqrt(1 - sinThetaM*sinThetaM);
			if(aDotN > 0) bDotN *= -1; //get direction opposite to incoming ray
			
			hit.reflectedOrdinary.dir = Util.reNorm(new double[]{
					sinThetaM * d[0] + aDotR * r[0] + bDotN * hit.normal[0], 
					sinThetaM * d[1] + aDotR * r[1] + bDotN * hit.normal[1],
					sinThetaM * d[2] + aDotR * r[2] + bDotN * hit.normal[2]
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
		r = 31 * r + Integer.hashCode(order);
		t = Double.doubleToLongBits(nonAbsorbedAmplitudeCoeff); r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		DiffractionGrating other = (DiffractionGrating) obj;
		return Double.doubleToLongBits(nonAbsorbedAmplitudeCoeff) == Double.doubleToLongBits(other.nonAbsorbedAmplitudeCoeff) 
				&& order == other.order;
			
	}	
	
}
