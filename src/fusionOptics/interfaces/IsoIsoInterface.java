/**
 * Copyright 2011 Oliver Ford
 *
 * This file is part of the minerva-optics 'RayTracer'.
 *
 *   RayTracer is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   RayTracer is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RayTracer.  If not, see <http://www.gnu.org/licenses/>.
 *   
 *   @author oliford <codes<at>oliford.co.uk>
 */
package fusionOptics.interfaces;

import jafama.FastMath;



import java.util.Arrays;
import java.util.List;

import fusionOptics.Util;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;



/** Snell's law based interface between two isotropic media. 
 * 
 * The base class is 100% Transmission only, no reflection, but this can be overridden
 * by derived classes which must provide complex Reflection and Transmission coefficients.
 * 
 */
public class IsoIsoInterface extends DualMediumInterface {
	private double nonAbsorbedAmplitudeCoeff;
	
	private final static IsoIsoInterface ideal = new IsoIsoInterface(0.0);
	
	public static IsoIsoInterface ideal(){ return ideal; }
	
	
	public IsoIsoInterface(double absorbedIntensityFraction) {
		this.nonAbsorbedAmplitudeCoeff = FastMath.sqrt(1.0 - absorbedIntensityFraction);
	}

	@Override
	public void checkCompatibility(Surface surface) {
		if(!ignoreInterfaceCompatibility && (
				(surface.getFrontMedium() != null && surface.getFrontMedium().getMaterial().getNAxes() > 0) || 
				(surface.getBackMedium() != null && surface.getBackMedium().getMaterial().getNAxes() > 0)) )
			throw new IllegalArgumentException("IsoIsoInterface and it's derivatives can only handle isotropic media");
	}

	protected class InterfaceCoeffs {
		double RsRe, RsIm;
		double RpRe, RpIm;
		double TsRe, TsIm;
		double TpRe, TpIm;
		
		public double magRs2(){ return RsRe*RsRe + RsIm*RsIm; }
		public double magRp2(){ return RpRe*RpRe + RpIm*RpIm; }
		public double magTs2(){ return TsRe*TsRe + TsIm*TsIm; }
		public double magTp2(){ return TpRe*TpRe + TpIm*TpIm; }
		
	}
	
	
	/** @return pure unaltered transmission
	 */
	public InterfaceCoeffs getTransmissionReflectionCoefficients(double nI, double nT, double cosThetaI, double cosThetaT, double wavelength){
		InterfaceCoeffs ret = new InterfaceCoeffs();
		ret.RsRe = 0; ret.RsIm = 0;
		ret.RpRe = 0; ret.RpIm = 0;
		ret.TsRe = 1; ret.TsIm = 0;
		ret.TpRe = 1; ret.TpIm = 0;
		return ret;
	}
	

	@Override
	public final void calcIntersection(Intersection hit, Medium incidentMedium,
			Medium transmissionMedium, double minIntensity) {

		//we need to use the object we hit's refractive index ratio
		double nI = (incidentMedium == null ? 1 : incidentMedium.getRefractiveIndex(0, hit.incidentRay.wavelength));
		double nT = (transmissionMedium == null ? 1 : transmissionMedium.getRefractiveIndex(0, hit.incidentRay.wavelength));
		double rIndexRatio = nI / nT;

		
		double cosThetaI = -Util.dot(hit.incidentRay.dir, hit.normal);
		
		if( (1 - cosThetaI*cosThetaI) >= 1/(rIndexRatio*rIndexRatio)){
			//total internal refraction
			Reflector.pureReflection(hit, minIntensity, nonAbsorbedAmplitudeCoeff, true);
			
			return;
		}

		//otherwise we're creating a transmitted ordinary ray
		
		double a = + rIndexRatio * cosThetaI;
		double b = - Math.sqrt( 1 - rIndexRatio*rIndexRatio*(1-cosThetaI*cosThetaI));
		double c = a + b;
		
		double transmissionDir[] = Util.reNorm(new double[]{
				c * hit.normal[0] + rIndexRatio * hit.incidentRay.dir[0], 
				c * hit.normal[1] + rIndexRatio * hit.incidentRay.dir[1], 
				c * hit.normal[2] + rIndexRatio * hit.incidentRay.dir[2], 
			});
		
		//calculate the direction of the 's' polaraistion (the component in the plane)
		double sPolDir[] = Util.cross(hit.incidentRay.dir, hit.normal);
		if(Util.dot(sPolDir, sPolDir) == 0) 
			sPolDir = hit.incidentRay.up; //s and p are the same, keep the original
		else
			sPolDir = Util.reNorm(sPolDir); //otherwise, normalise it, so it's just a direction
		
		//rotate incident ray's polarisation definition to the s/p frame
		hit.incidentRay.rotatePolRefFrame(sPolDir);		
		
		double cosThetaT = -Util.dot(transmissionDir, hit.normal);
		InterfaceCoeffs ic = getTransmissionReflectionCoefficients(nI, nT, cosThetaI, cosThetaT, hit.incidentRay.wavelength);
		
		if((ic.magTs2() + ic.magTp2()) > 0.0) {
			hit.transmittedOrdinary = new RaySegment();
			hit.transmittedOrdinary.startHit = hit;
			hit.transmittedOrdinary.startPos = hit.pos;
			
			hit.transmittedOrdinary.dir = transmissionDir;
			
			//the 's' of this intersection is new sense of 'up' for both rays. (since we have rotated to this frame anyway)
			hit.transmittedOrdinary.up = sPolDir;

			//Transfer the appropriate amount of the polarised electric field amplitudes to the new rays
			//The phases will be taken care of by the +/- in the coefficients.
			
			hit.transmittedOrdinary.E0 = Pol.complexMulAll(hit.incidentRay.E1, 
					ic.TsRe, ic.TsIm, ic.TpRe, ic.TpIm, nonAbsorbedAmplitudeCoeff);
			
			if(hit.transmittedOrdinary.startIntensity() < minIntensity){
				hit.transmittedOrdinary = null;
			} else {
				//set-up both rays for tracing on
				hit.transmittedOrdinary.length = Double.POSITIVE_INFINITY;
				hit.transmittedOrdinary.medium = transmissionMedium;
				hit.transmittedOrdinary.wavelength = hit.incidentRay.wavelength;
				hit.transmittedOrdinary.endHit = null;
			}
		}
		
		if((ic.magRs2() + ic.magRp2()) > 0.0) {
				
			hit.reflectedOrdinary = new RaySegment();
			hit.reflectedOrdinary.startHit = hit;
			hit.reflectedOrdinary.startPos = hit.pos;
			
			hit.reflectedOrdinary.dir = Util.reNorm(new double[]{
					hit.incidentRay.dir[0] + (2 * hit.normal[0] * cosThetaI),
					hit.incidentRay.dir[1] + (2 * hit.normal[1] * cosThetaI),
					hit.incidentRay.dir[2] + (2 * hit.normal[2] * cosThetaI),
				});
			
			hit.reflectedOrdinary.up = sPolDir.clone(); //don't link arrays, just in case
			
			hit.reflectedOrdinary.E0 = Pol.complexMulAll(hit.incidentRay.E1,
					ic.RsRe, ic.RsIm, ic.RpRe, ic.RpIm, nonAbsorbedAmplitudeCoeff);	
			
			if(hit.reflectedOrdinary.startIntensity() < minIntensity){
				hit.reflectedOrdinary = null;
			} else {
				hit.reflectedOrdinary.length = Double.POSITIVE_INFINITY;
				hit.reflectedOrdinary.medium = transmissionMedium;
				hit.reflectedOrdinary.wavelength = hit.incidentRay.wavelength;
				hit.reflectedOrdinary.endHit = null;
			}
		}
	}

	@Override
	public int hashCode() {
		long t; int r = 1;
		t = Double.doubleToLongBits(nonAbsorbedAmplitudeCoeff); r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		IsoIsoInterface other = (IsoIsoInterface) obj;
		return Double.doubleToLongBits(nonAbsorbedAmplitudeCoeff) == Double.doubleToLongBits(other.nonAbsorbedAmplitudeCoeff);
			
	}
}
