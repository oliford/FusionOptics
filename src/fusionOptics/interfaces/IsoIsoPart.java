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

import fusionOptics.types.Medium;
import net.jafama.FastMath;

/** Simple uncoated standard interface between two isotropic media
 * Reflection and transmission coefficients given by user
 * 
 * Simplified reflection, ignoring phase and polarisation
 *  
 * @author oliford
 *
 */
public class IsoIsoPart extends IsoIsoInterface {
	
	private final static IsoIsoPart ideal = new IsoIsoPart(1.0, 0.0);
	protected double reflectedIntensityFraction;

	public final static IsoIsoPart ideal(){ return ideal; }
	
	public IsoIsoPart(double reflectedIntensityFraction, double absorbedIntensityFraction) {
		super(absorbedIntensityFraction);
		this.reflectedIntensityFraction = reflectedIntensityFraction;
	}
	
	@Override
	public InterfaceCoeffs getTransmissionReflectionCoefficients(double nI, double nT, double cosThetaI, double cosThetaT, double wavelength) {
		InterfaceCoeffs ret = new InterfaceCoeffs();
		//calculate the field amplitude fraction of s and p that are transmitted and reflected (Fresnel ampltiude equations)
		// [ http://physics.tamuk.edu/~suson/html/4323/prop-em.html ] has a complete derivation and very careful discussion
		//ret.RsRe = (nI*cosThetaI - nT*cosThetaT) / (nI*cosThetaI + nT*cosThetaT);
		//ret.RpRe = (nT*cosThetaI - nI*cosThetaT) / (nI*cosThetaT + nT*cosThetaI);
		//ret.RsIm = 0; ret.RpIm = 0;
		
		//ret.TsRe = 2*nI*cosThetaI / (nI*cosThetaI + nT*cosThetaT);
		//ret.TpRe = 2*nI*cosThetaI / (nI*cosThetaT + nT*cosThetaI);
		//ret.TsIm = 0; ret.TpIm = 0;
		ret.RsRe = reflectedIntensityFraction; ret.RsIm = 0;
		ret.RpRe = reflectedIntensityFraction; ret.RpIm = 0;
		ret.TsRe = 1-reflectedIntensityFraction; ret.TsIm = 0;
		ret.TpRe = 1-reflectedIntensityFraction; ret.TpIm = 0;
		//E fields in the transmitted beam are weaker because it would be wider than the incident/reflected
		//Since we want to model them as integrated E over area, we need to correct for that change in cross-sectional area
		//Also, the power passing the interface in 1 unit time is held in a smaller volume in the medium
		//in which the light is travelling slower, so it has a higher energy density and hence a higher E field.
		double deltaCSA = FastMath.sqrt((nT * cosThetaT) / (nI * cosThetaI));
		ret.TsRe *= deltaCSA;
		ret.TpRe *= deltaCSA;
		
		return ret;
	}
	
}
