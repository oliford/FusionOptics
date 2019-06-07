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

import fusionOptics.types.Interface;
import fusionOptics.types.Medium;
import algorithmrepository.exceptions.NotImplementedException;
import net.jafama.FastMath;

/**  Interface with single layer anti-reflective coating.
 * 
 * Currently ignores overall phase shifts
 *  
 * @author oliford
 *
 */
public class IsoIsoAntiReflective extends IsoIsoInterface {
	
	
	private IsoIsoInterface iface1, iface2;
	private double nC;
	private double thickness;
	
	//not currently valid, because we do everything real
	private IsoIsoAntiReflective(IsoIsoInterface iface1, IsoIsoInterface iface2, double layerIndex, double layerThickness, double absorbedIntensityFraction) {
		super(absorbedIntensityFraction);
		this.nC = layerIndex;
		this.thickness = layerThickness;
		this.iface1 = iface1;
		this.iface2 = iface2;
	}
	
	/** Cannot be implemented. There is no ideal since we need to know the intended wavelength. This routine is only
	 * here to stop users accidently picking up IsoIsoInterface.ideal() and not getting a coating at all */
	public static IsoIsoAntiReflective ideal() {
		throw new NotImplementedException();
	}
	
	/**
	 * 
	 * @param layerIndex	Typical 
	 * @param layerThickness
	 * @param absorbedIntensityFraction
	 */
	public IsoIsoAntiReflective(double layerIndex, double layerThickness, double absorbedIntensityFraction) {
		super(absorbedIntensityFraction);
		this.nC = layerIndex;
		this.thickness = layerThickness;
		this.iface1 = IsoIsoStdFresnel.ideal();
		this.iface2 = IsoIsoStdFresnel.ideal();
	}
	
	@Override
	public InterfaceCoeffs getTransmissionReflectionCoefficients(double nI, double nT, double cosThetaI, double cosThetaT, double wavelength) {
		
		double cosThetaC = FastMath.sqrt(1.0 - (1.0 - cosThetaI)*(nI*nI/nC*nC));
		double cosThetaC2 = FastMath.sqrt(1.0 - (1.0 - cosThetaT)*(nT*nT/nC*nC));
		
		InterfaceCoeffs ic12 = iface1.getTransmissionReflectionCoefficients(nI, nC, cosThetaI, cosThetaC, wavelength);
		InterfaceCoeffs ic23 = iface2.getTransmissionReflectionCoefficients(nC, nT, cosThetaC, cosThetaT, wavelength);
		//returns  [Rs, Rp, Ts, Tp][Re/Im]
		
		
		/* See the SVG file for the full derivation.         
		 */
	
		double beta = 2 * Math.PI * nC * thickness * cosThetaC / wavelength;
		double cos2Beta = FastMath.cos(2 * beta); 
		
		//Formulae from "Optical Coatings" document by CVI Melles Griot 
		// [ https://cvimellesgriot.com/products/Documents/TechnicalGuide/Optical-Coatings.pdf ]
		
		double Rs2 = (ic12.RsRe*ic12.RsRe + ic23.RsRe*ic23.RsRe + 2*ic12.RsRe*ic23.RsRe*cos2Beta) /
						(1 + ic12.RsRe*ic12.RsRe*ic23.RsRe*ic23.RsRe + 2*ic12.RsRe*ic23.RsRe*cos2Beta);
		
		double Rp2 = (ic12.RpRe*ic12.RpRe + ic23.RpRe*ic23.RpRe + 2*ic12.RpRe*ic23.RpRe*cos2Beta) /
						(1 + ic12.RpRe*ic12.RpRe*ic23.RpRe*ic23.RpRe + 2*ic12.RpRe*ic23.RpRe*cos2Beta);
		
		InterfaceCoeffs ret = new InterfaceCoeffs();
		
		//unfortunately they are intensity. For now, we guess that the phases
		//are just +/- 180' and the same as the ideal Fresnel gave:
		ret.RsRe = (((ic12.RsRe + ic23.RsRe)/2 < 0) ? -1 : 1) * FastMath.sqrt(Rs2);
		ret.RpRe = (((ic12.RpRe + ic23.RpRe)/2 < 0) ? -1 : 1) * FastMath.sqrt(Rp2);
		ret.TsRe = FastMath.sqrt(1.0 - Rs2);
		ret.TpRe = FastMath.sqrt(1.0 - Rp2);
		
		//ignore overall phases
		ret.RsIm = 0; ret.RpIm = 0; ret.TsIm = 0; ret.TpIm = 0; 
		 
		
		//return new double[]{ signS*FastMath.sqrt(Rs2), signP*FastMath.sqrt(Rp2), Ts, Tp };
		return ret;
	}

	@Override
	public int hashCode() {
		long t; int r = super.hashCode();
		r = 31 * r + iface1.hashCode();
		r = 31 * r + iface2.hashCode();		
		t = Double.doubleToLongBits(nC); 		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(thickness);	r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj) || getClass() != obj.getClass())
			return false;
		IsoIsoAntiReflective other = (IsoIsoAntiReflective) obj;
		return super.equals(obj) 
				&& iface1.equals(other.iface1)
				&& iface2.equals(other.iface2)
				&& Double.doubleToLongBits(nC) == Double.doubleToLongBits(other.nC)
				&& Double.doubleToLongBits(thickness) == Double.doubleToLongBits(other.thickness);
	}
	
}
