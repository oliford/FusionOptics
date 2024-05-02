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
package fusionOptics.materials;

import fusionOptics.types.Material;

/** Ideal glass with refractive index that has a linear dependence on wavelength.
 *
 * @author oliford
 */
public class IsotropicLinearDispersiveGlass extends Material {
	public final static double lambda_d = 587.56e-9;
	public final static double lambda_F = 486.13e-9;
	public final static double lambda_C = 656.27e-9;
	
	private double refractiveIndex_d;
	private double abbeNumber_d;
	private double transmission; 
	private double dNdLambda;
	
	/** Create ideal glass with given refractive index and perfect transmission (no absorption).
	 * 
	 * @param refractiveIndexD	Refractive index at sodium d line (587.56nm)
	 * @param abbeNumberD	Abbe Number of Sodium d line = (nd - 1) / (nf - nc), nd = n(587.56nm), nf=n(486.13nm), nc=n(656.27nm)
	 * 
	 */
	public IsotropicLinearDispersiveGlass(double refractiveIndex_d, double abbeNumber_d) {
		set(refractiveIndex_d, abbeNumber_d);
	}
	
	/** Create ideal glass with given refractive index and given transmission. */
	public IsotropicLinearDispersiveGlass(double refractiveIndex_d, double abbeNumber_d, double transmission) {
		this.transmission = transmission;
		set(refractiveIndex_d, abbeNumber_d);
	}
	
	@Override
	public int getNAxes() { return 0; }

	@Override
	public double getRefractiveIndex(int modeNumber, double wavelength, double temperature) { 
		return refractiveIndex_d + dNdLambda * (wavelength - lambda_d);
	}

	@Override
	public double getTransmission(int modeNumber, double wavelength, double temperature) { return transmission; }

	@Override
	public double getVerdetConstant(int modeNumber, double wavelen, double temperature) { 
		throw new UnsupportedOperationException();
		}

	public void setRefractiveIndexD(double nd) { set(nd, abbeNumber_d); }

	public void set(double refractiveIndex_d, double abbeNumber_d){
		this.refractiveIndex_d = refractiveIndex_d;
		this.abbeNumber_d = abbeNumber_d;
		this.transmission = 1.0;
		this.dNdLambda = (abbeNumber_d == 0) ? 0 : (refractiveIndex_d - 1.0) / (abbeNumber_d * (lambda_F - lambda_C));
	}

	public double getRefractiveIndexD() { return refractiveIndex_d; }

	@Override
	public int hashCode() {
		int r = 1;
		long t;		
		t = Double.doubleToLongBits(abbeNumber_d); 		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(dNdLambda); 		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(refractiveIndex_d);	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(transmission); 		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(lambda_C); 			r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(lambda_d); 			r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(lambda_F);			r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		IsotropicLinearDispersiveGlass other = (IsotropicLinearDispersiveGlass) obj;
		return Double.doubleToLongBits(abbeNumber_d) == Double.doubleToLongBits(other.abbeNumber_d)
					&& Double.doubleToLongBits(dNdLambda) == Double.doubleToLongBits(other.dNdLambda)
					&& Double.doubleToLongBits(refractiveIndex_d) == Double.doubleToLongBits(other.refractiveIndex_d)
					&& Double.doubleToLongBits(transmission) == Double.doubleToLongBits(other.transmission);
	}
	
	
}
