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

/** Ideal glass with fixed refractive index and transmission 
 *  coefficient for all wavelengths and temperatures.
 *
 * @author oliford
 */
public class UniaxialFixedIndexGlass extends Material {
	private double refractiveIndexOrdinary, refractiveIndexExtraordinary;
	private double transmission; 
	
	/** Create ideal glass with given refractive index and perfect transmission (no absorption) */
	public UniaxialFixedIndexGlass(double refractiveIndexOrdinary, double refractiveIndexExtraordinary) {
		this.refractiveIndexOrdinary = refractiveIndexOrdinary;
		this.refractiveIndexExtraordinary = refractiveIndexExtraordinary;
		this.transmission = 1.0;
	}
	
	/** Create ideal glass with given refractive index and given transmission. */
	public UniaxialFixedIndexGlass(double refractiveIndexOrdinary, double refractiveIndexExtraordinary, double transmission) {
		this.refractiveIndexOrdinary = refractiveIndexOrdinary;
		this.refractiveIndexExtraordinary = refractiveIndexExtraordinary;
		this.transmission = transmission;
	}
	
	@Override
	public int getNAxes() { return 1; }

	@Override
	public double getRefractiveIndex(int modeNumber, double wavelength, double temperature) { 
		return modeNumber == 0 ? refractiveIndexOrdinary : refractiveIndexExtraordinary;
	}

	@Override
	public double getTransmission(int modeNumber, double wavelength, double temperature) { return transmission; }

	@Override
	public double getVerdetConstant(int modeNumber, double wavelen, double temperature) { 
		throw new UnsupportedOperationException();
		}

	@Override
	public int hashCode() {
		long t; int r = 1;		
		t = Double.doubleToLongBits(refractiveIndexExtraordinary);	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(refractiveIndexOrdinary);		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(transmission);					r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		UniaxialFixedIndexGlass other = (UniaxialFixedIndexGlass) obj;
		return Double.doubleToLongBits(refractiveIndexExtraordinary) == Double.doubleToLongBits(other.refractiveIndexExtraordinary)
					&& Double.doubleToLongBits(refractiveIndexOrdinary) == Double.doubleToLongBits(other.refractiveIndexOrdinary)
					&& Double.doubleToLongBits(transmission) == Double.doubleToLongBits(other.transmission);
			
	}

}
