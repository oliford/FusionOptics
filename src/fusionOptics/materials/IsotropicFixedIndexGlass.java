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
public class IsotropicFixedIndexGlass extends Material {
	private double refractiveIndex;
	private double transmission; 
	
	/** Create ideal glass with given refractive index and perfect transmission (no absorption) */
	public IsotropicFixedIndexGlass(double refractiveIndex) {
		this.refractiveIndex = refractiveIndex;
		this.transmission = 1.0;
	}
	
	/** Create ideal glass with given refractive index and given transmission. */
	public IsotropicFixedIndexGlass(double refractiveIndex, double transmission) {
		this.refractiveIndex = refractiveIndex;
		this.transmission = transmission;
	}
	
	public void setRefractiveIndex(double refractiveIndex) { this.refractiveIndex = refractiveIndex; }
	
	public void setTransmission(double transmission) { this.transmission = transmission; }
	
	@Override
	public int getNAxes() { return 0; }

	@Override
	public double getRefractiveIndex(int modeNumber, double wavelength, double temperature) { return refractiveIndex; }

	@Override
	public double getTransmission(int modeNumber, double wavelength, double temperature) { return transmission; }

	@Override
	public double getVerdetConstant(int modeNumber, double wavelen, double temperature) { 
		throw new UnsupportedOperationException();
		}


	@Override
	public int hashCode() {
		long temp = Double.doubleToLongBits(refractiveIndex);
		int  result = 31 + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(transmission);
		result = 31 * result + (int) (temp ^ (temp >>> 32));
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		IsotropicFixedIndexGlass other = (IsotropicFixedIndexGlass) obj;
		return Double.doubleToLongBits(refractiveIndex) == Double.doubleToLongBits(other.refractiveIndex)
				&& Double.doubleToLongBits(transmission) == Double.doubleToLongBits(other.transmission);		
	}

}
