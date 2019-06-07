package fusionOptics.materials;

import java.util.Arrays;

import fusionOptics.types.Material;

import binaryMatrixFile.BinaryMatrixFile;
import algorithmrepository.exceptions.NotImplementedException;
import net.jafama.FastMath;

/** For materials with a refractive index calculated from the Sellmeier equation */ 
public class SellmeierBased extends Material {


	protected double coeffs[][];
	protected double minWavelength = 0 , maxWavelength = Double.POSITIVE_INFINITY;
	protected double magneticAnomaly;
	
	public SellmeierBased(double coeffs[][], double magneticAnomaly) {
		this.coeffs = coeffs;
		this.magneticAnomaly = magneticAnomaly;
	}

	/** Sellmeier equation for wavelength and temperature dependence of refractive indices
	 * 
	 *	F = (T - T0) * (T + T0 + 546);
	 *	n^2 = A1 + (A2 + B1*F) / (L*L - (A3+B2*F)*(A3+B2*F)) + B3*F - A4*L*L;
	 * 
	 * [Equation 1 of Optical and Quantum Electronics 16 (1984) Short Communication,
	 *  "A Temperature Dependent Dispersion Equation For Congruently Grown Lithium Niobate" ] 
	 * 
	 * @param wavelen		in m
	 * @param temperature	in K
	 */
	private final double sellmeierEquation(double wavelength, double temperature, int modeNumber){
		if(wavelength < minWavelength || wavelength > maxWavelength)
			return Double.NaN;
		
		double L = wavelength / 1e-6; //in um
		double T = temperature - 273.15; // in def.C
			
		double A1 = coeffs[modeNumber][0];
		double A2 = coeffs[modeNumber][1];
		double A3 = coeffs[modeNumber][2];
		double A4 = coeffs[modeNumber][3];
		double B1 = coeffs[modeNumber][4];
		double B2 = coeffs[modeNumber][5];
		double B3 = coeffs[modeNumber][6];
		double T0 = coeffs[modeNumber][7];
		
		double F = (T - T0) * (T + T0 + 546);
		
		double n2 = A1 + (A2 + B1*F) / (L*L - (A3+B2*F)*(A3+B2*F)) + B3*F - A4*L*L;
		
		return FastMath.sqrt(n2);
	}
	
	/* Useful code if we want to turn this into the interpolation based 
	if(refractiveIndexVsWavelen[i] == null && sellmeierCoeffs[i] != null){
		//fill in data from Sellmeier equation, on the range for which it's valid
		double lp0 = FastMath.log10(sellmeierCoeffs[i][0][0]);
		double lp1 = FastMath.log10(sellmeierCoeffs[i][0][1]);
		int nL = 500;
		double l[] = new double[nL];
		for(int j=0; j < nL; j++){
			double exp = lp0 + j*(lp1-lp0)/((double)nL-1.0);
			l[j] = Math.pow(10, exp);
		}
		
		refractiveIndexVsWavelen[i] = new double[][]{ l, null, null };
		for(int nIdx=0; nIdx < 2; nIdx++){
			refractiveIndexVsWavelen[i][nIdx+1] = new double[nL];
			for(int j=0; j < nL; j++){
				refractiveIndexVsWavelen[i][nIdx+1][j] = sellmeierEquation(i, l[j], 300, (nIdx == 1));
			}
		}
		BinaryMatrixFile.mustWrite("/tmp/blah.bin", refractiveIndexVsWavelen[i], true);
	}*/
	
	@Override
	public double getRefractiveIndex(int modeNumber, double wavelength, double temperature) {
		return sellmeierEquation(wavelength, temperature, modeNumber);
		
	}

	@Override
	public int getNAxes() { return coeffs.length - 1;	}

	@Override
	public double getTransmission(int modeNumber, double wavelength, double temperature) {
		throw new NotImplementedException();
	}

	@Override
	public double getVerdetConstant(int modeNumber, double wavelen, double temperature) {
		throw new NotImplementedException();
	}
	

	@Override
	public int hashCode() {
		long t; int r = 1;
		r = 31 * r + Arrays.hashCode(coeffs);		
		t = Double.doubleToLongBits(magneticAnomaly); 	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(maxWavelength);		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(minWavelength);		r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		SellmeierBased other = (SellmeierBased) obj;
		return Arrays.deepEquals(coeffs, other.coeffs)
				&& Double.doubleToLongBits(magneticAnomaly) == Double.doubleToLongBits(other.magneticAnomaly)
				&& Double.doubleToLongBits(maxWavelength) == Double.doubleToLongBits(other.maxWavelength)
				&& Double.doubleToLongBits(minWavelength) == Double.doubleToLongBits(other.minWavelength);
	}
}
