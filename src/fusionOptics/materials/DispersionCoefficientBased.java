package fusionOptics.materials;

import java.util.Arrays;

import jafama.FastMath;
import fusionOptics.types.Material;

/** Based on the simpler dispersion formula (Sellmeier equation?)
 *  n² = 1 + c0 l² / (l² - c3) + c1 l² / (l² - c4) + c2 l² / (l² - c5)  
 * 
 * @author oliford
 *
 */
public abstract class DispersionCoefficientBased extends Material {
	
	private double C[][];

	public DispersionCoefficientBased(double C[][]) {
		this.C = C;
	}
	

	@Override
	public int getNAxes() { return C.length-1; }

	@Override
	public double getRefractiveIndex(int modeNumber, double wavelength, double temperature) {
		double l2 = wavelength*wavelength*1e6*1e6; //wavelength^2 in um
		return FastMath.sqrt(1 + C[modeNumber][0]*l2/(l2-C[modeNumber][3]) + C[modeNumber][1]*l2/(l2-C[modeNumber][4]) + C[modeNumber][2]*l2/(l2-C[modeNumber][5]));
	}

	@Override
	public double getTransmission(int modeNumber, double wavelength,
			double temperature) {
		return 1.0;
	}

	@Override
	public double getVerdetConstant(int modeNumber, double wavelen,
			double temperature) {
		return 0.0;
	}

	@Override
	public int hashCode() { 
		long t; int r = 1;
		r = 31 * r + Arrays.deepHashCode(C);
		
		//t = Double.doubleToLongBits(verdet_l0); 		r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	@Override
	public boolean equals(Object obj) { return (obj != null) && this.getClass().isInstance(obj); }

}
