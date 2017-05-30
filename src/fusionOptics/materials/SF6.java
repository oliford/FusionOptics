package fusionOptics.materials;

import java.util.Arrays;
import fusionOptics.types.Material;
import algorithmrepository.LinearInterpolation1D;
import jafama.FastMath;

/** This is the SF6 Schott glass
 * 
 *  Refractive index formula from  at http://refractiveindex.info *  
 *  (CC Attributation Share-Alike)  
 *  
 * Transmission from http://www.uqgoptics.com/materials_glasses_schott_SF6.aspx (unknown license)
 * 
 * <Andrea, E3, IPP>
 */

public class SF6 extends Material {
	
	private final static double C[] = new double[]{ 
			1.72448482, 0.0134871947, 0.390104889, 0.0569318095, 1.04572858, 118.557185 };
	
	/* [ G.Westenberger "Verdet constant and its dispersion in optical glasses" doi:10.1117/12.48307 ] 
	 * SF6*/
	private final static double verdet_a = 1682.15e-9; // 1/T
	private final static double verdet_b = -88.7082e-20; //m²/T
	private final static double verdet_l0 = 156.4e-9; //m
	
	private final static double tL[] = new double[]{
			320.0000, 365.0000, 370.0000, 380.0000, 390.0000, 400.0000, 
			404.7000, 420.0000, 435.8000, 460.0000, 500.0000, 
			546.1000, 580.0000, 620.0000, 660.0000, 700.0000,
			1060.0000,1529.6000, 1970.1000, 2325.4000, 2500.0000 };
	private final static double tV[] = new double[]{ 
			0.00, 0.04, 0.12, 0.40, 0.65, 0.77,
			0.82, 0.91, 0.950, 0.975, 0.988,
			0.994, 0.994, 0.994, 0.994, 0.995,
			0.997, 0.991, 0.93, 0.80, 0.73}; 
	
	
	/** This data from  http://www.uqgoptics.com/materials_glasses_schott_SF6.aspx and is for 25mm thickness */
	LinearInterpolation1D transmittance = new LinearInterpolation1D(tL, tV);
	
	@Override
	public int getNAxes() {
		return 0; //implied by Schott catalog
	}
		@Override
	public double getRefractiveIndex(int modeNumber, double wavelength, double temperature) {
		double l2 = wavelength*wavelength*1e6*1e6; //wavelength^2 in um
		return FastMath.sqrt(1 + C[0]*l2/(l2-C[1]) + C[2]*l2/(l2-C[3]) + C[4]*l2/(l2-C[5]));
	}
		@Override
	public double getTransmission(int modeNumber, double wavelength, double temperature) {
		return FastMath.pow(transmittance.eval(wavelength/1e-9), 1.000 / 0.025); //convert to: per 1m
	}
		@Override
	public double getVerdetConstant(int modeNumber, double wavelen, double temperature) {		
			/* [ G.Westenberger "Verdet constant and its dispersion in optical glasses" doi:10.1117/12.48307 ] 
			 * SF6*/
			return Math.PI / wavelen * (verdet_a + verdet_b / (wavelen*wavelen - verdet_l0*verdet_l0));
	}
	
	//Materials with no modifiable content are equal if they are the same type
	// but the static variables may have changed at compile and the comparison of hashcode might be against an old
	// code

	@Override
	public int hashCode() { 
		long t; int r = 1;
		r = 31 * r + Arrays.hashCode(C);
		r = 31 * r + Arrays.hashCode(tL);
		r = 31 * r + Arrays.hashCode(tV);
		t = Double.doubleToLongBits(verdet_a); 		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(verdet_b); 		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(verdet_l0); 	r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	
	@Override
	public boolean equals(Object obj) { return (obj != null) && obj instanceof SF6; }
}
