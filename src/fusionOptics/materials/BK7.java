package fusionOptics.materials;

import java.util.Arrays;

import algorithmrepository.ExtrapolationMode;
import algorithmrepository.Interpolation1D;
import algorithmrepository.InterpolationMode;
import fusionOptics.types.Material;

import net.jafama.FastMath;

/** This is the BK7 Schott glass
 * 
 *  Refractive index formula from  at http://refractiveindex.info *  
 *  (CC Attributation Share-Alike)  
 *  
 * Transmission from http://www.uqgoptics.com/materials_optical_schottBK7.aspx (unknown license)
 */
public class BK7 extends Material {

	private final static double C[] = new double[]{
		1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653 };
	
	private final static double tL[] = new double[]{ 
			310.0000, 320.0000, 334.1000, 350.0000, 365.0000,
			370.0000, 380.0000, 390.0000, 400.0000, 404.7000,
			420.0000, 435.8000, 460.0000, 500.0000, 546.1000,
			580.0000, 620.0000, 660.0000, 700.0000, 1060.0000,
			1529.6000, 1970.1000, 2325.4000 };
	private final static double t[] = new double[]{
			0.070,   0.350,   0.770,   0.930,   0.969,   
			0.974,   0.980,   0.989,   0.991,   0.993,   
			0.993,   0.994,   0.994,   0.996,   0.996,   
			0.996,   0.997,   0.997,   0.998,   0.998,   
			0.985,   0.850,   0.570, }; 
	
	
	/** This data from  http://www.uqgoptics.com/materials_optical_schottBK7.aspx and is for 25mm thickness */
	Interpolation1D transmittance = new Interpolation1D(tL, t, InterpolationMode.LINEAR, ExtrapolationMode.EXCEPTION);
	
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
		/*
		 * [ E.Munin "Faraday Effect And Energy Gap In Optical Matierals", J.Phys D (1992) ]
		 * 	claims that the Becquerel formula is valid for Schott BK7 
		 *  in the range 457.9 - 632.8nm with Magnetooptic anomaly of 0.72
		 */
		double magnetoOpticAnomaly = 0.72;
		
		if(wavelen < 457.9e-9 || wavelen > 660e-9) //32.8e-9)
			throw new IllegalArgumentException("Wavelength out of range for BK7.getVerdetConstant()");
		
		return verdetConstFromBecquerelRelation(modeNumber, wavelen, temperature, magnetoOpticAnomaly);
	}
	
	//Materials with no modifiable content are equal if they are the same type
	// but the static variables may have changed at compile and the comparison of hashcode might be against an old
	// code
	@Override
	public int hashCode() { 
		int result = 31 + Arrays.hashCode(tL);
		result = 31 * result + Arrays.hashCode(t);
		result = 31 * result + Arrays.hashCode(C);
		return result;
	}
	@Override
	public boolean equals(Object obj) { return (obj != null) && obj instanceof BK7; }	

}
