package fusionOptics.materials;

import java.util.Arrays;
import fusionOptics.types.Material;
import algorithmrepository.LinearInterpolation1D;
import jafama.FastMath;

/** This is the LAK9 Schott glass
 * 
 *  Refractive index formula from  at http://refractiveindex.info *  
 *  (CC Attributation Share-Alike)  
 *  
 * Transmission from (http://www.us.schott.com/advanced_optics/english/knowledge-center/technical-articles-and-tools/abbe-diagramm.html#) (unknown license)
 * 
 * <Andrea, E3, IPP>
 */

public class LAK9 extends Material {
	
	public static final double e = 1.60217657e-19;
	public static final double c = 299792458;
	public static final double m_e = 9.10938291e-31;
	
	private final static double C[] = new double[]{
			1.46231905, 0.00724270156, 0.344399589, 0.0243353131, 1.15508372, 85.4686868 };
	
	private final static double tL[] = new double[]{ 
			320.0000, 334.1000, 350.0000, 365.0000, 370.0000,
			380.0000, 390.0000, 400.0000, 404.7000,	420.0000,
			435.8000, 460.0000, 500.0000, 546.1000,	580.0000,
			620.0000, 660.0000, 700.0000, 1060.0000,1530.0000,
			1970.0000, 2325.0000, 2500.0000};
	private final static double t[] = new double[]{
			0.020, 	0.200,	0.550,	0.782,	0.830,   
			0.890,	0.930,	0.950,	0.957,	0.970,   
			0.977,	0.984,	0.992,	0.994,	0.994,   
			0.995,	0.995,	0.996,	0.995,	0.966,   
			0.860,	0.420,	0.140}; 
	
	/** This data from  http://www.us.schott.com/advanced_optics/english/knowledge-center/technical-articles-and-tools/abbe-diagramm.html#
	 *  and is for 25mm thickness */
	
	LinearInterpolation1D transmittance = new LinearInterpolation1D(tL, t);
	
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
		
		double magnetoOpticAnomaly = 0.7;
		/*
		double a = 714.57e-9;
		double b = 20.8496e-20;
		double l0 = 106.5e-9;
		*/
		/*
		double a = 701.16e-9;
		double b = 12.2336e-20;
		double l0 = 106.5e-9;
		*/
		//double magnetoOpticAnomaly = -Math.PI*(a+b/(wavelen*wavelen-l0*l0))/(e/(2 * m_e * c)*wavelen*wavelen*getLinearDispersion(modeNumber,wavelen,temperature));
		//System.out.println(magnetoOpticAnomaly);
		
		return verdetConstFromBecquerelRelation(modeNumber, wavelen, temperature, magnetoOpticAnomaly);	}
	
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
	public boolean equals(Object obj) { return (obj != null) && obj instanceof LAK9; }	

}
