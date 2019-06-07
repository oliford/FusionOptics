package fusionOptics.materials;

import java.util.Arrays;

import fusionOptics.types.Material;

import algorithmrepository.LinearInterpolation1D;
import net.jafama.FastMath;

/** This is the SFL6 Schott glass, which I'm assuming is the 
 * same as the 'SLF6' mentioned in the AUG MSE system design.
 * 
 * The Schott catalog implies there is almost no birefringence.
 * 
 *  According to the 'Handbook of Optics', SF6 is:
 *    Dense Flint, 805254 SF6, 33% SiO2, 62% PbO, 5% K2O 
 *    
 *  SFL6 has the same 'code' (805254) which gives the wavelength and
 *  dispersion (I think), but is much denser/heavier. Presumably this
 *  means the composition is somewhat different in order to get the 
 *  same index?
 *  
 *  So the refractive formula is taken from the 'SF6' entry
 *  at http://refractiveindex.info
 *  
 *  (CC Attributation Share-Alike)  
 *  
 */
public class SchottSFL6 extends Material {

	private final static double C[] = new double[]{
			1.72448482, 0.0134871947, 0.390104889, 0.0569318095, 1.04572858, 118.557185 
		};
	
	/* [ G.Westenberger "Verdet constant and its dispersion in optical glasses" doi:10.1117/12.48307 ] 
	 * SF L6 and SF L56*/
	private final static double verdet_a = 273.20e-9; // 1/T
	private final static double verdet_b = -7.4974e-20; //m²/T
	private final static double verdet_l0 = 154.3e-9; //m
	
	//private final static double t = 8.74767363010858e-08; // m^-1. Shcott catalog gives 0.850 per 10mm at 400nm
	/** This data from http://www.uqgoptics.com/materials_glasses_schott_SFL6.asp, which matches the Schott single value at 400nm 
	 * This data is per 5mm.
	 */
	private final static double tL[] = new double[]{ 
		365.0000, 370.0000, 380.0000, 390.0000, 400.0000,
		404.7000, 420.0000, 435.8000, 460.0000, 500.0000,
		546.1000, 580.0000, 620.0000, 660.0000, 700.0000,
		1060.0000, 1529.6000, 1970.1000, 2325.4000,
	};
	private final static double tV[] = new double[]{
					0.220, 0.450, 0.760, 0.870, 0.920,
					0.930, 0.956, 0.970, 0.979, 0.987,
					0.994, 0.996, 0.997, 0.998, 0.998,
					0.998, 0.998, 0.989, 0.965, }; 

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
		return FastMath.pow(transmittance.eval(wavelength/1e-9), 1.000 / 0.005); //convert to: per 1m
	}

	@Override
	public double getVerdetConstant(int modeNumber, double wavelen, double temperature) {
		/** V < 0.05rad/m/T, according to:
		 *  [ Mentioned in Alcator C-Mod MSE 'Mini-proposal', but no original source.
		 *     http://www.psfc.mit.edu/research/alcator/miniproposals/431.pdf ] */
		//return 0.05e-3; //radians per T
		/** FIXME: This needs checking */

		/* [ G.Westenberger "Verdet constant and its dispersion in optical glasses" doi:10.1117/12.48307 ] 
		 * SF L6 and SF L56*/
		return Math.PI / wavelen * (verdet_a + verdet_b / (wavelen*wavelen - verdet_l0*verdet_l0));
		// That gives V = 0.42 rad/m/T at 650nm which IS significant enough to cause us problems.
	}

	@Override
	public int hashCode() { 
		long t; int r = 1;
		r = 31 * r + Arrays.hashCode(C);
		r = 31 * r + Arrays.hashCode(tL);
		r = 31 * r + Arrays.hashCode(tV);
		t = Double.doubleToLongBits(verdet_a); 		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(verdet_b); 		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(verdet_l0); 		r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	
	@Override
	public boolean equals(Object obj) { return (obj != null) && obj instanceof SchottSFL6; }

}
