package fusionOptics.materials;

import java.util.Arrays;

import fusionOptics.types.Material;

/** Air
 * 
 *  Refractive index formula from  at http://refractiveindex.info *  
 *  (CC Attributation Share-Alike)  
 *  
 *  <Andrea, E3, IPP>
 */

public class Air extends Material {
	
	private final static double C[] = new double[]{
			0.05792105, 238.0185, 0.00167917, 57.362};

	@Override
	public double getRefractiveIndex(int modeNumber, double wavelength, double temperature) {
		double l2 = 1/(wavelength*wavelength*1e6*1e6); //wavelength^-2 in um
		return (1 + C[0]/(C[1]-l2) + C[2]/(C[3]-l2));
	}
	
	@Override
	public double getTransmission(int modeNumber, double wavelength, double temperature) {
		return 1;
	}

	@Override
	public int getNAxes() {
		return 0;
	}

	@Override
	public double getVerdetConstant(int modeNumber, double wavelen, double temperature) {
		return 0;
	}
	
	//Materials with no modifiable content are equal if they are the same type
	// but the static variables may have changed at compile and the comparison of hashcode might be against an old
	// code
	@Override
	public int hashCode() { 
		int result = 31 + Arrays.hashCode(C);
		return result;
	}
	@Override
	public boolean equals(Object obj) { return (obj != null) && obj instanceof Air; }	

}
