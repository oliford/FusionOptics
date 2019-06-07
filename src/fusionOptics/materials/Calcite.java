package fusionOptics.materials;

import java.util.Arrays;

import fusionOptics.types.Material;

import algorithmrepository.exceptions.NotImplementedException;
import net.jafama.FastMath;

/** Calcite Crystal CaCO3
 * 
 * Mostly from [ http://refractiveIndex.info (CC Atributation Share-Alike) ]
 * 
 * 
 * 
 * @author oliford
 */
public class Calcite extends Material {
	
	public static final double C[][] = new double[][]{
			{ 0.8559, 0.0588, 0.8391, 0.141, 0.0009, 0.197, 0.6845, 7.005 },	
			{ 1.0856, 0.07897, 0.0988, 0.142, 0.317, 11.468, 0, 0 },
	};
	
	@Override
	public int getNAxes() { return 1;	}

	@Override
	public double getRefractiveIndex(int m, double wavelength,
			double temperature) {
		double l2 = wavelength*wavelength*1e6*1e6; //squared wavelength in um
		return FastMath.sqrt(1 + C[m][0]*l2/(l2-C[m][1]*C[m][1])
								+ C[m][2]*l2/(l2-C[m][3]*C[m][3])
								+ C[m][4]*l2/(l2-C[m][5]*C[m][5])
								+ C[m][6]*l2/(l2-C[m][7]*C[m][7]) //is a bit necessary at large l, for the O mode
								);
	}

	@Override
	public double getTransmission(int modeNumber, double wavelength,
			double temperature) {
		throw new NotImplementedException();
	}

	@Override
	public double getVerdetConstant(int modeNumber, double wavelen, double temperature) {
		throw new NotImplementedException();
		//return verdetConstFromBecquerelRelation(modeNumber, wavelen, temperature);
	}
	
	//Materials with no modifiable content are equal if they are the same type
	// but the static variables may have changed at compile and the comparison of hashcode might be against an old
	// code
	@Override
	public int hashCode() { return Arrays.deepHashCode(C); }
	@Override
	public boolean equals(Object obj) { return (obj != null) && obj instanceof Calcite; }	

}
