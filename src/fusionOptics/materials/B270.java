package fusionOptics.materials;

import java.util.Arrays;

import fusionOptics.types.Material;

import algorithmrepository.LinearInterpolation1D;
import algorithmrepository.exceptions.NotImplementedException;
import net.jafama.FastMath;

/** This is the B270 Schott glass
 * 
 *  Refractive index and abbe from Edmund Optics website  
 *  
 * Transmission unknown
 */
public class B270 extends IsotropicLinearDispersiveGlass {

	// [Edmund optics]
	public final static double refractiveIndex_d = 1.523;
	public final static double abbeNumber_d = 58.50;
	
	public B270() {
		super(refractiveIndex_d, abbeNumber_d, 1.0);
	}

	@Override
	public int getNAxes() {
		return 0; //implied by Schott catalog
	}

	@Override
	public double getTransmission(int modeNumber, double wavelength, double temperature) {
		return 1.0; //unknown
	}

	//@Override
	//public double getVerdetConstant(int modeNumber, double wavelen, double temperature) {
		/*
		 * Unknown, using the same magneto optic anomaly of BK7
		 */
	//	double magnetoOpticAnomaly = 0.72;
		
	//	if(wavelen < 457.9e-9 || wavelen > 660e-9) //32.8e-9)
	//		throw new IllegalArgumentException("Wavelength out of range for BK7.getVerdetConstant()");
		
	//	return verdetConstFromBecquerelRelation(modeNumber, wavelen, temperature, magnetoOpticAnomaly);
	//}
	
	@Override
	public double getVerdetConstant(int modeNumber, double wavelen,
			double temperature) {
		
		if(wavelen < 640e-9 || wavelen > 670e-9){
			throw new RuntimeException("Not implemented. This was only measured at 656nm");
		}
		
		return 4.6913; // measured in Garching, Aug/Sept 2017, error < 0.2  wavelength ~651nm 
	}
	
	//Materials with no modifiable content are equal if they are the same type
	// but the static variables may have changed at compile and the comparison of hashcode might be against an old
	// code
	@Override
	public int hashCode() { 
		int result = 31 + new Double(refractiveIndex_d).hashCode();
		result = 31 * result + new Double(abbeNumber_d).hashCode();
		return result;
	}
	@Override
	public boolean equals(Object obj) { return (obj != null) && obj instanceof B270; }	

}
