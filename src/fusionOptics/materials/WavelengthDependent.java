package fusionOptics.materials;

import java.util.Arrays;

import algorithmrepository.ExtrapolationMode;
import algorithmrepository.Interpolation1D;
import algorithmrepository.InterpolationMode;
import fusionOptics.types.Material;
import net.jafama.FastMath;
import uk.co.oliford.jolu.BinaryMatrixFile;

/** Materials with wavelength dependent refractive indices given as a linear interpolatable
 * function.
 * 
 * This will use the Becquerel relation for the verdet constant, which will in general be wrong. 
 * 
 * @author oliford
 */
public class WavelengthDependent extends Material {
	
	private Interpolation1D interps[];
	public double minTemp, maxTemp;
	private double magneticAnomaly;
	
	/** 	Create of wavelength dependwnt refractive index interpolators.
	 *	@param magneticAnomaly Used in calculation of verdet constant 
	 *			(see Material.verdetConstFromBecquerelRelation())  
	 */
	public WavelengthDependent(double wavelen[], double n[][], double magneticAnomaly, double minTemp, double maxTemp) {
		this.minTemp = minTemp;
		this.maxTemp = maxTemp;
		this.magneticAnomaly = magneticAnomaly; 

		interps = new Interpolation1D[n.length];
		for(int i=0; i < n.length; i++) {
			interps[i] = new Interpolation1D(wavelen, n[i], InterpolationMode.LINEAR, ExtrapolationMode.EXCEPTION);
		}
	}
	
	@Override
	public double getRefractiveIndex(int modeNumber, double wavelength, double temperature) {
			return interps[modeNumber].eval(wavelength);
	}
	
	@Override
	public double getTransmission(int modeNumber, double wavelength, double temperature) {
		return 1;
	}

	@Override
	public int getNAxes() { return interps.length - 1; }

	@Override
	public double getVerdetConstant(int modeNumber, double wavelen,
			double temperature) {
		throw new UnsupportedOperationException();
	}


	@Override
	public int hashCode() {
		long t; int r = 1;		
		for(Interpolation1D interp : interps){
			r = 31 * r + Arrays.hashCode(interp.getX());
			r = 31 * r + Arrays.hashCode(interp.getF());
		}			
		t = Double.doubleToLongBits(magneticAnomaly);	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(maxTemp);			r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(minTemp);			r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		WavelengthDependent other = (WavelengthDependent) obj;
		if(interps.length != other.interps.length)
			return false;
		for(int i=0; i < interps.length; i++){
			if(!Arrays.equals(interps[i].getX(), other.interps[i].getX())) return false;
			if(!Arrays.equals(interps[i].getF(), other.interps[i].getF())) return false;
		}
		return Double.doubleToLongBits(magneticAnomaly) == Double.doubleToLongBits(other.magneticAnomaly)
				&& Double.doubleToLongBits(maxTemp) == Double.doubleToLongBits(other.maxTemp)
				&& Double.doubleToLongBits(minTemp) == Double.doubleToLongBits(other.minTemp);
	}
}
