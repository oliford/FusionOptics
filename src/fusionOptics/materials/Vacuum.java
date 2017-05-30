package fusionOptics.materials;

import fusionOptics.types.Material;

/** Ideal Vacuum, n = 1.000, v = 0, T = 1.
 * Actually, you don't need to use this because Tracer 
 * will treat any null medium as vacuum */
public class Vacuum extends Material {

	@Override
	public double getRefractiveIndex(int modeNumber, double wavelength, double temperature) {
		return 1.0;
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

	@Override
	public int hashCode() { return 0; }  //must be null-like hashcode

	// As far as the ray-tracer is concerned, Vacuum and 'null' are the same thing
	@Override
	public boolean equals(Object obj) { return (obj == null || obj instanceof Vacuum);	}
}
