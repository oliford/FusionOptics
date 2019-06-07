package fusionOptics.materials;

import java.util.Arrays;

import fusionOptics.types.Material;

import net.jafama.FastMath;


/** Fused Silica (SiO2). Like Crystal Quartz, but not birefringent. 
 * 
 * There is some confusion out there, over what is 'quartz', 'fused quartz' and 'fused silica' etc
 * 
 * In general it seems to be that 'fused' means it's a glass, not a crystal and is therefore not
 *   birefringent.
 *   
 * It seems mostly just a purity thing but [ http://www.azom.com/article.aspx?ArticleID=4766 ]
 * says this:
 * " This synthetic material, normally referred to as synthetic fused silica, has 
 *    better optical properties and is somewhat more expensive than the other type.
 *    In the UK, terms such as quartz, silica, fused quartz and fused silica tend to
 *    be used interchangeably. In the USA, quartz refers to material melted from 
 *    grains, silica refers to the synthetic material."
 */
public class FusedSilica extends Material {
	
	public static final double d[] = new double[]{
			1.0, 
			0.6961663, 0.0684043, 
			0.4079426, 0.1162414,
			0.8974794, 9.896161
		};
	
	public double magnetoOpticAnomaly = 0.70;

	@Override
	public int getNAxes() { return 0;	}

	@Override
	/** [ http://refractiveindex.info/?group=GLASSES&material=F_SILICA ] */
	public double getRefractiveIndex(int modeNumber, double wavelength,
			double temperature) {		

		double wavelengthSq = (wavelength*1e6)*(wavelength*1e6); //formula requires lambda in um

		return FastMath.sqrt( 
					d[0]
					+ d[1] * wavelengthSq / (wavelengthSq - d[2]*d[2]) 
					+ d[3] * wavelengthSq / (wavelengthSq - d[4]*d[4]) 
					+ d[5] * wavelengthSq / (wavelengthSq - d[6]*d[6]) 
				);
	}

	@Override
	public double getTransmission(int modeNumber, double wavelength,
			double temperature) {
		return 1;
	}

	@Override
	public double getVerdetConstant(int modeNumber, double wavelen, double temperature) {
		/*
		 * [ E.Munin "Faraday Effect And Energy Gap In Optical Matierals", J.Phys D (1992) ]
		 * 	claims that the Becquerel formula is valid for Fused Silica (specifically Dynasil 1001)
		 *  
		 *  in the range 457.9 - 632.8nm with Magnetooptic anomaly of 0.70#
		
		 */

		if(wavelen < 457.9e-9 || wavelen > 660e-9) // || wavelen > 632.8e-9)  - I've push this up because my normal use case (AUG IMSE) is just off the end
			System.err.println("WARNING: "+
					//throw new IllegalArgumentException(
					"Wavelength out of range for FusedSilica.getVerdetConstant()");

		//System.out.println("lambda: " + String.valueOf(wavelen));
		//System.out.println("verdet of SiO2: " + String.valueOf(verdetConstFromBecquerelRelation(modeNumber, wavelen, temperature, magnetoOpticAnomaly)));
		return verdetConstFromBecquerelRelation(modeNumber, wavelen, temperature, magnetoOpticAnomaly);
	}
	
	//Materials with no modifiable content are equal if they are the same type
	// but the static variables may have changed at compile and the comparison of hashcode might be against an old
	// code
	@Override
	public int hashCode() { 
		return 31 + Arrays.hashCode(d);
	}
	@Override
	public boolean equals(Object obj) { return (obj != null) && obj instanceof FusedSilica; }	

}
