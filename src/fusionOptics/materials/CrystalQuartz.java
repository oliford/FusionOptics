package fusionOptics.materials;

import fusionOptics.Util;

/** Crystal Quartz / Fused Quartz - The birefringent form of Fused Silica (SiO2) 
 * 
 * There is some confusion out there, over what is 'quartz', 'fused quartz' and 'fused Silica' etc
 * 
 * In general it seems to be that 'fused' means it's a glass, not a crystal and is therefore not
 *   birefringent. 
 * 
 * 
 * It seems mostly just a purity thing but [ http://www.azom.com/article.aspx?ArticleID=4766 ]
 * says this:
 * " This synthetic material, normally referred to as synthetic fused silica, has 
 *    better optical properties and is somewhat more expensive than the other type.
 *    In the UK, terms such as quartz, silica, fused quartz and fused silica tend to
 *    be used interchangeably. In the USA, quartz refers to material melted from 
 *    grains, silica refers to the synthetic material."
 * 
 * The main source for this data is a brochure for someone selling the stuff:
 * 		[ TYDEX J.S.Co, "Synthetic crystal quartz" Commercial Brochure, 
 *           http://www.tydexoptics.com/pdf/Crystal_quartz.pdf ]
 * 	
 * The Verdet constant at 589.3nm seems to be 0.01664 min cm^-1 gauss^-1
 *	 [S.Ramaseshan "Determination Of The Magnetooptic Anomaly Of Some Glasses" 1946]
 *  277.33 deg m^-1 T^-1, 4.84 rad m^-1 T^-1
 *  

	Temperature dependance (in 0 - 400 'C) is pretty small with 
	  |dn/dT| < 1.5e-5 Kelvin^-1 in the wavelength range 500 - 1500nm
	[  T Toyoda 1983
	  "The Temperature Dependence of the Refractive Indices of Fused Silica and Crystal Quartz" ]
	  
	  
	There is also this paper for 'fused silica', which gives a refractive index formula, but only 1  
	   (no extraordinary ray). It does however, give roughly (slightly lower) temperature dependence.
	   
*/
public class CrystalQuartz extends WavelengthDependent {
	
	 public CrystalQuartz() {
		super(new double[]{   //wavelength / m
					185e-9, 194e-9, 204e-9, 219e-9, 231e-9, 
					243e-9, 263e-9, 291e-9, 340e-9, 405e-9, 
					589e-9, 1083e-9, 1800e-9, 2500e-9, 3000e-9 
				},new double[][]{ 
					{ //ordinary index 
						1.676, 1.660, 1.643, 1.625, 1.614, 
						1.605, 1.593, 1.581, 1.567, 1.557, 
						1.544, 1.534, 1.524, 1.512, 1.500, 
					},{ //extraordinary index
						1.690, 1.673, 1.656, 1.637, 1.626, 
						1.617, 1.604, 1.591, 1.577, 1.567, 
						1.553, 1.543, 1.532, 1.520, 1.507, 
					}
				},
				273, // Min temp. This is 0'C to 400'C, which gives a 0.003 change in index 
				673, // Max Temp
				1.00); //magnetic anomaly for calculating verdet constant

	 }

	 public static void main(String[] args) {
		 System.out.println(Util.calcWaveplateFullWaveThickness(new CrystalQuartz(), 600e-9));
	 }

	 @Override
	 public double getVerdetConstant(int modeNumber, double wavelen, double temperature) {
		 /*
		  * [ E.Munin "Faraday Effect And Energy Gap In Optical Matierals", J.Phys D (1992) ]
		  * 	claims that the Becquerel formula is valid for Fused Silica (specifically Dynasil 1001)
		  *  
		  *  in the range 457.9 - 632.8nm with Magnetooptic anomaly of 0.70#
		  *  
		  *  At 630nm, this gives a verdet constant of about 150 °/T/m, which is much lower
		  *  than the often published ~ 270
		  *  
		  */
		 double magnetoOpticAnomaly = 0.70;

		 if(wavelen < 457.9e-9 || wavelen > 632.8e-9)
			System.err.println("WARNING: "+
			//throw new IllegalArgumentException(
			 		 "Wavelength out of range for Quartz.getVerdetConstant()");

		 return verdetConstFromBecquerelRelation(modeNumber, wavelen, temperature, magnetoOpticAnomaly);
	 }
}
