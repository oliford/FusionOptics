/**
 * Copyright 2011 Oliver Ford
 *
 * This file is part of the minerva-optics 'RayTracer'.
 *
 *   RayTracer is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   RayTracer is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RayTracer.  If not, see <http://www.gnu.org/licenses/>.
 *   
 *   @author oliford <codes<at>oliford.co.uk>
 */
package fusionOptics.materials;

import binaryMatrixFile.BinaryMatrixFile;
import oneLiners.OneLiners;

/** Alpha Barium Borate
 * 
 * Alpha means high temp phase
 * 
 * @author oliford
 */
public class AlphaBariumBorate extends SellmeierBased {
	
	public AlphaBariumBorate() {
		/* SellmeierBased is:
		 * 
		 * F = (T - T0) * (T + T0 + 546);
		 *	n^2 = A1 + (A2 + B1*F) / (L*L - (A3+B2*F)*(A3+B2*F)) + B3*F - A4*L*L;
		 *
		 *
		 */
		super(new double[][]{
				/* For the Sellmeier coeffs A1,A2,A3,A4
				 * 
				 * [ http://www.unitedcrystals.com/BBOProp.html United Crystals website ]
				 * Beta:	2.7359, 0.01878, 0.01822, 0.01354 
				 * 			2.3753, 0.01224, 0.01667, 0.01516
				 * 
				 * [ http://www.agoptics.com/Birefrigent-Crystal/alpha-BBO.htm ]
				 * Alpha:	2.7471, 0.01878, 0.01822, 0.01354 
 				 * 			2.3753, 0.01224, 0.01667, 0.01516 	<--- Looks like these are from Beta
				 * 
				 * [ ... via John Howard ANU ]
				 * Beta:	2.7359, 0.01878, 0.01822, 0.01354 	[Unknown]
				 * 			2.3753, 0.01224, 0.01667, 0.01516
				 * 
				 * Alpha:	2.7471, 0.01878, 0.01822, 0.01354])	[Scott Silburn]
				 *			2.3174, 0.01224, 0.01667, 0.01516])
				 *
				 * 
				 * The 'Beta' values agree with JH's measurements for alpha
				 * 
				 * 
				 * FIXME: No temperature dependence yet. 
				 *   dN/dT is given usually.
				 *   Here, B3 = dN^2/dF = 2.N.dN/dF
				 * 
				 *   but F is a bit of a weird function of T
				 *   work it out later
				 *  
				 */
					{ //ordinary						
						//2.7471, 0.01878, Math.sqrt(0.01822), 0.01354, //A1,A2,A3,A4
						2.7359, 0.01878, Math.sqrt(0.01822), 0.01354, 
						0, 0, //B1, B2
						0, //B3, we sould be able to do
						27, //T0 / deg.C - Not given by ref, although it does give dn/dT
						
					},{ //extraordinary
						//2.3174, 0.01224, Math.sqrt(0.01667), 0.01516, //A1,A2,A3,A4
						2.3753, 0.01224, Math.sqrt(0.01667), 0.01516,
						0, 0, //B1, B2 
						0, //B3, we sould be able to do
						27 //T0 / deg.C					
					},				
				},
				Double.NaN ); // no data
		
		
		this.minWavelength = 196e-9;
		this.maxWavelength = 2200e-9;
		
		
	}
	
	@Override
	public double getRefractiveIndex(int modeNumber, double wavelength, double temperature) {
		
		if(wavelength < 450e-9 || wavelength > 710e-9)
			throw new IllegalArgumentException("Experimentally, the αBBO Sellmeier are completely wrong outside 450nm < λ < 710nm.");
			
		return super.getRefractiveIndex(modeNumber, wavelength, temperature);
	}
	
	/* Also listed as just BBO in Handbook of Optics Chapter33 - Optical and physical props.. 
	 * Originally:
	 * [ D. Eimerl, "Optical, Mechanical, and Thermal Properties of Barium Borate" J. Appl. Phys. 62:1968 – 1983 (1987). ]

	 * Transparency (um) UV=0.205, IR=3.0, 
	 *  
	 *  no=1.540, ne=1.655,
	 *  
	 * at l=0.4047um: dno/dT=-16.6, dne/dT=-9.8
	 * at l=0.5790um  -16.4, -9.4
	 * at l=1.014um  -16.8 	-8.8
	 * 
	 */
	
	@Override
	public double getTransmission(int modeNumber, double wavelength,
			double temperature) {
		return 1.0;
	}
}
