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

/** Beta Barium Borate
 * 
 * Beta means low temp phase
 * 
 * @author oliford
 */
public class BetaBariumBorate extends SellmeierBased {
	
	public BetaBariumBorate() {
		/* SellmeierBased is:
		 * 
		 * F = (T - T0) * (T + T0 + 546);
		 *	n^2 = A1 + (A2 + B1*F) / (L*L - (A3+B2*F)*(A3+B2*F)) + B3*F - A4*L*L;
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
				 * [ refractiveindex.info, originally Handbook of Optics 3rd ed, Vol 4 2009 ]
				 * Beta: 	2.7405, 0.0184,  0.0179,  0.0155	<== Doesn't match, looks closer to alpha
				 * 			2.3730,	0.0128,  0.0156,  0.0044   <-- looks closer to beta than alpha, but doesn't match
				 * 
				 * [  http://www.lambdaoptics.com/Beta-Barium-Borate-β-BaB2O4-BBO.html ]
				 * Beta: 	2.7359, 0.01878, 0.01822, 0.01354 
 							2.3753,	0.01224, 0.01667, 0.01516
				 * 
				 * 
				 * FIXME: No temperature dependence yet. 
				 *   dN/dT is given usually.
				 *   See notes at bottom of AlphaBariumBorate
				 *   Here, B3 = dN^2/dF = 2.N.dN/dF
				 * 
				 *   but F is a bit of a weird function of T
				 *   work it out later
				 *  
				 */
					{ //ordinary						
						2.7359, 0.01878, Math.sqrt(0.01822), 0.01354, //A1,A2,A3,A4
						0, 0, //B1, B2
						0, //B3, we sould be able to do
						27, //T0 / deg.C - Not given by ref, although it does give dn/dT
						
					},{ //extraordinary
						2.3753, 0.01224, Math.sqrt(0.01667), 0.01516, //A1,A2,A3,A4
						0, 0, //B1, B2 
						0, //B3, we sould be able to do
						27 //T0 / deg.C					
					},				
				},
				Double.NaN ); // no data
		
		this.minWavelength = 196e-9;
		this.maxWavelength = 2200e-9;
	}
	
	
}
