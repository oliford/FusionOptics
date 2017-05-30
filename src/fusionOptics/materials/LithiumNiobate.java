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

import algorithmrepository.exceptions.NotImplementedException;

public class LithiumNiobate extends SellmeierBased {
	
	public LithiumNiobate() {
		super(new double[][]{
					/* [ Optical and Quantum Electronics 16 (1984) Short Communication 
					  A Temperature Dependent Dispersion Equation For Congruently Grown Lithium Niobate ]
					*/
					{ //ordinary
						4.9048, 0.11775, 0.21802, 0.027153, //A1,A2,A3,A4
						2.2314e-8, -2.9671e-8, 2.1429e-8, //B1,B2,B3
						24.5, //T0 / deg.C
					},{ //extraordinary
						4.5820, 0.09921, 0.21090, 0.021940, //A1,A2,A3,A4
						5.2716e-8, -4.9143e-8, 2.2971e-7, //B1,B2,B3
						24.5 //T0 / deg.C					
					},				
				},
				1.00 ); //This is the 'magnetic anomaly' for calculating the Verdet constant from the Becquerel formula
		
				//
		
		this.minWavelength = 0.4e-6;
		this.maxWavelength = 3.1e-6;
	}
	
	/* Lithium Niobate LiNbO3, very good summary info here: 
	[ U.Schlarb Refractive "Indices Of Lithium Niobate As A Function Of Temperature
	    Wavelength And Composition - A Generalized Fit", 1993, Phys Rev B. v48 n21 
	    10.1103/PhysRevB.48.15613
 	*/
	
	@Override
	public double getVerdetConstant(int modeNumber, double wavelen,
			double temperature) {

		/* I had written here "this seem to be OK for LithiumNiobate", but I've no idea what data I checked
		/ that against, or where it came from.
		 * 
		 * I have now found this:
		 * [ S.Kase "Optical Absorption And Interband Faraday Rotation In LiTaO3 And LiNbO3"
		 *   Assuming 'interband' faraday rotation isn't anything special.
		 *   This gives for 2eV photos (600nm) around 5e-5/ deg / cm / G which is 
		 *   	~ 0.5 deg / cm / T
		 *      0.05 deg / mm / T
		 *      0.005 deg /mm / 100mT
		 *   
		 *   Which is about the same as Quartz, apparently.
		 *   
		 *   The becquerel with anomaly = 1 gives 0.34 deg / mm / 100mT
		 *   which is much bigger.
		 *   
	*/
		throw new NotImplementedException();
	}
	
	@Override
	public double getTransmission(int modeNumber, double wavelength,
			double temperature) {
		return 1;
	}
}
