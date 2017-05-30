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
package fusionOptics.types;

import fusionOptics.materials.BK7;
import fusionOptics.materials.CrystalQuartz;
import fusionOptics.materials.LithiumNiobate;
import fusionOptics.materials.SchottSFL6;

/** Each instantiated object describes properties of a type of material (e.g. Quartz etc).
 * Each class might always initialised to only one matieral, or might provide
 * a whole series of materials from a common source or basis. e.g. from a 
 * manufacturer's glass catalogue.
 */ 
public abstract class Material {

	
	/** Returns the number of different directions/modes of the material
	 * 0=isotropic, 1=uniaxial, 2=biaxial */
	public abstract int getNAxes();
	
	/** Returns the refractive indicies along the given principal directions */
	public abstract double getRefractiveIndex(int modeNumber, double wavelength, double temperature);
	
	/** Returns the transmission coefficient m^-1 of the given mode at the given wavelength and temperature */ 
	public abstract double getTransmission(int modeNumber, double wavelength, double temperature);
	
	/** Returns the verdet constant for the given material at the given wavelength.
	 * 
	 * Implementors must override this, but can just call verdetConstFromBecquerelRelation
	 * and multiply by the appropriate magnetic anomaly factor.
	 * 
	 * However, verdetConstFromBecquerelRelation() doesn't seem to be at all valid for 
	 * some materials (e.g. SFL6), so the call has to be explicit in each case.	 *
	 * 
	 *  The definition of direction is taken from 
	 *  [ G.Westenberger "Verdet constant and its dispersion in optical glasses" doi:10.1117/12.48307 ]
	 *  
	 *  "The sign of V — and thus of the rotation angle — is positive (negative), if the rotation is clockwise 
	 *  (counter-clockwise) looking parallel to the vector of the magnetic flux density irrespective of the direction 
	 *  of the propagation of the wave being parallel or antiparallel to . Positive Verdet constants correspond to 
	 *  diamagnetic materials, whereas V is negative for paramagnetics."
	 *  
	 *  The paper then shows V for 'typical' glasses, of which all are +ve
	 *  

	 *  
	 * 
	 * @param material	
	 * @param wavelen	in m
	 * @param extraordinary		true for extraordinary wave, false for ordinary
	 * @return
	 */
	public abstract double getVerdetConstant(int modeNumber, double wavelen, double temperature);
	
	/** Returns the verdet constant for the given material at the given wavelength
	 * from the 'Becquerel' formula with a magneto-optic anomaly of 1.
	 * 
	 * V = lambda * dn/dLambda *(e / 2mc)
	 * 
	 * 
	 * 
	 * @param material	
	 * @param wavelen	in m
	 * @param extraordinary		true for extraordinary wave, false for ordinary
	 * @return
	 */
	
	/** Returns the local linear dispersion dN/dLambda at the given wavelength in m^-1.
	 * Default implementation does a numerical central derivative of N() at dl/l = 1e-5 */
	public double getLinearDispersion(int modeNumber, double wavelen, double temperature){
		double dLambda = wavelen * 1e-5;
		
		double Np = getRefractiveIndex(modeNumber, wavelen + dLambda/2, temperature);
		double Nm = getRefractiveIndex(modeNumber, wavelen - dLambda/2, temperature);
		return (Np - Nm) / dLambda;
		
	}
	
	/** Calc the Verdet constant (Faraday effect radians Tesla^-1 m^-1) from the Becquerel relations
	 * 
	 * This seems to be a common forumla that's valid for many, but not all, common glasses.
	 *  
	 *  See equation 1 of:
	 *  [ E.Munin "Faraday Effect And Energy Gap In Optical Matierals", J.Phys D (1992) ]
	 *  
	 *  From that paper, valid for:
	 *  	Fused Silica (Dynasil 1001)
	 *  	Schott BK-7 borosilicate glass 
	 *  	Undoped YAG crystal
	 *  
	 *  Also, a clear(er) derivataion from 
	 *  [ Kirk T. McDonald "Faraday Rotation", Princeton University (August 24, 2008)  
	 *	   http://www.physics.princeton.edu/~mcdonald/examples/faradayrotation.pdf (unpublished PDF)‎ ]
	 *
	 *  That states that the 
	 * 
	 */
	
	public static final double e = 1.60217657e-19;
	public static final double c = 299792458;
	public static final double m_e = 9.10938291e-31;
	
	public double verdetConstFromBecquerelRelation(int modeNumber, double wavelen, double temperature, 
								double magnetoOpticAnomaly) {

		double dNdLambda = getLinearDispersion(modeNumber, wavelen, temperature);
		
		//Equ1 in E.Munin1992 has wavelen in nm and dNdLambda in nm^-1
		// both at in m here, so no conversion is needed		
		//return -293.8 * magnetoOpticAnomaly * wavelen * dNdLambda;
		
		//from K.T.McDonald, he has c^2 in his final equation which has the wrong units, just e/mc earlier up
		// He says it's e/mc for 'spin', e/2mc for 'orbital' and half again in practice in many cases 
		// with e/2mc, it's the same as above (but provides the -ve sign that then agrees with normal glasses
		// (dN/dλ -ve), having +ve Verdet (clockwise rotation looking parallel to field)
		//
		return -e / (2 * m_e * c) * magnetoOpticAnomaly * wavelen * dNdLambda;
		//return -e / (m_e * c) * magnetoOpticAnomaly * wavelen * dNdLambda;
	}

	// Materials are required to provide hashCode() and equals()
	@Override
	public abstract int hashCode();

	@Override
	public abstract boolean equals(Object obj);
}
