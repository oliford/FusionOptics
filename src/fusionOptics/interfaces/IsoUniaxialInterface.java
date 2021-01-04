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
package fusionOptics.interfaces;

import java.text.DecimalFormat;

import fusionOptics.Util;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;

import binaryMatrixFile.BinaryMatrixWriter;

import algorithmrepository.exceptions.NotImplementedException;

import oneLiners.OneLiners;
import net.jafama.FastMath;


/**
 * Calculation of refracted rays and ray integrated amplitude transfer coefficients
 * for rays passing from an isotropic to a uniaxial material.
 * 
 * There is a nice overview of the situation here:
 * [ A. Weidlich and A.Wilkie "Realistic Rendering of Birefringency in Uniaxial Crystals"
 *   Institute of Computer Graphics and Algorithms, Vienna University of Technology ] 
 *   
 *  However, their maths is considerbly wrong in places and they really quite
 *  mess up their frame definitions. They took their maths mostly from two papers...
 *   
 *  The fresnel coefficients come from:  
 *     [ J.Lekner "Reflection and refraction by uniaxial crystals" 
 *   		J.Phys Condens Matter 3 (1991) 6121-6133 ] 

 * The extraordinary ray direction calculation comes from:
 *   [ M.Avendano-Alejo and O.N. Stavroudis "Huygens’s principle and rays 
 * 		in uniaxial anisotropic	media. II. Crystal axis orientation arbitrary"
 * 		J. Opt. Soc. Am. A Vol 19, No 8 (2002) ]
 * 
 * 
 * The errors in A.Weidlich's compsci paper are e.g. that they've
 * not noticed that there is both an 'a' and an 'alpha' in the final equs of 
 *  M.Avendano-Alejo's paper.
 *  
 *  The fresnel coefficients seem to make sense, and agree with figure 3 
 *   of A. Weidlich.
 *   
 * The ray direction code here now matches the test vectors given in
 * [ G.Beyerle "Ray tracing formulas for refraction and internal reflection in uniaxial crystals"
 *  Appl. Opt. 37, 7947-7953 (1998) ], which is used by M.Avendano 
 * 
 * @author oliford
 * 
 */
public class IsoUniaxialInterface extends DualMediumInterface {
	private double nonAbsorbedAmplitudeCoeff = 1.0;
	
	public boolean ignoreInterfaceCompatibility = false;
	
	private final static IsoUniaxialInterface ideal = new IsoUniaxialInterface();
	
	public static IsoUniaxialInterface ideal(){ return ideal; }
	
	/** If true, for rays that hit uniaxial glass perpendicular or parallel to the
	 * optic axis, the ray is not split and the birefringence is treated as a modification 
	 * to the polarisation states on the single ray as it propagates through the material.
	 * Rays that hit with a significant angle to the surface will throw an exception.
	 * 
	 * If false, the ray is split into separate E and O rays at the interface as it would be in the general case,
	 * even if the rays are actually parallel (i.e. for perp/para to optic axis).
	 * ** Currently, there are some problems with the phase shifts in this. **
	 */
	public static boolean simpleBirefringence = false;	
	
	/** Indices into fresnel coeffs arrays */
	private final static int F_Rss = 0;
	private final static int F_Rpp = 1;
	private final static int F_Rsp = 2;
	private final static int F_Rps = 3;
	private final static int F_Tso = 4;
	private final static int F_Tse = 5;
	private final static int F_Tpo = 6;
	private final static int F_Tpe = 7;
	
	@Override
	public void checkCompatibility(Surface surface) {
		int nAxesFront = surface.getFrontMedium() == null ? 0 : surface.getFrontMedium().getMaterial().getNAxes();
		int nAxesBack = surface.getBackMedium() == null ? 0 : surface.getBackMedium().getMaterial().getNAxes();
		
		if(ignoreInterfaceCompatibility || (nAxesFront == 0 && nAxesBack == 1) || (nAxesFront == 1 && nAxesBack == 0))
			return;
		
		throw new IllegalArgumentException("IsoUniaxialInterface and it's derivatives can only handle one isotropic medium and one uniaxial medium");
	}
	
	@Override
	/** Main entry from the raytracer. This really just splits between the
	 * uniaxial --> isotropic and isotropic --> uniaxial cases. */
	public void calcIntersection(Intersection hit, Medium incidentMedium,
			Medium transmissionMedium, double minIntensity) {
		
		if(simpleBirefringence){ //for the simple case, don't split the ray. Medium.rayPropagation() will deal with it.
			IsoIsoInterface.ideal().calcIntersection(hit, incidentMedium, transmissionMedium, minIntensity);
			return;
		}

		int nAxesIn = incidentMedium == null ? 0 : incidentMedium.getMaterial().getNAxes();
		int nAxesOut = transmissionMedium == null ? 0 : transmissionMedium.getMaterial().getNAxes();
		
		if(nAxesIn == 0 && nAxesOut == 1){
			isoToUni(hit, incidentMedium, transmissionMedium, minIntensity);
		}else if(nAxesIn == 1 && nAxesOut == 0){
			uniToIso(hit, incidentMedium, transmissionMedium, minIntensity);
		}else
			throw new RuntimeException("Can only handler isotropic-->uniaxial or uniaxial-->isotropic");
	}
	
	/** The standard, and more interesting, isotropic --> uniaxial case.
	 *  Here, we need to calculate:
	 *  	1) Ordinary and extraordinary transmitted ray directions.
	 *  	2) Extraordinary ray refractive index (which is a function of it's direction)
	 *  	3) Fresnel coefficents for reflected and both transmitted rays.
	 */
	public void isoToUni(Intersection hit, Medium incidentMedium,
			Medium transmissionMedium, double minIntensity) {
				
		//refractive indices
		double nI = incidentMedium == null ? 1.0 : incidentMedium.getRefractiveIndex(0, hit.incidentRay.wavelength);
		
		//refractive index of modes (ne is probably not the index for the transmitted extraordinary ray)
		double no = transmissionMedium.getRefractiveIndex(0, hit.incidentRay.wavelength);
		double ne = transmissionMedium.getRefractiveIndex(1, hit.incidentRay.wavelength);
		
		// setup our working frame, which has the incidence plane in (x,z) and surface in (x,y) 
		//Z points in the same direction to ray
		double z[] = new double[]{ -hit.normal[0], -hit.normal[1], -hit.normal[2] };
		double y[] = Util.cross(z, hit.incidentRay.dir);
		y = Util.dot(y,y) == 0 ? Util.createPerp(z) : Util.reNorm(y);
		double x[] = Util.cross(y, z);
			
		//Rotate the optic axis to that frame
		double mediumAxes[][] = transmissionMedium.getOpticAxes();
		double opticAxis[] = new double[]{
				Util.dot(mediumAxes[0], x),
				Util.dot(mediumAxes[0], y),
				Util.dot(mediumAxes[0], z),
		};
		
		// Incident angles 
		double cosThetaI = Util.dot(hit.incidentRay.dir, z);
		double sinThetaI = Util.dot(hit.incidentRay.dir, x);
	
		double rIndexRatio = nI / no;
	
		if( sinThetaI*sinThetaI >= 1/(rIndexRatio*rIndexRatio)){
			//total internal refraction
			//Reflector.pureReflection(hit, minIntensity, nonAbsorbedAmplitudeCoeff);
			throw new RuntimeException("Total internal reflection on entry to anisotropic media is not supported.");
		}
		
		//reflection direction
		double reflectedDir[] = new double[]{ sinThetaI, 0, -cosThetaI };
		
		//snells law for ordinary wave
		double sinThetaO = sinThetaI * nI / no;
		double ordinaryDir[] = (new double[]{
				sinThetaO, 
				0, 
				Math.sqrt(1 - sinThetaO*sinThetaO), 
			});
		
		//extra ordinary dir and fresnel coeffs from formulae in papers
		double ret[][] = calcTransmittedExtraordinaryDir(nI, no, ne, ordinaryDir, opticAxis);
		double extraordinaryDir[] = ret[0];
		double exWaveNormal[] = ret[1];
		double fresnel[][] = calcFresnelCoeffs(nI, no, ne, hit.incidentRay.wavelength, opticAxis, cosThetaI, sinThetaI);
		
		double cosThetaN = Util.dot(opticAxis, extraordinaryDir);
		//double cosThetaA = Util.dot(opticAxis, exWaveNormal);
		double sinSqThetaN = 1-cosThetaN*cosThetaN;
		//double sinSqThetaA = 1-cosThetaA*cosThetaA;
		
		//double nTENorm = no*ne / Math.sqrt(no*no*sinSqThetaA + ne*ne*cosThetaA*cosThetaA); //equ 31a
		double nTE = Math.sqrt(ne*ne*sinSqThetaN + no*no*cosThetaN*cosThetaN); //equ 31b
				
		//correct the E field amplitudes for the change of CSA by the refraction
		fresnel = correctForRayCrossSectionalArea(nI, no, nTE, cosThetaI, fresnel, ordinaryDir, extraordinaryDir);
		
		//calculate the direction of the 's' polaraistion (the incident component in the incidence plane)
		double sPolDir[] = Util.cross(hit.incidentRay.dir, hit.normal);
		if(Util.dot(sPolDir, sPolDir) == 0) 
			sPolDir = hit.incidentRay.up; //s and p are the same, keep the original
		else
			sPolDir = Util.reNorm(sPolDir); //otherwise, normalise it, so it's just a direction
		
		//rotate incident ray's polarisation definition to the s/p frame
		hit.incidentRay.rotatePolRefFrame(sPolDir);
		
		// Ordinary wave 
		hit.transmittedOrdinary = new RaySegment();
		hit.transmittedOrdinary.dir = reFrame(x, y, z, ordinaryDir);	
		hit.transmittedOrdinary.up = reFrame(x, y, z, fresnel[1]);
		hit.transmittedOrdinary.E0 = Pol.complexMulAll(hit.incidentRay.E1, 
				fresnel[0][F_Tso], 0, 0, 0, fresnel[0][F_Tpo], 0, 0, 0, nonAbsorbedAmplitudeCoeff);
		hit.transmittedOrdinary = completeRay(hit.transmittedOrdinary, hit, transmissionMedium, minIntensity);
		
		// Extraordinary wave	
		hit.transmittedExtraordinary = new RaySegment();
		hit.transmittedExtraordinary.dir = reFrame(x, y, z, extraordinaryDir);	
		hit.transmittedExtraordinary.up = reFrame(x, y, z, fresnel[2]);
		hit.transmittedExtraordinary.E0 = Pol.complexMulAll(hit.incidentRay.E1, 
				fresnel[0][F_Tse], 0, 0, 0, fresnel[0][F_Tpe], 0, 0, 0, nonAbsorbedAmplitudeCoeff);
		hit.transmittedExtraordinary.raySpecificRefractiveIndex = nTE; //and we need to store the refractive index for that one
		hit.transmittedExtraordinary = completeRay(hit.transmittedExtraordinary, hit, transmissionMedium, minIntensity);
		
		// Reflected wave
		hit.reflectedOrdinary = new RaySegment();
		hit.reflectedOrdinary.dir = reFrame(x, y, z, reflectedDir);	
		hit.reflectedOrdinary.up = sPolDir.clone();
		hit.reflectedOrdinary.E0 = Pol.complexMulAll(hit.incidentRay.E1, 
				fresnel[0][F_Rss], 0, fresnel[0][F_Rsp], 0, fresnel[0][F_Rps], 0, fresnel[0][F_Rpp], 0, nonAbsorbedAmplitudeCoeff);
		hit.reflectedOrdinary = completeRay(hit.reflectedOrdinary, hit, incidentMedium, minIntensity);		
	}
	
	/** Uniaxial --> isotropic, 'inverse' case.
	 * 
	 * For an exiting ordinary ray, this is simply snells law.
	 * For an exiting extraordinary ray, as in M.Avendano, we look 
	 *   at this as the reverse of this extraordinary ray being creating
	 *   by an incoming ray on this side, and solve first for the fictitious 
	 *   ordinary ray that doesn't really exist. The exiting ray is again
	 *   just snells law applied to that.			
	 *   
	 *   */
	public void uniToIso(Intersection hit, Medium incidentMedium,
			Medium transmissionMedium, double minIntensity) {
		
		//refractive indices, nX = exit, no = ordinary, ne = extraordinary
		double nX = transmissionMedium == null ? 1.0 : transmissionMedium.getRefractiveIndex(0, hit.incidentRay.wavelength);
		
		double no = incidentMedium.getRefractiveIndex(0, hit.incidentRay.wavelength);
		double ne = incidentMedium.getRefractiveIndex(1, hit.incidentRay.wavelength);
		
		double nI;
		
		/** The (possibly fictitious) ordinary ray to apply snells law to
		 * find the exiting ray */
		double ordinaryWorldDir[];
		
		//If the stored refractive index is zero or invalid, it  
		// means it was the ordinary ray.
		if(Double.isNaN(hit.incidentRay.raySpecificRefractiveIndex)  
				|| hit.incidentRay.raySpecificRefractiveIndex == 0
				|| hit.incidentRay.raySpecificRefractiveIndex == no){
		
			//in which case, we can just do snells law to find the existing ray
			ordinaryWorldDir = hit.incidentRay.dir;
			nI = no;
			
		}else{ //Extraordinary ray, we need to find the fictitious ordinary ray.
			 
				
			// setup our working frame, which has the extraordinary ray in the plane (x,z) and material surface in (x,y). 
			//Z points in opposite direction to ray. This isn't quite the inverse of the problem, since in the original
			// the incident ray and ordinary ray were in the (x,z) plane, but I think the math still works.
			double extraordinaryDir[] = new double[]{
					-hit.incidentRay.dir[0], //the extraordinary dir is actually opposite,					
					-hit.incidentRay.dir[1], //since we're almost inverting the whole problem
					-hit.incidentRay.dir[2],
				};
			
			double z[] = new double[]{ hit.normal[0], hit.normal[1], hit.normal[2] };
			double y[] = Util.cross(z, extraordinaryDir);
			
			if(Util.dot(y,y) == 0){ // extraordinary dir was (a)parallel with normal, we can just invent the y direction
				y = Util.createPerp(z);
			}else{
				y = Util.reNorm(y);
			}
			double x[] = Util.cross(y, z);
				
			//Rotate the optic axis to that frame
			double mediumAxes[][] = incidentMedium.getOpticAxes();
			double opticAxis[] = new double[]{
				Util.dot(mediumAxes[0], x),
				Util.dot(mediumAxes[0], y),
				Util.dot(mediumAxes[0], z),
			};
			
			//Rotate extraordinary dir to the frame
			extraordinaryDir = new double[]{
				Util.dot(extraordinaryDir, x),
				0,
				Util.dot(extraordinaryDir, z),				
			};
			
			//calculate the fictitious ordinary ray
			double ordinaryDir[] = calcExtraordinaryExitDir(nX, no, ne, extraordinaryDir, opticAxis);
			
			//and convert that back to world coordinartes
			ordinaryWorldDir = new double[]{
					-ordinaryDir[0] * x[0] - ordinaryDir[1] * y[0] - ordinaryDir[2] * z[0],
					-ordinaryDir[0] * x[1] - ordinaryDir[1] * y[1] - ordinaryDir[2] * z[1],
					-ordinaryDir[0] * x[2] - ordinaryDir[1] * y[2] - ordinaryDir[2] * z[2],
			};
			
			nI = hit.incidentRay.raySpecificRefractiveIndex;
		}
		
		//now we just need to refract the ordinary ray as normal.
		double rIndexRatio = no / nX;

		double cosThetaI = -Util.dot(ordinaryWorldDir, hit.normal);
		
		if( (1 - cosThetaI*cosThetaI) >= 1/(rIndexRatio*rIndexRatio)){
			//total internal refraction
			//Reflector.pureReflection(hit, minIntensity, nonAbsorbedAmplitudeCoeff);
			//return;
			System.err.println("Total internal reflection on exit from anisotropic media is not supported. Aborting ray");
			return;
			
		}

		double a = + rIndexRatio * cosThetaI;
		double cosThetaX = Math.sqrt( 1 - rIndexRatio*rIndexRatio*(1-cosThetaI*cosThetaI));
		double c = a - cosThetaX;
		
		double exitDir[] = Util.reNorm(new double[]{
				c * hit.normal[0] + rIndexRatio * ordinaryWorldDir[0], 
				c * hit.normal[1] + rIndexRatio * ordinaryWorldDir[1], 
				c * hit.normal[2] + rIndexRatio * ordinaryWorldDir[2], 
			});
		
		//calculate the direction of the 's' polarisation (the incident component in the incidence plane)
		double sPolDir[] = Util.cross(hit.incidentRay.dir, hit.normal);
		if(Util.dot(sPolDir, sPolDir) == 0) 
			sPolDir = hit.incidentRay.up; //s and p are the same, keep the original
		else
			sPolDir = Util.reNorm(sPolDir); //otherwise, normalise it, so it's just a direction
			
		//rotate incident ray's polarisation definition to the s/p frame
		hit.incidentRay.rotatePolRefFrame(sPolDir);
		
		//FIXME: I couldn't find a nice paper giving the Fresnel coefficients
		// for this case. In principal I guess I could just derive it using the boundary
		//conditions as in J.Lekner, but I've no idea what happens to the internally reflected
		// extraordinary component, since you need to find the new ray dir etc for it.
		
		// For now, we are just assuming it is 100% transmitted for the 
		// uniaxial -->  isotropic case since people will probably
		// be more interested in the entry case.
		
		// Ordinary wave 
		hit.transmittedOrdinary = new RaySegment();
		hit.transmittedOrdinary.dir = exitDir;	
		hit.transmittedOrdinary.up = sPolDir;
		//There are no electric field amplitude coeffs, so we don't need to do the deltaCSA compensation.
		hit.transmittedOrdinary.E0 = Pol.complexMulAll(hit.incidentRay.E1, 
				1, 0, 0, 0, 0, 0, 1, 0, nonAbsorbedAmplitudeCoeff);
		hit.transmittedOrdinary = completeRay(hit.transmittedOrdinary, hit, transmissionMedium, minIntensity);
	}
	
	/** Fresnel coefficients from J.Lekner 
	 * OpticAxis should already be in frame with incidence plane as (x,z) and surface as (x,y).
	 * 
	 * ** This is reproduction of maths in the paper, and is therefore not very readable. **
	 * 
	 * @param ni	Refractive index of incident medium 
	 * 		(I'm making the blind assumption that this approximates to the refractive index of the incident 
	 * 			ray in the incident medium, if the incident medium is ansiotropic and the ray extraordinary)
	 * @param no	Refractive index of ordinary ray in transmitting medium.	
	 * @param ne	Refractive index of (completely)extraordinary ray in transmitting medium, i.e. the 
	 * 						optic axis parallel polarisation when prop. dir. is perp to optic axis.
	 * @param wavelen	Wavelength of the incident ray.
	 * @param opticAxis		Optic axis vector w.r.t surface coord system.
	 * @param cosThetaI		Cosine of angle of incidence.
	 */
	private static double[][] calcFresnelCoeffs(double ni, double no, double ne, double wavelen, double opticAxis[], double cosThetaI, double sinThetaI) {
		
		//dielectric constants
		double epO = no * no;
		double epE = ne * ne;
		double dEp = epE - epO;
		
		//incidence angle trig
		double tanThetaI = sinThetaI / cosThetaI;
		
		double k = 2 * Math.PI / wavelen;	//wave normal number
		double K = k * ni * sinThetaI; 		// 'transverse wave vector for all waves'
		
		//cosines of optic axis
		double alpha = opticAxis[0];
		double beta = opticAxis[1];
		double gamma = opticAxis[2];
		
		double q1 = k*ni*cosThetaI;
		double qo = FastMath.sqrt(epO*k*k - K*K);
		double d = epO*(epE*(epO + gamma*gamma*dEp)*k*k - (epE - beta*beta*dEp)*K*K);
		double qe = (FastMath.sqrt(d) - alpha*gamma*K*dEp) / (epO + gamma*gamma*dEp);
		
		double Eox,Eoy,Eoz, Eex,Eey,Eez;
		
		
		if(cosThetaI == 1.0 && (alpha*alpha + beta*beta) <= 0){
			Eox = 1; Eoy = 0; Eoz = 0;
			Eex = 0; Eey = 1; Eez = 0;
		}else{
			//Electric field vector solutions in the medium
			Eox = -beta*qo;
			Eoy = alpha*qo - gamma*K;
			Eoz = beta*K;
			double Eol = FastMath.sqrt(Eox*Eox + Eoy*Eoy + Eoz*Eoz);
			Eox /= Eol; Eoy /= Eol; Eoz /= Eol;
			
			Eex = alpha*qo*qo - gamma*qe*K;
			Eey = beta*epO*k*k;
			Eez = gamma*(epO*k*k - qe*qe) - alpha*qe*K;
			double Eel = FastMath.sqrt(Eex*Eex + Eey*Eey + Eez*Eez);
			Eex /= Eel; Eey /= Eel; Eez /= Eel;
		}
		
		double A = (qo + q1 + K*tanThetaI)*Eox - K*Eoz;
		double B = (qe + q1 + K*tanThetaI)*Eex - K*Eez;
		double D = (q1 + qe)*A*Eey - (q1 + qo)*B*Eoy;
		
		double qt = q1 + K*tanThetaI;
		
		//reflection coeffs, s to s, p to p, etc
		double rss = ((q1 - qe)*A*Eey - (q1 - qo)*B*Eoy) / D;
		double rpp = 2*qt*((q1 + qe)*Eox*Eey - (q1 + qo)*Eex*Eoy)/D - 1.0;
		double rsp = 2*ni*k*(A*Eex - B*Eox) / D;
		double rps = 2*ni*k*(qe - qo)*Eoy*Eey / D;
		
		//transmission coeffs: s to ordinary, s to extraordinary etc
		// these are changes to the E field magnitude everywhere inside the ray
		// so do not include the compensation for the change of ray 
		// cross-sectional area!
		
		// These signs here were -,-,+,- from the paper but they've been
		//fiddled here to make sure the overall polarisation in the tracer
		// doesn't change on entry to a uniaxial medium, which I don't think
		// they should. ( See BirefingentRaySplittingTest.testEntryPhases() )
		// To pass that test they can be -+-+ or +-+-
		
		// They also effect the phase result of TestBirefringentRayMaths
		// but adding multiples of 180° to the overall phase difference.
		double tso = -2*q1*B / D;
		double tse = -2*q1*A / D;
		double tpo = 2*ni*k*(q1 + qe)*Eey / D;
		double tpe = -2*ni*k*(q1 + qo)*Eoy / D;
		
		return new double[][]{ 
				{ rss, rpp, rsp, rps,    tso, tse, tpo, tpe, },
				{ Eox, Eoy, Eoz }, 
				{ Eex, Eey, Eez } 
			};
	}
	
	
	/**E fields in the transmitted beam are weaker because it would be wider than the incident/reflected
	 *	Since we want to model them as per unit area, we need to correct for that change in cross-sectional area
	 * Also, the power passing the interface in 1 unit time is held in a smaller volume in the medium
	 * in which the light is travelling slower, so it has a higher energy density and hence a higher E field.
	 *
	 * Originally mine, also mentioned in A.Weidlich.
	 * 
	 * @param ni			Incident refractive index
	 * @param no			Ordinary transmitted refractive index (same as ordinary mode index of medium) 
	 * @param nTE			Refractive index of extraordinary ray. Not the same as ne - the extraordinary mode in the medium.
	 * @param cosThetaI		Incident angle.
	 * @param fresnel		Fresnel coefficients for Electric field
	 * @param dirs			Directions in medium of transmitted rays [O,E][x,y,z]
	 * @return				Fresnel coefficients for ray integrated electric field.
	 */
	private static double[][] correctForRayCrossSectionalArea(double ni, double no, double nTE, double cosThetaI, double[][] fresnel,  double[] ordinaryDir, double[] extraordinaryDir) {
		double cosThetaTO = ordinaryDir[2];
		double cosThetaTE = extraordinaryDir[2];
		
		
		double deltaCSAO = FastMath.sqrt((no * cosThetaTO) / (ni * cosThetaI));
		double deltaCSAE = FastMath.sqrt((nTE * cosThetaTE) / (ni * cosThetaI));
		
		fresnel[0][F_Tso] *= deltaCSAO; // tso
		fresnel[0][F_Tse] *= deltaCSAE; // tse
		fresnel[0][F_Tpo] *= deltaCSAO; // tpo
		fresnel[0][F_Tpe] *= deltaCSAE; // tpe
		
		return fresnel;
	}

	/** Calculate the extraordinary ray direction given
	 * the ordinary ray direction.
	 * 
	 *  From M.Avendano. Matches test cases perfectly.
	 *  Assumes surface in (x,y) and incidence in (x,z). 
	 *  
	 *  ** This is reproduction of maths in the paper, and is therefore not very readable. **
	 * 
	 *  This is a bit like in A.Weidlich, but that is wrong. 
	 * */
	public static double[][] calcTransmittedExtraordinaryDir(double ni, double no, double ne, double ordinaryDir[], double opticAxis[]){

		//cosines of optic axis
		double alpha = opticAxis[0];
		double beta = opticAxis[1];
		double gamma = opticAxis[2];
		
		double N = no*no - ne*ne;
		double Gamma = ne*ne + N*(1.0 - gamma*gamma);
		double deltaOSq = ne*ne - no*no*(1 - ordinaryDir[2]*ordinaryDir[2]);
		double a = alpha*ordinaryDir[0] + beta * ordinaryDir[1];
		double Delta = FastMath.sqrt(Gamma*deltaOSq + no*no*N*a*a);
		double deltaEG = FastMath.sqrt(ne*ne*(no*no*Gamma*Gamma - N*(gamma*Delta + no*no*a)*(gamma*Delta + no*no*a)));
		
		double XE = (no*no*(ne*ne + N*beta*beta) * ordinaryDir[0] -no*no*N*alpha*beta * ordinaryDir[1] - gamma*N*Delta*alpha) / deltaEG;
		double YE = (-no*no*N*alpha*beta*ordinaryDir[0] + no*no*(ne*ne + N*alpha*alpha)*ordinaryDir[1] - gamma*N*Delta*beta) / deltaEG;
		double ZE = (Gamma*Delta - gamma*N*Delta*0) / deltaEG;
		double lE = FastMath.sqrt(XE*XE + YE*YE + ZE*ZE);
		XE /= lE; YE /= lE; ZE /= lE;
		
		//while we have the maths, also calc the extraordinary wave normal vector
		
		double Ngx = (ordinaryDir[0]*(N*a*gamma - Delta));
		double Ngy = (ordinaryDir[1]*(N*a*gamma - Delta));
		double Ngz = (ordinaryDir[2]*-(deltaOSq + N*a*a));
		double lN = FastMath.sqrt((1 - ordinaryDir[2]*ordinaryDir[2])*(N*a*gamma - Delta)*(N*a*gamma - Delta) + (deltaOSq + N*a*a)*(deltaOSq + N*a*a));
		//double lN = FastMath.sqrt(Ngx*Ngx + Ngy*Ngy + Ngz*Ngz);
		Ngx /= lN; Ngy /= lN; Ngz /= lN;
		
		//incidence is in (x,z) +ve in z, so the wave normal should also be +ve in z
		if(Ngz < 0){
			Ngx = -Ngx;
			Ngy = -Ngy;
			Ngz = -Ngz;
		}
		
		return new double[][]{  {XE, YE, ZE }, { Ngx, Ngy, Ngz } };
	}
	
	/** Calculates the ordinary ray direction inside the medium, given the extraordinary
	 * ray direction (also inside). From M.Avendano. 
	 * 
	 * Assumes surface in (x,y) and incidence in (x,z). 
	 * 	
	 * ** This is reproduction of maths in the paper, and is therefore not very readable. **
	 */ 
	public static double[] calcExtraordinaryExitDir(double ni, double no, double ne, double exDir[], double opticAxis[]){

		//cosines of optic axis
		double alpha = opticAxis[0];
		double beta = opticAxis[1];
		double gamma = opticAxis[2];
		
		double N = no*no - ne*ne;
		
		double mu = FastMath.sqrt(ne*ne + N * FastMath.pow2(alpha*exDir[0] + beta*exDir[1] + gamma*exDir[2]));
		
		double SeDotA = Util.dot(exDir, opticAxis);
		double DeltaE = FastMath.sqrt(
				ne*ne*(N + ne*ne*exDir[2]*exDir[2]) +
				N*(SeDotA) * ( (N*gamma*gamma - ne*ne)*(SeDotA) + 2*ne*ne*gamma*exDir[2])
				);
		
		double XO = ((ne*ne + N*alpha*alpha)*exDir[0] + N*alpha*beta*exDir[1] + exDir[2]*gamma*N*alpha) / (no*mu);
		double YO = ((N*alpha*beta)*exDir[0] + (ne*ne + N*beta*beta)*exDir[1] + exDir[2]*gamma*N*beta) / (no*mu);
		double ZO = DeltaE / (no*mu);
		
		return new double[]{  XO, YO, ZO  };
	}


	/** Create a world coord vector by projecting v onto x,y,z **/
	public static final double[] reFrame(double x[], double y[], double z[], double v[]){
		return new double[]{
				v[0] * x[0] + v[1] * y[0] + v[2] * z[0], 
				v[0] * x[1] + v[1] * y[1] + v[2] * z[1], 
				v[0] * x[2] + v[1] * y[2] + v[2] * z[2], 
			};
	}
	
	/** Fill in the common ray info, and delete the ray if the intensity is too low */
	private RaySegment completeRay(RaySegment ray, Intersection hit, Medium medium, double minIntensity){
		if(ray.startIntensity() < minIntensity){
			return null;
		} else {
			ray.startHit = hit;
			ray.startPos = hit.pos;
			
			//set-up both rays for tracing on
			ray.length = Double.POSITIVE_INFINITY;
			ray.medium = medium;
			ray.wavelength = hit.incidentRay.wavelength;
			ray.endHit = null;
			return ray;
		}
	}

	@Override
	public int hashCode() {
		long t; int r = 1;
		t = Double.doubleToLongBits(nonAbsorbedAmplitudeCoeff); r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		IsoUniaxialInterface other = (IsoUniaxialInterface) obj;
		return Double.doubleToLongBits(nonAbsorbedAmplitudeCoeff) == Double.doubleToLongBits(other.nonAbsorbedAmplitudeCoeff);
			
	}
}
