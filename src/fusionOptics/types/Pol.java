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

import java.text.DecimalFormat;
import java.util.concurrent.ConcurrentLinkedQueue;

import fusionOptics.Util;
import fusionOptics.surfaces.Plane;
import fusionOptics.surfaces.Square;

import otherSupport.Profiler;

import net.jafama.FastMath;


/** Describes a state of polarisation and provides various 
 * simple conversion and maths routine, avoiding trigonometry 
 * and handling infinities where possible.
 * 
 * Each individual polarisation is treated as a completely coherent
 * single wave through the system, since the ray tracing is a
 * 'moment in time' thing. 
 * 
 * The concept of 'unpolarised light' is really a short-hand to say
 * that a one moment in time, the initial waves may destructively interfere
 * while at the next they might constructively interfere etc.
 * Of course, rays that are split and recombined will always
 * be coherent w.r.t other paths of themselves, regardless of
 * the initial phases.
 *
 * So, in order to treat the behaviour of 'unpolarised light', simply
 * start with every ray with two orthogonal polarisations.
 * At each point, plot the sum of those two, with varying phases 
 * between them. You don't need to re-trace to recalculate the polarisation
 * with different initial phases, because the phases of each Polarisation
 * object are relative to the corresponding one at the start of the 
 * first ray segment. 
 * 
 * The PSF subsystem has an optimise way of dealing with incoherent light,
 * by adding 4 states of polaristion to each starting ray, and from the final
 * state of each of them, solving for the Muller matrix of the whole system.
 * 
 * @author oliford
 */
public abstract class Pol {
	public final static int uRe = 0;
	public final static int uIm = 1;
	public final static int rRe = 2;
	public final static int rIm = 3;
	
	private static final DecimalFormat strFmt = new DecimalFormat("#.##");
	
	private static double pols[][][] = null;
	private static int nextPol = 0;
	
	/** Creates a new polarisation object by multiplying by ghe given complex coefficients 
	 * 
	 * @param CuRe	Up multiplier, real part
	 * @param CuIm	Up multiplier, imaginary part
	 * @param CrRe	Right multiplier, real part
	 * @param CrIm	Right multiplier, imaginary part
	 * @param commonMul common multiplier for reduction of whole magnitude
	 */
	public static final double[] complexMul(double E[], double CuRe, double CuIm, double CrRe, double CrIm, double commonMul){
		return new double[]{
			commonMul * (CuRe*E[uRe] - CuIm*E[uIm]), 
			commonMul * (CuRe*E[uIm] + CuIm*E[uRe]),
			commonMul * (CrRe*E[rRe] - CrIm*E[rIm]), 
			commonMul * (CrRe*E[rIm] + CrIm*E[rRe]) 
		};
	}

	/** Creates a new polarisation object by multiplying by the given complex coefficients 
	 * 
	 * @param CuRe	Up multiplier, real part
	 * @param CuIm	Up multiplier, imaginary part
	 * @param CrRe	Right multiplier, real part
	 * @param CrIm	Right multiplier, imaginary part
	 * @param commonMul common multiplier for reduction of whole magnitude
	 */
	public static long t0 = System.nanoTime();;
	public static long tNowt;
	public static long tAlloc;
	public static long tMul;
	public static long tOut;
	
	public static final double[][] complexMulAll(double E[][], double CuRe, double CuIm, double CrRe, double CrIm, double commonMul){
		//double retE[][] = new double[E.length][4];
		double retE[][] = alloc(E.length);
		for(int i=0; i < E.length; i++){
			retE[i][uRe] =  commonMul * (CuRe*E[i][uRe] - CuIm*E[i][uIm]); 
			retE[i][uIm] =  commonMul * (CuRe*E[i][uIm] + CuIm*E[i][uRe]);
			retE[i][rRe] =  commonMul * (CrRe*E[i][rRe] - CrIm*E[i][rIm]); 
			retE[i][rIm] =  commonMul * (CrRe*E[i][rIm] + CrIm*E[i][rRe]); 
		}
		return retE;
	}
	
	public static final double[][] complexMulAll(double E[][], 
			double CuuRe, double CuuIm, double CurRe, double CurIm, 
			double CruRe, double CruIm, double CrrRe, double CrrIm, 
			double commonMul){
		//double retE[][] = new double[E.length][4];
		double retE[][] = alloc(E.length);
		for(int i=0; i < E.length; i++){
			retE[i][uRe] =  commonMul * (CuuRe*E[i][uRe] - CuuIm*E[i][uIm]   + CruRe*E[i][rRe] - CruIm*E[i][rIm]); 
			retE[i][uIm] =  commonMul * (CuuRe*E[i][uIm] + CuuIm*E[i][uRe]   + CruRe*E[i][rIm] + CruIm*E[i][rRe]);
			retE[i][rRe] =  commonMul * (CrrRe*E[i][rRe] - CrrIm*E[i][rIm]   + CurRe*E[i][uRe] - CurIm*E[i][uIm]); 
			retE[i][rIm] =  commonMul * (CrrRe*E[i][rIm] + CrrIm*E[i][rRe]   + CurRe*E[i][uIm] + CurIm*E[i][uRe]); 
		}
		return retE;
	}
	
	/** Allocation of pols arrays, since they are used so much
	 * we don't want the GC to continually allocate and deallocate them.
	 * This actually only helps a little.
	 * 
	 * @param l
	 * @return
	 */
	public final static double[][] alloc(int l){
		if(pols == null){
			pols = new double[10][l][4];
			nextPol = 0;
			//System.out.println("Pols preallocing for length " + l + " x100");
		}else if(pols[0].length != l){
			pols = new double[pols.length][l][4];
			nextPol = 0;
			//System.out.println("Pols reallocing for change of length to " + l + " x"+pols.length);
		}else if(nextPol >= pols.length){
			pols = new double[pols.length*2][l][4];
			nextPol = 0;
			//System.out.println("Pols reallocing because too short length " + l + " x"+pols.length);
		}
		return pols[nextPol++];
	}
	
	/** Deallocates all the polarisations arrays.
	 * TODO: Find out if all this is really necessary.  */
	public final static void recoverAll() {
		nextPol = 0;
	}
	
	public static final double[][] scaleAll(double[][] E, double scale) {
		//double retE[][] = new double[E.length][4];
		double retE[][] = alloc(E.length);
		for(int i=0; i < E.length; i++) {
			retE[i][uRe] = E[i][uRe] * scale; 
			retE[i][uIm] = E[i][uIm] * scale;
			retE[i][rRe] = E[i][rRe] * scale; 
			retE[i][rIm] = E[i][rIm] * scale; 
		}
		return retE;
	}

	public static final double[][] copyAll(double[][] E) {
		return scaleAll(E, 1);
	}

	/** Some polarisation helpers with trig avoidance and NaN mitigation */
	public static final double intensity(double E[]) { return E[uRe]*E[uRe] + E[uIm]*E[uIm] + E[rRe]*E[rRe] + E[rIm]*E[rIm] /*2*Enopol*Enopol*/;	}

	public static final double intensity(double[][] E) {
		double sumI = 0;
		for(int j=0; j < E.length; j++)
			sumI += intensity(E[j]);
		return sumI;
	}
	
	public static final double mag2Eu(double E[]){ return E[uRe]*E[uRe] + E[uIm]*E[uIm]; }
	public static final double mag2Er(double E[]){ return E[rRe]*E[rRe] + E[rIm]*E[rIm]; }

	public static final double tanPhaseEu(double E[]){ return (E[uRe] == 0) ? Double.POSITIVE_INFINITY : (E[uIm] / E[uRe]); }
	public static final double tanPhaseEr(double E[]){ return (E[rRe] == 0) ? Double.POSITIVE_INFINITY : (E[rIm] / E[rRe]); }
	
	public static final double cosPhaseEu(double E[]){ double mag = FastMath.sqrt(E[uRe]*E[uRe] + E[uIm]*E[uIm]); return (mag == 0) ? 0 : (E[uRe] / mag); }
	public static final double sinPhaseEu(double E[]){ double mag = FastMath.sqrt(E[uRe]*E[uRe] + E[uIm]*E[uIm]); return (mag == 0) ? 0 : (E[uIm] / mag); }
	
	public static final double cosPhaseEr(double E[]){ double mag = FastMath.sqrt(E[rRe]*E[rRe] + E[rIm]*E[rIm]); return (mag == 0) ? 0 : (E[rRe] / mag); }
	public static final double sinPhaseEr(double E[]){ double mag = FastMath.sqrt(E[rRe]*E[rRe] + E[rIm]*E[rIm]); return (mag == 0) ? 0 : (E[rIm] / mag); }
	
	public static final double polarisedIntensityFrac(double E[]){ return (mag2Eu(E) + mag2Er(E)) / intensity(E); }
	
	/** tan(phaseR - phaseU) */
	public static final double tanPhaseDiff(double E[]){
		if(E[uRe] == 0){
			return (E[rIm] == 0) ? Double.POSITIVE_INFINITY : (-E[rRe] / E[rIm]);			
		}else if(E[rRe] == 0){
			return (E[uIm] == 0) ? Double.POSITIVE_INFINITY : (E[uRe] / E[uIm]);
		}else{
			double tanPhsU = E[uIm] / E[uRe];
			double tanPhsR = E[rIm] / E[rRe];
			return (tanPhsR - tanPhsU) / (1 + tanPhsR*tanPhsU);
		}
		
	}

	/** cos(phaseR - phaseU) */
	public static final double cosPhaseDiff(double E[]){
		double magEu = FastMath.sqrt(E[uRe]*E[uRe] + E[uIm]*E[uIm]);
		double magEr = FastMath.sqrt(E[rRe]*E[rRe] + E[rIm]*E[rIm]);
		return (magEu==0 || magEr==0) ? 0 : ((E[rRe]*E[uRe] + E[rIm]*E[uIm]) / (magEu*magEr));
	}
	
	/** sin(phaseR - phaseU) */
	public static final double sinPhaseDiff(double E[]){
		double magEu = FastMath.sqrt(E[uRe]*E[uRe] + E[uIm]*E[uIm]);
		double magEr = FastMath.sqrt(E[rRe]*E[rRe] + E[rIm]*E[rIm]);
		return (magEu==0 || magEr==0) ? 0 : ((E[rIm]*E[uRe] - E[rRe]*E[uIm]) / (magEu*magEr));
	}
	
	
	/** 4 component stokes vector describing the polarisation state of the ray, at the start.
	 * s[0] = I 
	 * s[1] = I * p * cos(2.chi) * cos(2.psi)
	 * s[2] = I * p * cos(2.chi) * sin(2.psi)  
	 * s[3] = I * p * sin(2.chi)
	 * 
	 * chi is the ellipticity angle
	 * psi is the angle of the ellipse major axis away from 'up', +ve for clockwise in the perp plane when viewed 
	 * looking along the dir direction, so that chi=0,psi=+90 deg is polarised parallel to the 'right' direction.   
	 * 
	 * I is the total intensity and p the polarisation fraction
	 */
	public static final double s0(double E[]){ return (E[uRe]*E[uRe] + E[uIm]*E[uIm]) + (E[rRe]*E[rRe] + E[rIm]*E[rIm]); } //same as intensity
	public static final double s1(double E[]){ return (E[uRe]*E[uRe] + E[uIm]*E[uIm]) - (E[rRe]*E[rRe] + E[rIm]*E[rIm]); }
	public static final double s2(double E[]){ return 2 * (E[uRe]*E[rRe] + E[uIm]*E[rIm]); }
	public static final double s3(double E[]){ return 2 * (E[uIm]*E[rRe] - E[uRe]*E[rIm]); }
	
	public static final double sin2Chi(double E[]) {
		double s1=s1(E), s2=s2(E), s3=s3(E);
		double pI2 = s1*s1 + s2*s2 + s3*s3;
		return (pI2 == 0) ? 0 : (s3 / FastMath.sqrt(pI2));
	}
	
	public static final double cos2Chi(double E[]) {
		double s1=s1(E), s2=s2(E), s3=s3(E);
		double pI2 = s1*s1 + s2*s2 + s3*s3;
		return (pI2 == 0) ? 0 : FastMath.sqrt((s1*s1 + s2*s2) / pI2);
	}
	
	public static final double tan2Chi(double E[]) {
		double s1=s1(E), s2=s2(E), s3=s3(E);
		double pIcos2ChiSq = s1*s1 + s2*s2;
		return (pIcos2ChiSq == 0) ? Double.POSITIVE_INFINITY : (s3 / FastMath.sqrt(pIcos2ChiSq));
	}


	/** Slow helpers */
	public static final double chi(double E[]) {
		double s1=s1(E), s2=s2(E), s3=s3(E);		
		return FastMath.atan2(s3, FastMath.sqrt(s1*s1 + s2*s2)) / 2;
	}

	public static final double psi(double E[]) {
		double s1=s1(E), s2=s2(E);		
		return FastMath.atan2(s2, s1) / 2;
	}

	/** Modifies the electic field amplitudes so that the actual polarisation
	 * is maintained, if the 'up' vector has been changed.
	 *  
	 * @param oUnU	(oldUp . newUp) == (oldRight . newRight)
	 * @param oUnR  (oldUp . newRight) == -(oldRight . newUp)
	 */
	public static final double[] rotateFrame(double E[], double oUnU, double oUnR) {
		return new double[]{ 
			oUnU * E[uRe] - oUnR * E[rRe],
			oUnU * E[uIm] - oUnR * E[rIm],
			oUnR * E[uRe] + oUnU * E[rRe],
			oUnR * E[uIm] + oUnU * E[rIm]
		};
	}
	
	
	
	public static final String toString(double E[]) {
		return "I="+strFmt.format(intensity(E))+
				",ψ="+strFmt.format(psi(E)*180/Math.PI) +
				"°,χ="+strFmt.format(chi(E)*180/Math.PI) +
				"°,φu="+strFmt.format(FastMath.atan2(E[Pol.uIm],E[Pol.uRe])*180/Math.PI) +
				"°,φr="+strFmt.format(FastMath.atan2(E[Pol.rIm],E[Pol.rRe])*180/Math.PI) + "°";
	}

	
	/** Projects the final polarisation of incident ray of given intersection, into the up/right frame
	 * of that plane. The returned frame is looking along the PLANE's normal (NOT the intersection normal)
	 * , with up aligned to it's up.
	 * 
	 * NB: This IS sensible, if you want your final polarisation to look like the incident light, simply point the
	 * final surface normal in the same direction (so the light hits the non-normal side)
	 *  
	 * @param reduceP If false, the full amplitude of the P polarisation is transferred, if true, it is really
	 * 					projected, i.e. only component of E parallel to both the surface and incidence planes
	 * 					is transferred. Setting true (project) makes some kind of geometric sense but setting
	 * 					false (full transfer) seems to be more physical since it replicates the effect of putting
	 * 					a perfect collimating lens just before the surface.
	 */
	public static double[][] projectToPlanesView(Intersection hit, boolean reduceP) {
		Plane plane = (Plane)hit.surface;

		//normal of plane 
		double n[] = plane.getNormal();
		
		double a = Util.dot(hit.incidentRay.dir, n);
		
		//polarisation direction in plane surface, perp to incidence plane
		double sPolDir[] = Util.cross(hit.incidentRay.dir, n);
		double sLen = Util.dot(sPolDir, sPolDir);
		if(sLen == 0){ //if completely perp, just use plane up
			sPolDir = plane.getUp();
		}
		
		Util.reNorm(sPolDir);
		
		hit.incidentRay.rotatePolRefFrame(sPolDir);
		
		
		//s pol (now 'up') is already in plane, p pol (now 'right') needs to be projected, or not, I'm not sure
		if(!reduceP)
			a /= Math.abs(a); //we still need the sign though
		double EInPlane[][] = Pol.complexMulAll(hit.incidentRay.E1, 1, 0, a, 0, 1.0);
		
		//now we need to rotate the projected frame so that 'up' is the plane's up
		//we need our 'right' for this, which is not actually the same as Plane.getRight(): (it's -ve)
		double right[] = Util.cross(n, plane.getUp());
		
		//rotate each Pol
		for(int k=0; k < EInPlane.length; k++){
			EInPlane[k] = Pol.rotateFrame(EInPlane[k], Util.dot(sPolDir, plane.getUp()), Util.dot(sPolDir, right));
		}

		return EInPlane;
	}
	
	
}
