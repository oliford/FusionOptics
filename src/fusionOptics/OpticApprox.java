package fusionOptics;

import fusionOptics.types.Material;
import net.jafama.FastMath;

/** Static calc utils for approximations and ideal cases */
public abstract class OpticApprox {

	/** Calculates the optical path difference between the ordinary and extraordinary waves
	 * leaving a wave/displacer plate.
	 * 
	 *  [F.E.Veiras "Phase shift formulas in uniaxial media: an application to waveplates" 
	 *  	 APPLIED OPTICS (49) 15 2010 ] - equation 12
	 *
	 * @param theta		Angle between plate surface and optic axis
	 * @param delta		Angle between plane of incidence and optic axis projected onto plate surface.
	 * @param alpha		Angle of incidence on first surface.
	 * @param no 		Ordinary refractive index.
	 * @param ne		Extraordinary refractive index.
	 * @param n			Refractive index of external medium.
	 * @param L			Thickness of plate.
	 * @return
	 */
	public final static double waveplateOPD(double n, double no, double ne, double theta, double delta, double alpha, double L){
		return waveplateOPD(n, no, ne, 
				FastMath.sin(theta), FastMath.cos(theta), 
				FastMath.sin(delta), FastMath.cos(delta), 
				FastMath.sin(alpha), L);
	}
	
	/** See other waveplateOPD() */
	public final static double waveplateOPD(double n, double no, double ne, 
							double sinTheta, double cosTheta, 
							double sinDelta, double cosDelta, 
							double sinAlpha, double L){
				
		double X = (ne*ne*sinTheta*sinTheta + no*no*cosTheta*cosTheta);
		
		return L*( 
						  FastMath.sqrt(no*no - n*n*sinAlpha*sinAlpha)
						  //The next line is + in the paper but - from my derivation and - fixes it to match
						  //the basic geometrical calc and also Avendano's calc
						- n*(no*no - ne*ne)*sinTheta*cosTheta*cosDelta*sinAlpha / X
						- no*FastMath.sqrt(
								ne*ne*X - (ne*ne - (ne*ne - no*no)*cosTheta*cosTheta*sinDelta*sinDelta)*n*n*sinAlpha*sinAlpha 
							) / X 
				 );
				
	}
}
