package fusionOptics.types;

import java.util.Arrays;

import fusionOptics.Util;
import fusionOptics.interfaces.IsoUniaxialInterface;

import algorithmrepository.exceptions.NotImplementedException;
import net.jafama.DoubleWrapper;
import net.jafama.FastMath;

/** Describes an instance of a material (in space). 
 * i.e. it has a material, orientation, optic axis info
 * 
 *  
 *  
 * @author oliford
 *
 */
public class Medium {
	/** Maximum angle at which to to consider a ray is propagating exactly parallel or perpendicular
	 * to the optic axis.  */
	private double sinMaxAxisAlignAngle = FastMath.sin(0.4 * Math.PI / 180); //0.3 degrees
	
	protected Material material;
	
	/** in Kelvin */
	protected double temperature;
	
	/** The optic axes of the respective refractive indices. Can be null for isotropic media */
	protected double opticAxes[][];
	
	/** Creates a medium of the given material at room temperature (300K)
	 * If the material is uniaxial or biaxial, the first optic axis will be in the x direction.
	 * If the matieral is biaxial, the second optic axis will be in the y direction.
	 * 
	 * @param material
	 */
	public Medium(Material material){
		this.material = material;
		this.temperature = 300;
		switch(material.getNAxes()){			
			case 1: opticAxes = new double[][]{ { 1, 0, 0 } }; break;
			case 2: opticAxes = new double[][]{ { 1, 0, 0 }, { 0, 1, 0 } }; break;
			default: opticAxes = null; 
		}	
	}

	/** Create a medium with the given medium, optic axes and temperature */
	public Medium(Material material, double opticAxes[][], double temperature) {
		this.material = material;
		this.opticAxes = opticAxes;
		this.temperature = temperature;
	}
	
	/** Return the material object for this */
	public Material getMaterial(){ return material; }
	
	/** Returns the medium temperature in Kelvin */
	public double getTemperature(){ return temperature; }
	
	/** Sets the medium temperature in Kelvin */
	public void setTemperature(double temperature){ this.temperature = temperature; }
	
	/** Rotate the medium optical axes / tensors etc */
	public void rotate(double point[], double matrix[][]) {
		if(opticAxes == null) //isotropic
			return;
		double newAxes[][] = new double[opticAxes.length][3];
		for(int i=0;i<3;i++) {
			for(int j=0;j<3;j++) {
				for(int k=0;k<opticAxes.length; k++) {
					newAxes[k][i] += matrix[i][j] * opticAxes[k][j];
				}
			}
		}
		opticAxes = newAxes;
	}

	/* All the 'material' functions are copied through here, taking our temperature (etc, that we might invent later */
	
	/** Returns the refractive index of the given mode at the medium's current temperature */
	public double getRefractiveIndex(int modeNumber, double wavelength) {		
		return material.getRefractiveIndex(modeNumber, wavelength, temperature);
	}
	
	/** Returns the refractive index at the medium's current temperature for a ray of the given wavelength 
	 * and direction, polarised perpendicular to the optic axis.
	 * 
	 * @param dir
	 * @param wavelength
	 * @return
	 */
	/* public double getBirefringentExtraordinaryIndex(double dir[], double wavelength) {		
		double no = material.getRefractiveIndex(0, wavelength, temperature);
		double ne = material.getRefractiveIndex(1, wavelength, temperature);
		
		double
		blerg
		
		 Angle would be ray direction
		
	} */
	
	/** Returns the transmission coefficient at the medium's current temperature */
	public double getTransmission(int modeNumber, double wavelength) {
		return material.getTransmission(modeNumber, wavelength, temperature);
	}
	
	/** Returns the verdet constant for the given material at the given wavelength and the medium's current temperature */
	public double getVerdetConstant(int modeNumber, double wavelen) {
		return material.getVerdetConstant(modeNumber, wavelen, temperature);
	}
	
	/** Fills in the ending polarisation of a ray given the medium's effect on it's propagation */
	public void rayPropagation(RaySegment ray){
		
		//for isotropic media, the refractive indices are the same, so the polarisation frame
		//can stay as it is
		if(material == null)
			rayPropagation(ray, 1.0, 1.0, 1.0, 1.0);
		else if(material.getNAxes() == 0) { // Isotropic Media
			
			double n = material.getRefractiveIndex(0, ray.wavelength, temperature);
			double t = material.getTransmission(0, ray.wavelength, temperature)+0;
			t = FastMath.pow(t, ray.length);
			rayPropagation(ray, n, n, t, t);
			
		}else if(material.getNAxes() == 1) { // Uniaxial Media
			
			if(IsoUniaxialInterface.simpleBirefringence){
				// If the optic axis is perp or parallel to the ray direction,
				// we can handle it as the simpler specific cases of two polarisations in one 
				// RaySegment, since the ray directions are the same
				double no = material.getRefractiveIndex(0, ray.wavelength, temperature);
				double to = material.getTransmission(0, ray.wavelength, temperature)+0;
			
				double dirDotOptic = FastMath.abs(Util.dot(opticAxes[0], ray.dir));
				if((1.0 - dirDotOptic) < sinMaxAxisAlignAngle){
					//Transmission almost exactly parallel to optic axis is the 
					// same as isotropic, since both polarisations are perp to it
					double n = material.getRefractiveIndex(0, ray.wavelength, temperature);
					rayPropagation(ray, no, no, to, to);
					
				}else if(dirDotOptic < sinMaxAxisAlignAngle){
					//Transmission almost exactly perp to optic axis is also quite simple
					// but we need to rotate the polarisation states to the optic axis frame
					//first
					double ne = material.getRefractiveIndex(1, ray.wavelength, temperature);
					double te = material.getTransmission(0, ray.wavelength, temperature)+0;
					te = FastMath.pow(te, ray.length);
					
					ray.rotatePolRefFrame(opticAxes[0]); //up is now extraordinary
					
					rayPropagation(ray, ne, no, te, to);
					
				}else{
					//we can't deal with the general case properly.
					throw new RuntimeException("IsoUniaxialInterface.simpleBirefringence is true, but we found at ray at angle "+
									FastMath.acos(dirDotOptic) *180/Math.PI
									+"° to the optic axis. The max from perp/para is " +
									FastMath.asin(sinMaxAxisAlignAngle)*180/Math.PI + "°");
				}
				
			}else{				
				//For the completely general case, the ray should have been split
				// by IsoUniaxialInterface into a separate RaySegement for the E and O rays
				// which each have only one index
				double n,t;
				if(Double.isNaN(ray.raySpecificRefractiveIndex)
						|| ray.raySpecificRefractiveIndex == 0){
					//ordinary ray
					n = material.getRefractiveIndex(0, ray.wavelength, temperature);
					t = material.getTransmission(0, ray.wavelength, temperature)+0;
				
				}else{ //extraordinary
					n = ray.raySpecificRefractiveIndex;
					t = material.getTransmission(1, ray.wavelength, temperature)+0;
				
				}
				t = FastMath.pow(t, ray.length);
				rayPropagation(ray, n, n, t, t);
			}
				
		}else{ // Biaxial Media

			throw new RuntimeException("Biaxial medium not yet supported.");
		}
				
		
	}
	
	/** Set a ray's 'nWaves' variable and rotate the 'up' and 'right' phases according to
	 * the given refractive indices.
	 * 
	 * @param ray
	 * @param nU	Refractive index for 'up' polarisaion
	 * @param nR	Refractive index for 'right' polarisaion
	 */
	public static void rayPropagation(RaySegment ray, double nU, double nR, double tU, double tR) {
		
		double wavelenInMediumU = ray.wavelength / nU; 
		double wavelenInMediumR = ray.wavelength / nR; 
		
		//use the average for the total wave count, and put the difference in the phases
		//
		ray.nWaves = (int)(ray.length * 2 / (wavelenInMediumU + wavelenInMediumR));		
		
		double partWaveU = ray.length / wavelenInMediumU - ray.nWaves;
		double partWaveR = ray.length / wavelenInMediumR - ray.nWaves;
		
		//and change the phases of the polarisations
		double cosPhaseU = FastMath.cos(2 * Math.PI * partWaveU);
		double sinPhaseU = FastMath.sin(2 * Math.PI * partWaveU);
		double cosPhaseR = FastMath.cos(2 * Math.PI * partWaveR);
		double sinPhaseR = FastMath.sin(2 * Math.PI * partWaveR);
		
		ray.E1 = Pol.complexMulAll(ray.E0, cosPhaseU * tU, sinPhaseU * tU, cosPhaseR * tR, sinPhaseR * tR, 1);
		
		
		/*
		double phaseU0 = FastMath.atan2(ray.E0[0][1], ray.E0[0][0]);
		double phaseR0 = FastMath.atan2(ray.E0[0][3], ray.E0[0][2]);
		double phaseU1 = FastMath.atan2(ray.E1[0][1], ray.E1[0][0]);
		double phaseR1 = FastMath.atan2(ray.E1[0][3], ray.E1[0][2]);
		
		double deltaPhaseUT = 2 * Math.PI * partWaveU;
		double deltaPhaseRT = 2 * Math.PI * partWaveR;
		double deltaPhaseUC = phaseU1 - phaseU0;
		double deltaPhaseRC = phaseR1 - phaseR0;
		
		System.out.println("--------------------------------------");
		System.out.println("phaseU0 = " + phaseU0*180/Math.PI + "\tphaseU1=" + phaseU1*180/Math.PI);
		System.out.println("phaseR0 = " + phaseR0*180/Math.PI + "\tphaseR1=" + phaseR1*180/Math.PI);
		System.out.println("ΔU: Targ = " + deltaPhaseUT*180/Math.PI + "\tCalc =" + deltaPhaseUC*180/Math.PI);
		System.out.println("ΔR: Targ = " + deltaPhaseRT*180/Math.PI + "\tCalc =" + deltaPhaseRC*180/Math.PI);
		*/
		
	}
	
	@Override
	public Medium clone() {
		Medium newMed = new Medium(material, opticAxes, temperature);
		newMed.sinMaxAxisAlignAngle = sinMaxAxisAlignAngle;
		return newMed;
	}

	public void setOpticAxes(double[][] opticAxes) {
		this.opticAxes = opticAxes;
	}

	public double[][] getOpticAxes() { return opticAxes; }

	@Override
	public int hashCode() {
		int result = 31 + material.hashCode();
		result = 31 * result + Arrays.hashCode(opticAxes);
		long temp = Double.doubleToLongBits(temperature);
		result = 31 * result + (int) (temp ^ (temp >>> 32));
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		Medium other = (Medium) obj;
		return material.equals(other.material)
				&& Arrays.deepEquals(opticAxes, other.opticAxes)
				&& Double.doubleToLongBits(temperature) == Double.doubleToLongBits(other.temperature);
	}
	
}
