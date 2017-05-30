package fusionOptics.media;

import java.util.Arrays;

import fusionOptics.Util;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;

import jafama.FastMath;

/** Medium immersed in a homogeneous magnetic field.
 * 
 * Applied faraday rotation
 * 
 * @author oliford
 */

public class MediumInMagneticField extends Medium {

	private double field[];
	
	public MediumInMagneticField(Material material) {
		super(material, null, 300);
		this.field = null;
	}
	
	public MediumInMagneticField(Material material, double[][] opticAxes, double temperature, double field[]) {
		super(material, opticAxes, temperature);
		this.field = field;
	}
	
	public void setField(double[] field) { this.field = field; }
	
	public double[] getField() { return this.field; }

	/** Indecies for E[] in terms of left/right circ polarised waves */
	private final static int lcRe = 0, lcIm = 1, rcRe = 2, rcIm = 3;
    
	@Override
	public void rayPropagation(RaySegment ray) {
		
		if(material == null){
			Medium.rayPropagation(ray, 1.0, 1.0, 1.0, 1.0);
			return;
		}
			
		if(material.getNAxes() != 0) {
			throw new IllegalArgumentException("MediumInMagneticFiled (Faraday rotation) can currently only cope with Isotropic materials");
			// It might actually be ok to let IsoUniaxialInterface split the rays and then just apply this for each of the E and O 
		    // rays ignoring the birefringence.
		}
		
		double Bpara = (field == null) ? 0 : Util.dot(ray.dir, field);
		if(Bpara == 0){
			super.rayPropagation(ray);
			return;
		}
		
		double n,t;
		if(Double.isNaN(ray.raySpecificRefractiveIndex) || ray.raySpecificRefractiveIndex == 0){
			//ordinary ray
			n = material.getRefractiveIndex(0, ray.wavelength, temperature);
			t = material.getTransmission(0, ray.wavelength, temperature)+0;
		
		}else{ //extraordinary
			n = ray.raySpecificRefractiveIndex;
			t = material.getTransmission(1, ray.wavelength, temperature)+0;
		
		}
		t = FastMath.pow(t, ray.length);
		
		
		double V = material.getVerdetConstant(0, ray.wavelength, temperature);
		
		//difference in L, R refractive indecies
		//V.B = Δn . ω / 2 c  
		//Δn = 2 c V B / ω
		//Δn = V B λ0 / π
		
		double dN = V * Bpara * ray.wavelength / Math.PI; 
				
		
		//propagate the L,R circ polarised waves with the repsective wavelength
		//use the average for the total wave count, and put the difference in the phases
		//
		ray.nWaves = (int)(ray.length * n / ray.wavelength);
		double partWaveCommon = ray.length * n / ray.wavelength - ray.nWaves;
		
		double partWaveL = partWaveCommon + ray.length * dN/2 / ray.wavelength;
		double partWaveR = partWaveCommon - ray.length * dN/2 / ray.wavelength;
	
		//and change the phases of the polarisations
		double cosPhaseL = FastMath.cos(2 * Math.PI * partWaveL);
		double sinPhaseL = FastMath.sin(2 * Math.PI * partWaveL);
		double cosPhaseR = FastMath.cos(2 * Math.PI * partWaveR);
		double sinPhaseR = FastMath.sin(2 * Math.PI * partWaveR);
		
		
		
		ray.E1 = Pol.alloc(ray.E0.length);
		// Break into L and R circ waves
        double Elr[] = new double[4];
		for(int i=0; i < ray.E0.length; i++){               
            Elr[lcRe] = (ray.E0[i][Pol.uRe] - ray.E0[i][Pol.rIm])/2; //Re(L)
            Elr[lcIm] = (ray.E0[i][Pol.uIm] + ray.E0[i][Pol.rRe])/2; //Im(L)
            Elr[rcRe] = (ray.E0[i][Pol.uRe] + ray.E0[i][Pol.rIm])/2; //Re(R)
            Elr[rcIm] = (ray.E0[i][Pol.uIm] - ray.E0[i][Pol.rRe])/2; //Im(R)            
        //}
        
		//Elr = Pol.complexMulAll(Elr, cosPhaseL * t, sinPhaseL * t, cosPhaseR * t, sinPhaseR * t, 1);
        Elr = Pol.complexMul(Elr, cosPhaseL * t, sinPhaseL * t, cosPhaseR * t, sinPhaseR * t, 1);
		
		//and convert back to Exy
		//ray.E1 = new double[ray.E0.length][4];
		//for(int i=0; i < ray.E0.length; i++){						
			ray.E1[i][Pol.uRe] = Elr[lcRe] + Elr[rcRe]; //Re(X)
			ray.E1[i][Pol.uIm] = Elr[lcIm] + Elr[rcIm]; //Im(X)
			ray.E1[i][Pol.rRe] = Elr[lcIm] - Elr[rcIm]; //Re(Y)
			ray.E1[i][Pol.rIm] = Elr[rcRe] - Elr[lcRe]; //Im(Y)			
		}	
		
		
	}

	@Override
	public int hashCode() {
		return super.hashCode() * 31 + Arrays.hashCode(field);
	}
	
	@Override
	public boolean equals(Object obj) {		
		return super.equals(obj) && Arrays.equals(field, ((MediumInMagneticField)obj).field);
	}
}
