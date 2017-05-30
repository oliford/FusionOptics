package fusionOptics.collection;

import fusionOptics.surfaces.Plane;
import fusionOptics.types.Intersection;
import fusionOptics.types.Pol;

/** Collects simple hit position stats (mean, covariance) with a surface */ 
public class HitPositionAverage implements IntersectionProcessor {

	public boolean byIntensity;
		
	public HitPositionAverage(boolean byIntensity) { this.byIntensity = byIntensity; }
	public HitPositionAverage() { this.byIntensity = false; }
	
	public double sumIR, sumIU, sumIRR, sumIRU, sumIUU, sumI;
	
	@Override
	public void nextIntersection(Intersection hit) {
		double posRU[] = ((Plane)hit.surface).posXYZToPlaneRU(hit.pos);
		
		double I = byIntensity ? Pol.intensity(hit.incidentRay.E1) : 1;
		
		sumI += I;
		sumIR += I * posRU[0];
		sumIU += I * posRU[1];
		sumIRR += I * posRU[0] * posRU[0];
		sumIRU += I * posRU[0] * posRU[1];
		sumIUU += I * posRU[1] * posRU[1];
		
	}
	
	public double getSumI(){ return sumI; }
	
	public double getMeanR(){ return sumIR / sumI; }
	public double getMeanU(){ return sumIU / sumI; }
	
	public double[] getMeanPosRU(){ 
		return new double[]{ sumIR / sumI, sumIU / sumI }; 
	}
	public void reset() {
		sumI = 0; sumIR = 0; sumIU = 0;
		sumIRR = 0; sumIRU = 0; sumIUU = 0;
	}
	
	public double getSigmaRR(){ return (sumIRR - sumIR*sumIR/sumI)/sumI; } 
	public double getSigmaUU(){ return (sumIUU - sumIU*sumIU/sumI)/sumI; }
		
}
