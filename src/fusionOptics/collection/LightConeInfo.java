package fusionOptics.collection;

import java.util.ArrayList;

import oneLiners.OneLiners;
import fusionOptics.Util;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Plane;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Optic;
import algorithmrepository.Algorithms;
import algorithmrepository.LinearInterpolation1D;
import binaryMatrixFile.BinaryMatrixWriter;
import jafama.FastMath;

/** Collects all vectors and positions of the light hitting a surface.
 * Can be used to find the approx optimal focus point and max/average light cones.
 * 
 * Useful for dealing with optic fibres
 * */
public class LightConeInfo implements IntersectionProcessor {
	private Plane targetPlane;
	
	private static class IncidentRay {
		double I;
		double pos[];
		double dir[];
	};
	
	/** vector is [ Intensity, hitPos[0,1,2], incidenceVec[0, 1, 2] ] */
	public ArrayList<IncidentRay> rays = new ArrayList<IncidentRay>();  
		
	public LightConeInfo(Plane target) {
		this.targetPlane = target;
		reset();
	}
	
	@Override
	public void nextIntersection(Intersection hit) {
		Intersection targetHit = (hit.surface == targetPlane) ? hit : hit.incidentRay.findFirstEarlierIntersection(targetPlane);
		if(targetHit == null)
			return;
		
		IncidentRay iRay = new IncidentRay();
		iRay.I = hit.incidentRay.endIntensity();
		iRay.pos = hit.pos.clone();
		iRay.dir = hit.incidentRay.dir.clone();
		
		rays.add(iRay);
	}
	
	public void dump() {
		
	}
	
	public double[] getRayAngles(){
		double meanVec[] = getMeanVector();
		double angles[] = new double[rays.size()];
		for(int i=0; i < rays.size(); i++){
			angles[i] = FastMath.acos(Util.dot(rays.get(i).dir, meanVec));
		}
		return angles;
	}
	
	public double getMeanAngleToPlane() {  return FastMath.acos(Util.dot(getMeanVector(), targetPlane.getNormal())); }	
	public double getMinAngleToPlane() {
		double planeNormal[] = targetPlane.getNormal();
		double minAngle = Double.POSITIVE_INFINITY;
		for(IncidentRay iRay : rays){
			double angle = FastMath.acos(Util.dot(iRay.dir, planeNormal));
			if(angle < minAngle)
				minAngle = angle;
		}
		return minAngle;
	}	
	
	public double getMaxAngleToPlane() {  
		double planeNormal[] = targetPlane.getNormal();
		double maxAngle = Double.NEGATIVE_INFINITY;
		for(IncidentRay iRay : rays){
			double angle = FastMath.acos(Util.dot(iRay.dir, planeNormal));
			if(angle > maxAngle)
				maxAngle = angle;
		}
		return maxAngle;
	}
	
	/** @return the mean vector at the intersection, in cartesian space */
	public double[] getMeanVector() {
		double dir[] = new double[3];
		for(IncidentRay iRay : rays){
			dir[0] += iRay.dir[0];
			dir[1] += iRay.dir[1];
			dir[2] += iRay.dir[2];
		}
		int n = rays.size();
		return Util.reNorm(new double[]{ dir[0] / n, dir[1] / n, dir[2] / n });
	}
	
	/** Finds the acceptance angle (from the mean vector) required to collect the given fraction of rays */  
	public double getCapturingAngle(double rayFraction){
		double angles[] = getRayAngles();
		Algorithms.quicksort(angles, null);
		
		double fraction[] = OneLiners.linSpace(0.0, 1.0, angles.length);
		LinearInterpolation1D interp = new LinearInterpolation1D(fraction, angles);
		return interp.eval(rayFraction);		
	}
	
	/** Gets the average position of the closest approach of all pairs of rays */
	public double[] getApproxFocusPos(){
		double focusPos[] = new double[3];
		int n=0;
				
		for(IncidentRay iRay1 : rays){
			for(IncidentRay iRay2 : rays){
				if(iRay1 == iRay2)
					continue;
		
				double s = Algorithms.pointOnLineNearestAnotherLine(iRay1.pos, iRay1.dir, iRay2.pos, iRay2.dir);
				if(Double.isNaN(s)){
					//what??
					continue;
				}

				focusPos[0] += iRay1.pos[0] + s * iRay1.dir[0];
				focusPos[1] += iRay1.pos[1] + s * iRay1.dir[1];
				focusPos[2] += iRay1.pos[2] + s * iRay1.dir[2];
				n++;
				
			}
		}

		return new double[]{ focusPos[0]/n, focusPos[1]/n, focusPos[2]/n };
		
	}
	
	public double getFractionInMeanCone(double maxAngle){
		return getFractionInCone(getMeanVector(), maxAngle);
	}
	
	public double getFractionInPlanePerpCone(double maxAngle){
		return getFractionInCone(targetPlane.getNormal(), maxAngle);
	}
	
	public double getFractionInCone(double coneAxis[], double coneAngle){
		int n = 0;
		for(IncidentRay iRay : rays){
			double angle = FastMath.acos(FastMath.abs(Util.dot(iRay.dir, coneAxis)));
			if(angle < coneAngle)
				n++;
		}
		return (double)n / rays.size();
	}
	
	public Optic makeFibreCylinder(double length, double radius, double fixVector[]){
		double focusPos[] = getApproxFocusPos();
		double meanVec[] = (fixVector != null) ? fixVector : getMeanVector();			
				
		double pos[] = Util.plus(focusPos, Util.mul(meanVec, length/2));
		Cylinder cyld = new Cylinder("fibreCyld", pos, meanVec, radius, length, NullInterface.ideal());
		
		double angle90 = getCapturingAngle(0.90);
		
		Optic optic = new Optic("LightCone", new Element[]{ cyld } );
		for(int i=1; i < 4; i++){
			double coneRadius = FastMath.tan(angle90) * i * length;
			Disc disc = new Disc("disc"+i, Util.plus(focusPos, Util.mul(meanVec, -i*length)),
										meanVec, coneRadius, NullInterface.ideal());
			optic.addElement(disc);
		}
		
		return optic;
	}
	
	
	public void reset(){
		rays.clear();
	}
	
}
