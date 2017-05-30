package fusionOptics.collection;

import fusionOptics.Util;
import fusionOptics.surfaces.Plane;
import fusionOptics.types.Intersection;
import binaryMatrixFile.BinaryMatrixWriter;
import jafama.FastMath;

/** Collects statistics of the angle of rays through a given surface */
public class PlaneAngleInfo implements IntersectionProcessor {
	private Plane targetPlane;
	
	private double sumN, sumU, sumR, sumI;
	private double minAng, maxAng;
	BinaryMatrixWriter dbgOut = null;
	
	public PlaneAngleInfo(Plane target) {
		this.targetPlane = target;
		reset();
	}
	
	@Override
	public void nextIntersection(Intersection hit) {
		Intersection targetHit = (hit.surface == targetPlane) ? hit : hit.incidentRay.findFirstEarlierIntersection(targetPlane);
		if(targetHit == null)
			return;
		
		double I = hit.incidentRay.endIntensity();
		double N = Util.dot(targetHit.incidentRay.dir, targetPlane.getNormal());
		double U = Util.dot(targetHit.incidentRay.dir, targetPlane.getUp());
		double R = Util.dot(targetHit.incidentRay.dir, targetPlane.getRight());
		sumN += I * N;
		sumU += I * U;
		sumR += I * R;
		sumI += I;
		
		double ang = FastMath.atan2(FastMath.sqrt(U*U + R*R), N);
		if(ang < minAng)
			minAng = ang;
		if(ang > maxAng)
			maxAng = ang;
		
		if(dbgOut != null)
			dbgOut.writeRow(I, N, U, R, targetHit.incidentRay.dir);
		
	}
	
	public void dump() {
		System.out.println("Angles through plane '"+targetPlane.getName()+"': mean = " 
						+ getMeanAngle()*180/Math.PI + "°, min = " 
						+ minAng*180/Math.PI + "°, max = " + maxAng*180/Math.PI + "°");			
	}
	
	public double getMeanAngle() {
		return FastMath.atan2(FastMath.sqrt(sumU*sumU + sumR*sumR), sumN);
	}
	
	public double getMinAngle() { return minAng; }	
	public double getMaxAngle() { return maxAng; }
	
	/** @return the mean vector at the intersection, in cartesian space */
	public double[] getMeanVector() {
		double N[] = targetPlane.getNormal();
		double U[] = targetPlane.getUp();
		double R[] = targetPlane.getRight();
		
		return Util.reNorm(new double[]{
				sumN/sumI * N[0] + sumU/sumI * U[0] + sumR/sumI * R[0],
				sumN/sumI * N[1] + sumU/sumI * U[1] + sumR/sumI * R[1],
				sumN/sumI * N[2] + sumU/sumI * U[2] + sumR/sumI * R[2]  
		});
		
	}
	
	public void reset(){
		sumN = 0;
		sumU = 0;
		sumR = 0;
		sumI = 0;
		minAng = Double.POSITIVE_INFINITY;
		maxAng = Double.NEGATIVE_INFINITY;
	}
	
	public void setDbgOut(String fileName) {
		if(this.dbgOut != null)
			this.dbgOut.close();
		this.dbgOut = new BinaryMatrixWriter(fileName, 7);
	}
	
}
