package fusionOptics.collection;

import fusionOptics.surfaces.Plane;
import fusionOptics.types.Intersection;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import binaryMatrixFile.BinaryMatrixWriter;

/** Writes full hit information straight to a Binary matrix file */
public class HitsCollector implements IntersectionProcessor {
	
	public int iL;
	public int iA, iB;
	public int nHits;
	
	/** Colour table for rays */
	public double col[];
	
	/** Output for individual hit data (useful if only 1 start point) */
	public BinaryMatrixWriter hitsOut = null;
	
	public Plane polarisationPlane;
	
	public HitsCollector(String fileName, Plane polarisationPlane) {
		this.polarisationPlane = polarisationPlane;
		this.hitsOut = new BinaryMatrixWriter(fileName, 20);
	}
	
	@Override
	public void nextIntersection(Intersection imgHit) {
		
		//Walk backwards until we find the last time it went through the polarisation sensitive plane
		Intersection polHit = imgHit.incidentRay.findFirstEarlierIntersection(polarisationPlane);
		if(polHit == null)
			return;
		
		RaySegment startRay = polHit.incidentRay;
		while(startRay.startHit != null && startRay.startHit.incidentRay != null)
			startRay = startRay.startHit.incidentRay;
		
		Plane imagePlane = (Plane)polHit.surface;
		
		//make sure the sense of polarisation on the pol plane incident ray 
		//matches the pol plane's sense of up
		polHit.incidentRay.rotatePolRefFrame(polarisationPlane.getUp());
		
		double imgPos[] = imagePlane.posXYZToPlaneUR(polHit.pos);
		hitsOut.writeRow(
				iL,	iA, iB,									// 1: startIdxes (L,A,B)
				startRay.startPos, 								// 4: start posXYZ 
				imgHit.pos, 								// 7: hit posXYZ 
				imgPos,  									//10: image position
				(col != null) ? col : new double[]{0,0,0}, //12: colour, (to match drawing)
				Pol.intensity(imgHit.incidentRay.E1[0]),	//15: total intensity
				Pol.psi(polHit.incidentRay.E1[0]),			//16: rotation dir of light originally 'up'
				Pol.chi(polHit.incidentRay.E1[0]),			//17: ellipticity of light originally 'up'
				Pol.intensity(imgHit.incidentRay.E1[1]),	//18: total intensity
				Pol.psi(polHit.incidentRay.E1[1]),			//19: rotation dir of light originally 'right'						
				Pol.chi(polHit.incidentRay.E1[1])			//20: ellipticity of light originally 'right'				
				
			);
		nHits++;
	}

	public void destroy() {
		hitsOut.close();
	}
	
}

