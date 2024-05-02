package fusionOptics.collection;

import fusionOptics.surfaces.Plane;
import fusionOptics.types.Intersection;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import uk.co.oliford.jolu.BinaryMatrixWriter;

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
		this.hitsOut = new BinaryMatrixWriter(fileName, 14);
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
				(col != null) ? col : new double[]{0,0,0} //12: colour, (to match drawing)
				
			);
		nHits++;
	}

	public void destroy() {
		hitsOut.close();
	}
	
}

