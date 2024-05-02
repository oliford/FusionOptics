package fusionOptics.collection;

import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Plane;
import fusionOptics.types.Intersection;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import net.jafama.FastMath;
import uk.co.oliford.jolu.BinaryMatrixFile;
import uk.co.oliford.jolu.BinaryMatrixWriter;
import uk.co.oliford.jolu.OneLiners;

/** Collects an image of the average polarisation at a given polarisation plane,
 *  w.r.t. the position at another given plane.
 *  
 * @author oliford
 */
public class LensDepolarisationInfoCollector implements IntersectionProcessor {

	private Plane posPlane;
	private Plane polarisationPlane;
	private double x[], y[], sumPsi[][], sumI[][];
	private double polUpVec[];
	private int polIdx;
	private String prefix;
	
	BinaryMatrixWriter allPointsOut = null;
	
	public LensDepolarisationInfoCollector(Plane posPlane, Plane polarisationPlane, int nCells, int polIdx, String prefix) {
		this(posPlane, polarisationPlane, ((Plane)polarisationPlane).getUp(), nCells, polIdx, prefix); 
	}
	
	public LensDepolarisationInfoCollector(Plane posPlane, Plane polarisationPlane, double polUpVec[], int nCells, int polIdx, String prefix) {
		this.posPlane = posPlane;
		this.polarisationPlane = polarisationPlane;
		this.polUpVec = polUpVec;
		this.polIdx = polIdx;
		this.prefix = prefix;
		
		this.allPointsOut = new BinaryMatrixWriter(prefix + "-points.bin", 4);
		
		double maxRad = posPlane.getBoundarySphereRadius();
		x = OneLiners.linSpace(-maxRad, maxRad, nCells+1);
		y = OneLiners.linSpace(-maxRad, maxRad, nCells-1);
		sumI = new double[y.length][x.length];
		sumPsi = new double[y.length][x.length];
	}
	
	@Override
	public void nextIntersection(Intersection imgHit) {
		//Walk backwards until we find the last time it went through the polarisation sensitive plane
		RaySegment startRay = imgHit.incidentRay;
		Intersection posPlaneHit = null;
		Intersection polPlaneHit = null;
		
		do {
			if(startRay.endHit.surface == posPlane)
				posPlaneHit = startRay.endHit;
			if(startRay.endHit.surface == polarisationPlane)
				polPlaneHit = startRay.endHit;
			if(startRay.startHit == null)
				break;
			startRay = startRay.startHit.incidentRay;
		}while(true);
		
		if(posPlaneHit == null || polPlaneHit == null)
			return; //never found it and got to start of ray
		
		double xy[] = posPlane.posXYZToPlaneUR(posPlaneHit.pos);
		
		//put final ray in sense of plane
		polPlaneHit.incidentRay.rotatePolRefFrame(polUpVec);

		double dx = x[1] - x[0], dy = y[1] - y[0];
		int iX = (int)((xy[0] - x[0]) / dx), iY = (int)((xy[1] - y[0]) / dy);
		
		// This seems a bit weird, but we want the polarisation angle at the polarisation
		// plane, weighted in averaging by the intensity at the image plane, on the assumption,
		// as always, that the polarisation got recoded as something else (like time or spatial variation)
		// at the pol plane, but was actually counted at the imaging plane (CCD, PMT, etc)
		double I = Pol.intensity(imgHit.incidentRay.E1[polIdx]);
		double psi = Pol.psi(Pol.projectToPlanesView(polPlaneHit, false)[polIdx]);
		//double psi = Pol.psi(polPlaneHit.incidentRay.E1[polIdx]);
		
		sumI[iY][iX] += I;
		//sumPsi[iY][iX] += Pol.psi(polPlaneHit.incidentRay.E1[polIdx]) * I;
		sumPsi[iY][iX] += psi * I;
		if(allPointsOut != null){
			allPointsOut.writeRow(xy, I, psi);
		}
	
	}

	public void write() {
		for(int i=0; i < y.length; i++)
			for(int j=0; j < x.length; j++)
				sumPsi[i][j] /= sumI[i][j];
					
		BinaryMatrixFile.mustWrite(prefix + "-I.bin", x, y, sumI, false);
		BinaryMatrixFile.mustWrite(prefix + "-avgPsi.bin", x, y, sumPsi, false);
		
		if(allPointsOut != null)
			allPointsOut.close();
	}
}