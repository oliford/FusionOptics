package fusionOptics.pointSpread;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import fusionOptics.Util;
import fusionOptics.collection.IntersectionProcessor;
import fusionOptics.surfaces.Plane;
import fusionOptics.types.Intersection;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import uk.co.oliford.jolu.BinaryMatrixFile;
import uk.co.oliford.jolu.BinaryMatrixWriter;
import uk.co.oliford.jolu.OneLiners;



/**
 * Handles both collection of points when build a PSF and the collection
 * of PSFs from different source points.
 * 
 * @author oliford
 */
public class PointSpreadBuilder implements IntersectionProcessor {
	
	private class LightCollection {
		public double pos[];
		public double startDir[];
		public double startUp[];		
		public double E[][];
		@Override
		public LightCollection clone() {
			LightCollection lc = new LightCollection();
			lc.pos = pos.clone();
			lc.startDir = startDir.clone();
			lc.startUp = startUp.clone();
			lc.E = new double[E.length][];
			for(int i=0; i < E.length; i++)
				lc.E[i] = E[i].clone();
			return lc;
		}
	}
	
	/** List of points in current PSF building set. Already coherently summed. */
	private LinkedList<LightCollection> currentPSFPoints = new LinkedList<LightCollection>();
	private LinkedList<LightCollection> currentCoherencySetPoints = new LinkedList<LightCollection>();
	
	private double currentSource[];
	
	private PointSpreadFunction psf;
	private String dataOutFileName;
	private BinaryMatrixWriter psfDataOut;
	private Plane polarisationPlane;
	
	/** Radius of circle on image plane in which things are coherently added.
	 * Default behaviour is to assume that everything in a single PSF collection
	 * will add coherently */
	private double maxCoherentIntegrationRadius = Double.POSITIVE_INFINITY;
		
	public PointSpreadBuilder(Plane polarisationPlane) {
		this.polarisationPlane = polarisationPlane;
		this.dataOutFileName = null;
		
	}
	public PointSpreadBuilder(Plane polarisationPlane, String dataOutFileName) {
		this.polarisationPlane = polarisationPlane;
		this.dataOutFileName = dataOutFileName;
	}
	
	/** Start collecting of points for characterisation of a PSF */
	public void startNewPSF(double sourcePos[], PointSpreadFunction newPSF){
		this.currentSource = sourcePos;
		this.psf = newPSF;
		currentPSFPoints.clear();
	}
	
	/** Completes the current PSF */
	public PointSpreadFunction psfDone(boolean addToCache) {
		
		double pos[][] = new double[currentPSFPoints.size()][];
		double E[][][] = new double[currentPSFPoints.size()][][];
		double startDirs[][] = new double[currentPSFPoints.size()][];
		double startUps[][] = new double[currentPSFPoints.size()][];
		
		int i=0;
		for(LightCollection posData : currentPSFPoints){
			pos[i] = posData.pos;
			startDirs[i] = posData.startDir;
			startUps[i] = posData.startUp;
			E[i] = posData.E;
			//coherencySet = (Integer)posData[4];
			i++;
		}
		
		psf.setPoints(pos, E);
		psf.setSourceInfo(startDirs, startUps, E);
		psf.calcAverageMuellerMatrix(E);
		
		if(dataOutFileName != null){
			double psfData[] = psf.getCharacterisationData();
			if(psfData == null)
				throw new RuntimeException("fail");
			if(psfDataOut == null)
				psfDataOut = new BinaryMatrixWriter(dataOutFileName, psfData.length);
			psfDataOut.writeRow(psfData);
		}
		
		return psf;
	}
	
	/** Dumps the recorded points, the completed PSF's PDF and some regenerated points */ 
	public void dumpPSFInfo(String outPath, int nPoints){
		//points
		double pos[][] = new double[currentPSFPoints.size()][];		
		int i=0;
		for(LightCollection posData : currentPSFPoints) {
			pos[i] = new double[posData.E.length*4+2];
			pos[i][0] = posData.pos[0];
			pos[i][1] = posData.pos[1];
			pos[i][2] = posData.startDir[0];
			pos[i][3] = posData.startDir[1];
			pos[i][4] = posData.startDir[2];
			pos[i][5] = posData.startUp[0];
			pos[i][6] = posData.startUp[1];
			pos[i][7] = posData.startUp[2];
			for(int j=0; j < posData.E.length; j++){
				pos[i][8+j*4+0] = posData.E[j][0];
				pos[i][8+j*4+1] = posData.E[j][1];
				pos[i][8+j*4+2] = posData.E[j][2];
				pos[i][8+j*4+3] = posData.E[j][3]; 
			}
			i++;
		}
		BinaryMatrixFile.mustWrite(outPath + "/points.bin", pos, false);
		
		if(psf == null){
			System.err.println("WARNING: PointSpreadCollector.dumpPSFInfo(): Points only, PSF not completed.");
			return;
		}
			
		// PDF
		double x[] = OneLiners.linSpace(psf.getMinX(), psf.getMaxX(), 105);
		double y[] = OneLiners.linSpace(psf.getMinY(), psf.getMaxY(), 95);
		BinaryMatrixFile.mustWrite(outPath + "/pdf.bin", y, x, psf.getGridProbability(x, y), true);
		
		// regen points
		BinaryMatrixFile.mustWrite(outPath + "/regenPoints.bin", psf.generatePoints(nPoints), true);
	}

	/** Call when starting to trace a ray tree to signal that all the rays between calls of this
	 * can be added coherently. i.e, they came from a coherent source (because they were a single ray)
	 */
	public void nextCoherentSet() {
		double maxR2=maxCoherentIntegrationRadius*maxCoherentIntegrationRadius;
		
		if(Double.isInfinite(maxR2)) { // add them all incoherently
			currentPSFPoints.addAll(currentCoherencySetPoints);
			currentCoherencySetPoints.clear();
			return;
		}
		
		//double totalPowerA = 0;
		//for(int i=0; i < currentCoherencySetPoints.size(); i++){
		//	totalPowerA = Pol.intensity( ((double[][])(currentCoherencySetPoints.get(i)[3])) );
		//}
		
		//double totalPowerB = 0;
		
		LightCollection lc1;
		//get the next point from the collected set
		while((lc1 = currentCoherencySetPoints.poll()) != null){
			lc1 = lc1.clone(); //need copies because of the Pol allocs
			//totalPowerB += Pol.intensity(E);
			
			//scan the whole set (excluding this point because we polled it already), 
			//  for points that can be added coherently
			for(int i=0; i < currentCoherencySetPoints.size(); i++){
				LightCollection lc2 = currentCoherencySetPoints.get(i);
								
				double d2 = (lc2.pos[0]-lc1.pos[0])*(lc2.pos[0]-lc1.pos[0]) +
									(lc2.pos[1]-lc1.pos[1])*(lc2.pos[1]-lc1.pos[1]);
				
				//add them
				if(d2 < maxR2){
					for(int j=0; j < lc2.E.length; j++){
						lc1.E[j][0] += lc2.E[j][0];
						lc1.E[j][1] += lc2.E[j][1];
						lc1.E[j][2] += lc2.E[j][2];
						lc1.E[j][3] += lc2.E[j][3];
					}
					//and remove from the source set
					currentCoherencySetPoints.remove(i);
				}
			}
			
			//add that sum
			currentPSFPoints.add(lc1);
		}
		
	}

	@Override
	/** Adds ray/plane intersections (with polarisation) to the current PSF */
	public void nextIntersection(Intersection imgHit) {
		
		LightCollection lc = new LightCollection();
		
		Plane imagePlane = (Plane)imgHit.surface;
		lc.pos = imagePlane.posXYZToPlaneUR(imgHit.pos);
		
		//Walk backwards until we find the last time it went through the polarisation sensitive plane
		RaySegment startRay = imgHit.incidentRay;
		Intersection polHit = null;
		
		do {
			if(startRay.endHit.surface == polarisationPlane)
				polHit = startRay.endHit;
			if(startRay.startHit == null)
				break;
			startRay = startRay.startHit.incidentRay;
		}while(true);
		
		if(polHit == null)
			return; //never found it and got to start of ray
		
		//make sure the sense of polarisation on the pol plane incident ray 
		//matches the pol plane's sense of up
		//polHit.incidentRay.rotatePolRefFrame(polarisationPlane.getUp());
		
		double EInPlane[][] = Pol.projectToPlanesView(polHit, false);
		
		lc.startDir = startRay.dir.clone();
		lc.startUp = startRay.up.clone();
		lc.E = EInPlane;
		
		currentCoherencySetPoints.add(lc);
	}
	
	public int getNPointsCollected() {
		return currentPSFPoints.size();
	}

	public void setMaxCoherentIntegrationRadius(double maxCoherentIntegrationRadius){
		this.maxCoherentIntegrationRadius = maxCoherentIntegrationRadius;
	}
	
	public double[][] getCollectedPoints() {
		double pos[][] = new double[currentPSFPoints.size()][];
		int i = 0;
		for(LightCollection lc : currentPSFPoints) {
			pos[i] = lc.pos.clone();
			i++;
		}		
		return pos;
	}
	
}
