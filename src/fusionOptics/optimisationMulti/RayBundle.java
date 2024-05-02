package fusionOptics.optimisationMulti;


import java.util.ArrayList;
import java.util.List;

import net.jafama.FastMath;
import uk.co.oliford.jolu.BinaryMatrixWriter;
import uk.co.oliford.jolu.RandomManager;
import fusionOptics.Util;
import fusionOptics.collection.IntersectionProcessor;
import fusionOptics.drawing.RayDrawer;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.surfaces.Plane;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;


/** A Bundle of rays fired in the optimiser */
public class RayBundle implements IntersectionProcessor {
	private int minimumHitsForStatistics = 10;
	 
	
	/** The optimiser main class, for us to get info during run */
	protected OptimiseMulti optimiser;
	
	/** Plane on which to assess focus of the ray bundle */
	private Plane imagePlane;	
	
	/** All rays are fired at the same wavelength */
	private double wavelength;
	
	/** Initial rays used for evaluation - the same ones are used for each eval, to remove random noise */
	protected double rayStarts[][];
	protected double rayDirs[][];	
	protected int nRays;
	
	/** Where they want the rays to hit the image plane (R,U) */
	private double targetPos[] = null;
		
	/** Sharpness cost weighting for this bundle */
	private double sharpnessWeightX = 1.0, sharpnessWeightY = 1.0; 
	/** Intensity cost weighting for this bundle */	
	private double intensityWeight = 0.3;
	/** Target cost weighting for this bundle */	
	private double targetWeight = 0.0;
		
	/** Our bundle ID */
	private int bundleID;
	
	/** Sums of hits and positions and squares on image plane (R,U) for last evaluation*/
	private int evalNHits;	
	private double lastEvalSumR, lastEvalSumU;
	private double lastEvalSumR2, lastEvalSumU2;
	
	/** Spot statistics for last evaluation */
	private double lastEvalFWHMx, lastEvalFWHMy;		
	private double lastEvalTargetDist; 
	private double lastEvalIntensityLoss;


	private ArrayList<Surface> sequentialSurfaceList;

	private IntersectionProcessor extraISecProc;

	public double getLastEvalFWHM() {	return FastMath.sqrt(lastEvalFWHMx*lastEvalFWHMx + lastEvalFWHMy*lastEvalFWHMy);}
	public double getLastEvalFWHMx() {	return lastEvalFWHMx;	}
	public double getLastEvalFWHMy() {	return lastEvalFWHMy;	}
	public double getLastEvalTargetDist() { return lastEvalTargetDist;	}
	public double getLastEvalIntensityLoss() { return lastEvalIntensityLoss;	}

	public RayBundle(int bundleID) {
		this.bundleID = bundleID;
	}
	
	/** Initialise as rays diverging from a series of points across the given vector
	 * @param rayTarget		Element to fire rays at
	 * @param nPoints		Number of points from which to fire rays
	 * @param centralPoint	The central imaging point
	 * @param pointsDelta	Direction and spacing of points 
	 * @param wavelength	Wavelength of rays
	 * @param nRaysPerPoint			Number of rays from each point			
	 */
	public static RayBundle[] initRaysImagingAlongLine(Element rayTarget, Plane imagePlane, double centralPoint[], double pointsDelta[], 
													int nPoints, double wavelength, int nRaysPerPoint, int bundleID0){
		RayBundle[] rayBundles = new RayBundle[nPoints];		
		for(int iP=0; iP < nPoints; iP++){
			double objectPoint[] = new double[]{
				centralPoint[0] + (iP - (nPoints-1)/2.0) * pointsDelta[0],
				centralPoint[1] + (iP - (nPoints-1)/2.0) * pointsDelta[1],
				centralPoint[2] + (iP - (nPoints-1)/2.0) * pointsDelta[2] 
			};
			rayBundles[iP] = new RayBundle(bundleID0 + iP);
			rayBundles[iP].initRaysImaging(rayTarget, objectPoint, wavelength, nRaysPerPoint);
			rayBundles[iP].setImagePlane(imagePlane);
		}
		return rayBundles;
	}
	
	/** Initialise the evaluation ray tracing using rays diverging from a point 
	 * @param rayTarget		Element to fire rays at
	 * @param objectPoints	Points from which to fire rays - double[posIdx][x/y/z]
	 * @param wavelength	Wavelength of rays
	 * @param nRaysPerPoint			Number of rays from each point			
	 */
	public void initRaysImaging(Element rayTarget, double sourcePoint[], double wavelength, int nRays){
		this.wavelength = wavelength;
		this.nRays = nRays;
		
		rayStarts = new double[nRays][];
		rayDirs = new double[nRays][];
		for(int iR = 0; iR < nRays; iR++){
			rayStarts[iR] = sourcePoint.clone();
			rayDirs[iR] = Tracer.generateRandomRayTowardSurface(sourcePoint, rayTarget);
		}
	}
	
	/** Initialise the evaluation ray tracing using sets of parallel rays at various angles
	 * @param rayTarget		Element to fire rays at.
	 * @param rayDir0		Direction of rays at 0 angle - usually the optic axis.
	 * @param maxAngle		Maximum angle away from rayDir0
	 * @param wavelength	Wavelength of rays
	 * @param nSets			Number of ray bunches, at different angles
	 * @param nRaysPerSet	Number of rays at each angle
	 */
	public static RayBundle[] initRaysParallelMultiAngle(Element rayTarget, Plane imagePlane, double rayDir0[], double maxAngle, double rayLength, 
													double wavelength, int nBundles, int nRaysPerSet, int bundleID0){

		double cosMaxTheta = Math.cos(maxAngle);
		double aU[] = Util.createPerp(rayDir0);
		double bU[] = Util.cross(rayDir0, aU);		
		Util.reNorm(aU);
		Util.reNorm(bU);
		
		RayBundle rayBundles[] = new RayBundle[nBundles];
		for(int iS=0; iS < nBundles; iS++){

			//double cosTheta=1.0;
			//if(nSets > 1) //generate an incoming direction, uniformly distributed over the specified angle range
				//cosTheta = 1 - RandomManager.instance().nextUniform(0, 1) * (1 - cosMaxTheta);
			
			//double theta = -maxAngle + iS * 2 * maxAngle / (nSets - 1);
			double theta = iS * maxAngle / (nBundles - 1);
			 			
			//double sinTheta = FastMath.sqrt(1 - cosTheta*cosTheta);
			double cosTheta = FastMath.cos(theta);
			double sinTheta = FastMath.sin(theta);
			
			double phi = 0;//RandomManager.instance().nextUniform(0, 1) * 2 * Math.PI;
			
			//generate in coord sys (a,b,c) with c as axis toward target 
			double a = sinTheta * FastMath.cos(phi);
			double b = sinTheta * FastMath.sin(phi);
			double c = cosTheta;			
			 
			double rayDir[] = new double[]{
				c * rayDir0[0] + a * aU[0] + b * bU[0],
				c * rayDir0[1] + a * aU[1] + b * bU[1],
				c * rayDir0[2] + a * aU[2] + b * bU[2]
			};
			
			rayBundles[iS] = new RayBundle(bundleID0 + iS);
			rayBundles[iS].initRaysParallel(rayTarget, rayDir, aU, bU, rayLength, wavelength, nRaysPerSet);
			rayBundles[iS].setImagePlane(imagePlane);
		}
		return rayBundles;
	}
	
	/** Initialise the rays at a fixed angle
	 * @param rayTarget		Element to fire rays at.
	 * @param rayDir0		Direction of rays at 0 angle - usually the optic axis.
	 * @param maxAngle		Maximum angle away from rayDir0
	 * @param wavelength	Wavelength of rays
	 * @param nSets			Number of ray bunches, at different angles
	 * @param nRaysPerSet	Number of rays at each angle
	 * @param aU, bU		Unit vectors defining a plane perp to the central ray angle (not necc. this one) (see initRaysParallelMultiAngle() )
	 */
	public void initRaysParallel(Element rayTarget, double rayDir[], double aU[], double bU[], double rayLength, double wavelength, int nRays){
		this.wavelength = wavelength;
		this.nRays = nRays;
		
		rayStarts = new double[nRays][];
		rayDirs = new double[nRays][];
		
		//generate multiple rays in this direction, distributed randomly accross the bouding
		//sphere of the target element
		for(int iR=0; iR < nRays; iR++){
			rayDirs[iR] = rayDir.clone();
			
			//generate a random point on the circle made by the projection of the target bounding sphere
			// on to the perpendicular plane to rayDir0
			
			double r = rayTarget.getBoundarySphereRadius() * FastMath.sqrt(RandomManager.instance().nextUniform(0, 1));
			double ang = RandomManager.instance().nextUniform(0, 1) * 2*Math.PI;
			double lA = r * FastMath.sin(ang);
			double lB = r * FastMath.cos(ang);
			double cc[] = rayTarget.getBoundarySphereCentre();
			
			rayStarts[iR] = new double[]{
					cc[0] + lA * aU[0] + lB * bU[0] - rayLength * rayDir[0],
					cc[1] + lA * aU[1] + lB * bU[1] - rayLength * rayDir[1],
					cc[2] + lA * aU[2] + lB * bU[2] - rayLength * rayDir[2],
			};
		}		
	}	
	
	public void setExtraIntersectionProcessor(IntersectionProcessor extraISecProc){
		this.extraISecProc = extraISecProc;
	}
	
	/** Trace this ray bundle and evaluate the cost function */
	public double traceAndEvaluate(List<RayDrawer> rayDrawers) {
		
		if(imagePlane == null)
			throw new IllegalArgumentException("imagePlane not set for ray bundle '"+bundleID+"'");
		
		evalNHits = 0;
		lastEvalSumR = 0;
		lastEvalSumR2 = 0;
		lastEvalSumU = 0;
		lastEvalSumU2 = 0;
		
		if(rayDrawers != null)
			for(RayDrawer rayDrawer : rayDrawers)
				rayDrawer.startGroup("bundle_" + bundleID);				 
		
		for(int rayNo=0; rayNo < rayStarts.length; rayNo++){ //for each ray
			RaySegment ray = new RaySegment();
			
			ray.startPos = rayStarts[rayNo];
			ray.dir = rayDirs[rayNo];

			ray.E0 = new double[][]{{ 1, 0, 0, 0 }};
			ray.length = Double.POSITIVE_INFINITY;
			ray.up = Util.createPerp(ray.dir);
			
			ray.wavelength = wavelength;
			
			if(sequentialSurfaceList != null){
				Tracer.trace(sequentialSurfaceList, ray, Integer.MAX_VALUE, 1e-3);
			}else{
				Tracer.trace(optimiser.getAllElements(), ray, Integer.MAX_VALUE, 1e-3, false);
			}
			
			int n = ray.processIntersections(imagePlane, this, extraISecProc);
			
			if(n > 0 && rayDrawers != null)
				for(RayDrawer rayDrawer : rayDrawers)
					rayDrawer.drawRay(ray);
			
			Pol.recoverAll();
		}
		
		if(rayDrawers != null)
			for(RayDrawer rayDrawer : rayDrawers)
				rayDrawer.endGroup();
		
		// process stats		
		double meanR = lastEvalSumR / evalNHits;
		double meanU = lastEvalSumU / evalNHits;
		double stdR2 = (lastEvalSumR2 - lastEvalSumR*lastEvalSumR/evalNHits) / evalNHits;
		double stdU2 = (lastEvalSumU2 - lastEvalSumU*lastEvalSumU/evalNHits) / evalNHits;
		
		lastEvalFWHMx = 2.35 * FastMath.sqrt(stdR2);	
		lastEvalFWHMy = 2.35 * FastMath.sqrt(stdU2);
		double lastEvalFWHM = 2.35 * FastMath.sqrt(stdR2 + stdU2); 
		lastEvalIntensityLoss = ((double)nRays - evalNHits) / nRays;
		if(targetWeight != 0){
			if(targetPos == null)
				throw new RuntimeException("No targetPos but targetWeight is nonzero. Did you forget to call findTarget() or setTargetPos()?");
			lastEvalTargetDist = FastMath.sqrt(
										(Double.isNaN(targetPos[0]) ? 0 : FastMath.pow2(meanR - targetPos[0])) + 
										(Double.isNaN(targetPos[1]) ? 0 : FastMath.pow2(meanU - targetPos[1])));
		}
				
		double cost = sharpnessWeightX * lastEvalFWHMx 
						+ sharpnessWeightY * lastEvalFWHMy 
						+ targetWeight * lastEvalTargetDist 
						+ intensityWeight * lastEvalIntensityLoss;
				
		if(evalNHits <= minimumHitsForStatistics){
			System.err.println("WARNING: Ray bundle "+bundleID+" had too few hits ("+evalNHits+"), not adding to statistics.");
			cost =  ((evalNHits == 0) ? 10000 : cost) + 100.0 * FastMath.exp(- evalNHits / minimumHitsForStatistics);
		}
		
		if(Double.isNaN(cost)){
			System.out.println("RayBundle.traceAndEvaluate(): NaN after checks");
		}
		
		return cost;
	}
	
	
	/** Traces this ray bundle to the image plane and sets the target to the average position. */
	public void findTarget() {
		
		double targetWeightSave = targetWeight;
		targetWeight = 0;
		
		/*ArrayList<RayDrawer> drawers  = new ArrayList<RayDrawer>();
		VRMLDrawer drawer = new VRMLDrawer("/tmp/wtf.vrml");
		drawers.add(drawer);
		//drawer.drawElement(optimiser.getAllElements());
		traceAndEvaluate(drawers);
		drawer.destroy();*/
		traceAndEvaluate(null);
	
		if(evalNHits <= minimumHitsForStatistics)
			throw new RuntimeException("Ray bundle "+bundleID+" had too few hits ("+evalNHits+") while trying to find initial target.");			
		
		// process stats		
		double meanR = lastEvalSumR / evalNHits;
		double meanU = lastEvalSumU / evalNHits;
		
		targetPos = new double[]{ meanR, meanU };
		targetWeight = targetWeightSave;
		
	}

	@Override
	public void nextIntersection(Intersection hit) {
		
		double posRU[] = imagePlane.posXYZToPlaneRU(hit.pos);
		
		
		BinaryMatrixWriter hitsOut = optimiser.getHitsOutWriter();
		if(hitsOut != null){
			hitsOut.writeRow(optimiser.getIterationNo(), bundleID, evalNHits, posRU);
		}
		
		evalNHits++;
		lastEvalSumR += posRU[0];
		lastEvalSumU += posRU[1];
		lastEvalSumR2 += posRU[0]*posRU[0];
		lastEvalSumU2 += posRU[1]*posRU[1];
	}
	
	
	@Override
	public String toString() {
		
		double meanR = lastEvalSumR / evalNHits;
		double meanU = lastEvalSumU / evalNHits;
		
		return "RayBundle[ID="+bundleID+
					"]: start[0]={" + rayStarts[0][0] + ", " +rayStarts[0][1] + ", " +rayStarts[0][2] + 
					" }, lastEval: n="+evalNHits + ", mean={" + meanR + ", " + meanU + "}, fwhm = " + getLastEvalFWHM() + ", Iloss = " + lastEvalIntensityLoss + ", targDist = " + lastEvalTargetDist;   
	}
	
	/** Sets the sharpness weight. 
	 * The FWHM of the spot is multiplied by this and added to the cost */
	public void setSharpnessWeight(double sharpnessWeight) { this.sharpnessWeightX = sharpnessWeight; this.sharpnessWeightY = sharpnessWeight;	}
	public void setSharpnessWeights(double sharpnessWeightX, double sharpnessWeightY) { this.sharpnessWeightX = sharpnessWeightX; this.sharpnessWeightY = sharpnessWeightY;	}
	
	/** Sets the intensity weight. 
	 * The number of rays that do NOT hit the target is multiplied by this and added to the cost */
	public void setIntensityWeight(double intensityWeight) { this.intensityWeight = intensityWeight;	}
	
	/** Sets the target weight. 
	 * The distance of the mean hit from the target position is multiplied by this and added to the cost */
	public void setTargetWeight(double targetWeight) { this.targetWeight = targetWeight;	}
	
	/** Sets the plane on which to test the imaging/intesity/target */
	public void setImagePlane(Plane imagePlane) { this.imagePlane = imagePlane; }
	
	/** Sets the target position on the image plane (R,U) */
	public void setTargetPos(double posRU[]){ 
		this.targetPos = posRU;
	}

	/** Returns the source position of the first ray */
	public double[] getSourcePos() { return rayStarts[0]; }
	
	public void setOptimiser(OptimiseMulti optim) { this.optimiser = optim; }

	public void setMinimumHitsForStatistics(int minimumHitsForStatistics) {
		this.minimumHitsForStatistics = minimumHitsForStatistics;
	}
	public void setSequentialSurfaceList(ArrayList<Surface> sequentialSurfaceList) {
		this.sequentialSurfaceList = sequentialSurfaceList;
				
	}
	
}
