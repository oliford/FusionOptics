package fusionOptics.optimisation;

import java.text.DecimalFormat;

import fusionOptics.Util;
import fusionOptics.collection.IntersectionProcessor;
import fusionOptics.drawing.RayDrawer;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.surfaces.Plane;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;

import net.jafama.FastMath;
import algorithmrepository.Algorithms;
import algorithmrepository.ConjugateGradient;
import seed.digeom.Function;
import seed.digeom.FunctionND;
import seed.digeom.IDomain;
import seed.digeom.IFunction;
import seed.digeom.RectangularDomain;
import seed.optimization.BracketingByParameterSpace;
import seed.optimization.ConjugateGradientDirectionFR;
import seed.optimization.CoordinateDescentDirection;
import seed.optimization.GoldenSection;
import seed.optimization.GradientDescentDirection;
import seed.optimization.HookeAndJeeves;
import seed.optimization.IStoppingCondition;
import seed.optimization.LineSearchOptimizer;
import seed.optimization.MaxIterCondition;
import seed.optimization.Optimizer;
import seed.optimization.SearchDirectionMethod;
import seed.optimization.genetic.Genetic;
import seed.optimization.genetic.MSHGPPConfig;
import seed.optimization.genetic.MegaSuperHyperGeneticProblemPacifier;
import seed.optimization.genetic.SuperGeneticMk2;
import uk.co.oliford.jolu.BinaryMatrixWriter;
import uk.co.oliford.jolu.ColorMaps;
import uk.co.oliford.jolu.OneLiners;
import uk.co.oliford.jolu.RandomManager;

/** Optimises properties of an optical component for the best focus (smallest PSF)
 * 
 * Abstract - Implementors decide what to optimise
 * 
 * After construction, call:
 *  OptimseOptic.setTracingElements()
 *	OptimseOptic.initRaysImaging or OptimseOptic.initRaysImaging()
 *  OptimseOptic.setParameterSpace()		
 * 
 * @deprecated Replaced by a newer generalised optimiser: see fusionOptics.optimisationMulti
 * @author oliford
 */
public abstract class OptimiseOptic extends FunctionND implements IntersectionProcessor {

	/* --- Setup --- */
	/** Plane on which to assess focus */
	private Plane imagePlane;	
	/** Minimum and maximum of parameter space */
	private double pMin[], pMax[];		
	/** All rays are fired at the same wavelength */
	private double wavelength;	
	/** Value add to RMS for rays that don't hit the image plane */
	private double distForOut;	
	/** Initial rays used for evaluation - the same ones are used for each eval, to remove random noise */  
	protected double rayStarts[][][];
	protected double rayDirs[][][];	
	protected int nRaysTotal;
	/** Elements to be ray traced */
	private Element all;	
	/** Debugging output */
	private BinaryMatrixWriter hitsOut = null;	
	/** Graphics output */
	private RayDrawer rayDrawer = null;	
	private boolean enableDrawing = false;
	
	/* --- State --- */
	/** Current evaluation number */
	protected int evalNo;	
	/*** Current ray bundle number */
	protected int setNo;
	/** Current ray number */
	private int rayNo;		
	/** Number of rays that hit the imaging plane on the last evaluation */
	public int nRaysInImageLast;
	/** Colelcted hit info for the stats */
	protected double hits[][];
	/** FWHM of last evaluation */
	private double lastFWHM;
	
	private String hitsOutputFile = null;
	
	/** Sets the image plane and the elements to trace 
	 * 
	 * @param all		All elements that might be involved in the tracing
	 * @param imagePlane	Plane on which to assess focus.
	 */
	public void setTracingElements(Element all, Plane imagePlane) {
		this.all = all;
		this.imagePlane = imagePlane;		
		this.distForOut = imagePlane.getBoundarySphereRadius()*1.1;		
	}
	
	/** Initialise the evaluation ray tracing using rays diverging from a series of points across the given vector
	 * @param rayTarget		Element to fire rays at
	 * @param nPoints		Number of points from which to fire rays
	 * @param centralPoint	The central imaging point
	 * @param pointsDelta	Direction and spacing of points 
	 * @param wavelength	Wavelength of rays
	 * @param nRaysPerPoint			Number of rays from each point			
	 */
	public void initRaysImaging(Element rayTarget, double centralPoint[], double pointsDelta[], int nPoints, double wavelength, int nRaysPerPoint){
		double objectPoints[][] = new double[nPoints][3];
		for(int iP=0; iP < nPoints; iP++){
			objectPoints[iP][0] = centralPoint[0] + (iP - nPoints/2.0) * pointsDelta[0];
			objectPoints[iP][1] = centralPoint[1] + (iP - nPoints/2.0) * pointsDelta[1];
			objectPoints[iP][2] = centralPoint[2] + (iP - nPoints/2.0) * pointsDelta[2]; 
		}
		initRaysImaging(rayTarget, objectPoints, wavelength, nRaysPerPoint);
	}
	
	/** Initialise the evaluation ray tracing using rays diverging from from some points 
	 * @param rayTarget		Element to fire rays at
	 * @param objectPoints	Points from which to fire rays - double[posIdx][x/y/z]
	 * @param wavelength	Wavelength of rays
	 * @param nRaysPerPoint			Number of rays from each point			
	 */
	public void initRaysImaging(Element rayTarget, double objectPoints[][], double wavelength, int nRaysPerPoint){
		this.wavelength = wavelength;
		
		int nPoints = objectPoints.length;
		
		rayStarts = new double[nPoints][nRaysPerPoint][];
		rayDirs = new double[nPoints][nRaysPerPoint][];
		for(int iP=0; iP < nPoints; iP++){
			for(int iR = 0; iR < nRaysPerPoint; iR++){
				rayStarts[iP][iR] = objectPoints[iP].clone();
				rayDirs[iP][iR] = Tracer.generateRandomRayTowardSurface(objectPoints[iP], rayTarget);
			}
		}
		nRaysTotal = nPoints * nRaysPerPoint;
	}
	
	/** Initialise the evaluation ray tracing using sets of parallel rays at various angles
	 * @param rayTarget		Element to fire rays at.
	 * @param rayDir0		Direction of rays at 0 angle - usually the optic axis.
	 * @param maxAngle		Maximum angle away from rayDir0
	 * @param wavelength	Wavelength of rays
	 * @param nSets			Number of ray bunches, at different angles
	 * @param nRaysPerSet	Number of rays at each angle
	 */
	public void initRaysParallel(Element rayTarget, double rayDir0[], double maxAngle, double rayLength, 
									double wavelength, int nSets, int nRaysPerSet){
		this.wavelength = wavelength;
		
		double cosMaxTheta = Math.cos(maxAngle);
		double aU[] = Util.createPerp(rayDir0);
		double bU[] = Util.cross(rayDir0, aU);		
		Util.reNorm(aU);
		Util.reNorm(bU);
		
		rayStarts = new double[nSets][nRaysPerSet][];
		rayDirs = new double[nSets][nRaysPerSet][];
		for(int iS=0; iS < nSets; iS++){

			//double cosTheta=1.0;
			//if(nSets > 1) //generate an incoming direction, uniformly distributed over the specified angle range
				//cosTheta = 1 - RandomManager.instance().nextUniform(0, 1) * (1 - cosMaxTheta);
			
			//double theta = -maxAngle + iS * 2 * maxAngle / (nSets - 1);
			double theta = iS * maxAngle / (nSets - 1);
			 			
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
			
			//generate multiple rays in this direction, distributed randomly accross the bouding
			//sphere of the target element
			for(int iR=0; iR < nRaysPerSet; iR++){
				rayDirs[iS][iR] = rayDir.clone();
				
				//generate a random point on the circle made by the projection of the target bounding sphere
				// on to the perpendicular plane to rayDir0
				
				double r = rayTarget.getBoundarySphereRadius() * FastMath.sqrt(RandomManager.instance().nextUniform(0, 1));
				double ang = RandomManager.instance().nextUniform(0, 1) * 2*Math.PI;
				double lA = r * FastMath.sin(ang);
				double lB = r * FastMath.cos(ang);
				double cc[] = rayTarget.getBoundarySphereCentre();
				
				rayStarts[iS][iR] = new double[]{
						cc[0] + lA * aU[0] + lB * bU[0] - rayLength * rayDir[0],
						cc[1] + lA * aU[1] + lB * bU[1] - rayLength * rayDir[1],
						cc[2] + lA * aU[2] + lB * bU[2] - rayLength * rayDir[2],
				};
			}
		}
		nRaysTotal = nRaysPerSet * nSets;
	}	
	
	/** Sets the parameter space searched 
	 * 
	 * @param pMin
	 * @param pMax
	 */
	public void setParameterSpace(double pMin[], double pMax[]){
		this.pMin = pMin;
		this.pMax = pMax;		
	}
	
	public double[] optimise(Optimizer opt, double pInit[], int nIters){
		if(all == null) throw new IllegalArgumentException("No elements to trace, please call OptimseOptic.setTracingElements() first");
		if(rayStarts == null) throw new IllegalArgumentException("No rays, please call OptimseOptic.initRaysXXX() first");
		if(pMin == null) throw new IllegalArgumentException("No parameter space, please call OptimseOptic.setParameterSpace() first");
		
		if(opt == null){
			opt = new HookeAndJeeves(this); //simple default, otherwise use some of the following from outside:
		}
		//HookeAndJeeves opt = new HookeAndJeeves(this);
		//SuperGeneticMk2 opt = new SuperGeneticMk2(4, 1, true);
		
		/*MSHGPPConfig cfg = new MSHGPPConfig();
		cfg.populationSize = 20;
		cfg.numChildren = 2;
		cfg.initFromCurrent = true;
		MegaSuperHyperGeneticProblemPacifier opt = new MegaSuperHyperGeneticProblemPacifier(cfg);
		//*/
		
		/*IStoppingCondition lineSearchStop = new MaxIterCondition(50);
		GoldenSection lineOptimizer = new GoldenSection(lineSearchStop);
		lineOptimizer.setInitialBracketMethod(new BracketingByParameterSpace());
		Optimizer opt = new LineSearchOptimizer(null, new ConjugateGradientDirectionFR(), lineOptimizer);
		//*/
		opt.setObjectiveFunction(this);
		
		opt.init(pInit);
		
		//do the optimiser init before the SVG setup, so it doesnt draw the init exploration
		if(hitsOutputFile != null)
			hitsOut = new BinaryMatrixWriter(hitsOutputFile, pInit.length + 5);
		
		if(rayDrawer != null){
			if(rayDrawer instanceof SVGRayDrawing){
				//((SVGRayDrawing)rayDrawer).setSkipRays((rayStarts.length * rayStarts[0].length) / 100);
				((SVGRayDrawing)rayDrawer).generateLineStyles(ColorMaps.jet((nIters+1)*rayStarts.length), all.getBoundarySphereRadius() / 5000);
				((SVGRayDrawing)rayDrawer).setRayColour(0);
			}
			
			rayDrawer.drawElement(all);
			enableDrawing=true;
		}
				
		evalNo = 0;
		double initCost = eval(pInit);
		double initFWHM = lastFWHM;
		
		for(int i=0; i < nIters; i++){
			if(rayDrawer != null && rayDrawer instanceof SVGRayDrawing)
				((SVGRayDrawing)rayDrawer).setRayColour((i+1)*rayStarts.length+1);
				
			opt.refine();
			double p[] = opt.getCurrentPos(); 
			//cost is FWHM			
			eval(p); //re-run the eval, because the optimisers don't guarantee the optimal was eval'd last
			System.out.print("i="+ i + "\tFWHM=" + lastFWHM + "\tnRaysInImageLast=" + nRaysInImageLast + "\tp=");
			for(int j=0; j < p.length; j++){
				System.out.print("\t" + p[j]);
			}
			System.out.println();
			opt.dumpStatus();
		}
		hits = null;
		
		double p[] = opt.getCurrentPos();
		setParams(p);
		
		System.out.print("Init: FWHM=" + initFWHM + ", p = "); OneLiners.dumpArray(pInit);
		System.out.print("End:  FWHM=" + eval(p) + ", p = "); OneLiners.dumpArray(p);
		
		if(rayDrawer != null){
			rayDrawer.destroy();
			rayDrawer = null;
		}
		enableDrawing=false;
		return p;
	}

	@Override
	/** Default cost function evaluation - Spot size variance for each ray bundle / source point */
	public double eval(double p[]) {
		setParams(p);
		
		double varSumAllSets = 0; //variance for all sets
		nRaysInImageLast = 0;
		
		if(enableDrawing && rayDrawer != null)
			rayDrawer.startGroup("eval_" + evalNo);
		
		for(setNo=0; setNo < rayStarts.length; setNo++){  // for each ray bundle
			if(enableDrawing && rayDrawer != null)
				rayDrawer.startGroup("set_" + evalNo);
			
			if(enableDrawing && rayDrawer instanceof SVGRayDrawing)
				((SVGRayDrawing)rayDrawer).setRayColour(setNo);
			
			int nRaysInSet = rayStarts[setNo].length; 
			hits = new double[nRaysInSet][];
			
			for(rayNo=0; rayNo < rayStarts[setNo].length; rayNo++){ //for each ray
				RaySegment ray = new RaySegment();
				
				ray.startPos = rayStarts[setNo][rayNo];
				ray.dir = rayDirs[setNo][rayNo];
	
				ray.E0 = new double[][]{{ 1, 0, 0, 0 }};
				ray.length = Double.POSITIVE_INFINITY;
				ray.up = Util.createPerp(ray.dir);
				
				ray.wavelength = wavelength;
				
				Tracer.trace(all, ray, Integer.MAX_VALUE, 1e-3, false);
				
				if(enableDrawing && rayDrawer != null)
					rayDrawer.drawRay(ray);
				
				int n = ray.processIntersections(imagePlane, this);
				
				if(n < 1){
					hits[rayNo] = null; 
				}
				
				Pol.recoverAll();
			}
			
			if(enableDrawing && rayDrawer != null)
				rayDrawer.endGroup();
		
			varSumAllSets += calcSetSpotSizeVariance(nRaysInSet, p);
		
			calcSetExtraCostFunction(nRaysInSet, p);
		}
		
		varSumAllSets /= nRaysInImageLast;
		
		if(enableDrawing && rayDrawer != null)
			rayDrawer.endGroup();
		
		evalNo++;
		
		lastFWHM = 2.35 * FastMath.sqrt(varSumAllSets);		
		return lastFWHM + extraCostFunction(); 
	}
	
	protected void calcSetExtraCostFunction(int nRaysInSet, double p[]){ }
	protected double extraCostFunction(){ return 0; }
	
	/** calc intended image position for the current set */
	protected double[] targetImagePos(){ return null; }
	
	protected double calcSetSpotSizeVariance(int nRaysInSet, double p[]) {

		double targetPos[] = targetImagePos();
		
		// Calc mean of this set's PSF
		double meanX=0, meanY=0;
		int nRaysHitInThisSet = 0;
		for(int j=0; j < nRaysInSet; j++){
			if(hits[j] != null){
				meanX += hits[j][0];
				meanY += hits[j][1];
				if(hitsOut != null)
					hitsOut.writeRow(evalNo, setNo, j, hits[j], p);
				nRaysHitInThisSet++;
			}
		}	

		if(targetPos != null){
			meanX = targetPos[0];
			meanY = targetPos[1];
		}else{
			meanX /= nRaysHitInThisSet; 
			meanY /= nRaysHitInThisSet;
		}
		
		if(nRaysHitInThisSet <= 10){
			System.err.println("WARNING: Ray set "+setNo+" had too few hits ("+nRaysHitInThisSet+"), not adding to statistics.");
			return 0;
		}else{
		
			// Calc variance of this set's PSF
			double setVarianceSum=0;
			for(int j=0; j < nRaysInSet; j++){
				if(hits[j] != null){
					setVarianceSum += FastMath.pow2(hits[j][0] - meanX) + FastMath.pow2(hits[j][1] - meanY);
				}else{
					//now ignoring points outside image
					//v += distForOut * distForOut;
				}
			}
			//setVariance /= nRaysHitInThisSet;
			nRaysInImageLast += nRaysHitInThisSet;
			return setVarianceSum;
		}
		
	}
	
	@Override
	public void nextIntersection(Intersection hit) {
		hits[rayNo] = imagePlane.posXYZToPlaneUR(hit.pos);
	}
	
	/** Set the parameters, called by base for implementors to update the parameters */
	protected abstract void setParams(double p[]);
	
	/** Get the current parameters, called by base for implementors to return the parameters */
	protected abstract double[] getParams();
	
	@Override
	public IDomain getDomain() { return new RectangularDomain(pMin, pMax);	}
	
	/** If set, an svg files are created containing the tracing for all evaluations */  
	public void setSVGOut(String fileNamePrefix, int nSkipRays){
		if(fileNamePrefix == null){
			rayDrawer = null;
		}else{
			rayDrawer = new SVGRayDrawing(fileNamePrefix, Util.getBoundingBox(all), true);
			rayDrawer.setSkipRays(nSkipRays);
		}
		
		enableDrawing=false; //dont start drawing until optimise() has finished it's init
	}
	

	public void singleParamScan(int pIdx, double p0, double p1, int nSteps) {
		double p[] = getParams();
		
		double v0 = eval(p);
		double dp = (p1 - p0) / (nSteps - 1);
		System.out.println("v0=\t" + v0);
		for(int i=0; i < nSteps; i++){
			p[pIdx] = p0 + i * dp;
			double v = eval(p);
			System.out.println(i + "\t" + p[pIdx] + "\t" + v);
		}
		
	}


	/** Sets fileName to write the optimisation progress to, the file format will be:
	 * { evaluation No, ray set No, ray No, imgX, imgY, {params} }
	 * @param hitsOutputFile		Filename, or null to not write debugging info	 */
	public void setHitsOutputFile(String hitsOutputFile){ this.hitsOutputFile = hitsOutputFile;	}
}
