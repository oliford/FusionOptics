package fusionOptics.optimisationMulti;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import javax.print.attribute.standard.Media;

import oneLiners.OneLiners;
import otherSupport.ColorMaps;

import binaryMatrixFile.BinaryMatrixWriter;

import seed.digeom.FunctionND;
import seed.digeom.IDomain;
import seed.digeom.RectangularDomain;
import seed.optimization.HookeAndJeeves;
import seed.optimization.Optimizer;

import fusionOptics.Util;
import fusionOptics.collection.IntersectionProcessor;
import fusionOptics.drawing.RayDrawer;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.optimisation.OptimiseOptic;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Plane;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Material;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;

/** New optimisation system for optimising lots of parameters of lots of optics against lots of costs. */
public class OptimiseMulti extends FunctionND {

	private ArrayList<Parameter> parameters = new ArrayList<Parameter>();
	private ArrayList<RayBundle> rayBundles = new ArrayList<RayBundle>();
	private Element allElements;
	
	
	/** Debugging output */
	private String outputPrefix = null;
	private boolean hitsWriteEnable = false;
	private BinaryMatrixWriter hitsOutWriter = null;	
	private BinaryMatrixWriter progressWriter = null;
	/** Current iteration number */
	protected int iterationNo = -1;
	/** Current eval params */
	protected double currentEvalParams[];
	
	/** Graphics output */
	private class IterationDrawings { 
		public int skipRays;
		public String fileNamePrefix;
		public Class type = null;
		public String extraVRMLStart = null;
		public String extraVRMLEnd = null;
		public RayDrawer rayDrawer;
	}
	private ArrayList<IterationDrawings> drawings = new ArrayList<OptimiseMulti.IterationDrawings>();
	private int outputIterationPeriod;
	
	/** The optimiser */
	private Optimizer opt;
	
	/** Number of rays that hit their image planes on the last evaluation */
	public int nRaysInImageLast;
	
	private double lastBestFWHM, lastBestTargetDist, lastBestIntensityLoss;
	
	public OptimiseMulti() {
		opt = new HookeAndJeeves(this); //default, not great at more than a few dimensions
		
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
	}
	
	/** Sets the MegaOptimisation optimiser (seed.optimiser...) */
	public void setOptimiser(Optimizer opt) {
		this.opt = opt;
	}
	
	public double[] optimise(int nIters){
		return optimise(getParams(), nIters);
	}
	
	public double[] optimise(double pInit[], int nIters){
		if(allElements == null) throw new IllegalArgumentException("No elements to trace, please call OptimiseMulti.setTracingElements() first");
		if(rayBundles.size() <= 0) throw new IllegalArgumentException("No rays, please call OptimiseMulti.addRayBundle() at least once ");
		if(parameters.size() <= 0) throw new IllegalArgumentException("No parameters, please call OptimiseMulti.addParameter() at least once ");
		
		if(opt == null){
			opt = new HookeAndJeeves(this); //simple default, otherwise use some of the following from outside:
		}
		
		opt.setObjectiveFunction(this);
		
		opt.init(pInit);
		
		//do the optimiser init before the SVG setup, so it doesnt draw the init exploration
		if(outputPrefix != null){
			hitsOutWriter = new BinaryMatrixWriter(outputPrefix + "-hits.bin", 5);
			progressWriter = new BinaryMatrixWriter(outputPrefix + "-progress.bin", pInit.length + 5);
		}
		
		//iterationNo = -1;
		//hitsWriteEnable = true;
		double initCost = eval(pInit);
		//hitsWriteEnable = false; //only enable hits writing on evals we do out here, not that the optimiser does, other we get all the 'attempts', no matter how bad
		
		double t0 = System.currentTimeMillis();
		for(iterationNo=-1; iterationNo < nIters; iterationNo++){
				
			if(iterationNo >= 0)
				opt.refine();
			double p[] = opt.getCurrentPos();
			
			//re-run the eval, because the optimisers don't guarantee the optimal was eval'd last
			double cost;
			
			if(outputIterationPeriod > 0 && ((iterationNo+1) % outputIterationPeriod) == 0){
				hitsWriteEnable = true;
				
				cost = evalAndDraw(p);
				for(Parameter param : parameters)
					System.out.println(param.toString());
				
				hitsWriteEnable = false;
			}else
				cost = eval(p); 
				
			double itsPerSec = 1000.0 * (iterationNo+1) / (System.currentTimeMillis() - t0);
			System.out.print("i="+ iterationNo + "\tcost=" + cost + "\tnRaysInImageLast=" + nRaysInImageLast + "\titsPerSec=" + itsPerSec + "\tp=");
			for(int j=0; j < p.length; j++){
				System.out.print("\t" + p[j]);
			}
			System.out.println();
			
			dumpRayStats();	
			opt.dumpStatus();
			
			if(progressWriter != null){
				progressWriter.writeRow(iterationNo, nRaysInImageLast, lastBestFWHM, lastBestTargetDist, lastBestIntensityLoss, p);
			}

		}
		
		double p[] = opt.getCurrentPos();
		setParams(p);
		
		System.out.print("Init: cost=" + initCost + ", p = "); OneLiners.dumpArray(pInit);
		System.out.print("End:  cost=" + eval(p) + ", p = "); OneLiners.dumpArray(p);

		return p;
	}

	private void dumpRayStats() {
		double sumFWHM=0, sumTargDist=0, sumIntensityLoss=0;
		int nBundles = rayBundles.size();
		for(RayBundle rayBundle : rayBundles){
			sumFWHM += rayBundle.getLastEvalFWHM();
			sumTargDist += rayBundle.getLastEvalTargetDist();
			sumIntensityLoss += rayBundle.getLastEvalIntensityLoss();
		}
		lastBestFWHM = (sumFWHM/nBundles);
		lastBestTargetDist = (sumTargDist/nBundles);
		lastBestIntensityLoss = (sumIntensityLoss/nBundles);
		System.out.println("i=" + iterationNo 
				+ " <fwhm>=" + lastBestFWHM
				+ " <dist>=" + lastBestTargetDist
				+ " <Iloss>=" + lastBestIntensityLoss);

		
	}

	private void setParams(double[] p) {
		int n = parameters.size();
		if(p.length != n)
			throw new IllegalArgumentException("Params are the wrong length");
			
		int i=0;
		for(Parameter parameter : parameters){
			parameter.set(p[i]);
			i++;
		}
	}

	private double[] getParams() {
		int n = parameters.size();
		double p[] = new double[n];
			
		int i=0;
		for(Parameter parameter : parameters){
			p[i] = parameter.get();
			i++;
		}
		return p;
	}

	@Override
	public double eval(double[] p) {
		currentEvalParams = p;
		setParams(p);
		
		double cost = 0;		
		for(RayBundle rayBundle : rayBundles){
			
			cost += rayBundle.traceAndEvaluate(null);
		}
		
		return cost;
	}
	
	public double evalAndDraw(double p[]){
		
		double colorMap[][] = ColorMaps.jet(rayBundles.size());
		ArrayList<RayDrawer> drawers = new ArrayList<RayDrawer>();
		for(IterationDrawings drawing : drawings){
			
			if(drawing.type == SVGRayDrawing.class){
				drawing.rayDrawer = new SVGRayDrawing(drawing.fileNamePrefix + "-iter_" + iterationNo + ".svg", Util.getBoundingBox(allElements), true);			
				
			}else if(drawing.type == VRMLDrawer.class){
				drawing.rayDrawer = new VRMLDrawer(drawing.fileNamePrefix + "-iter_" + iterationNo + ".vrml");
				if(drawing.extraVRMLStart != null){
					((VRMLDrawer)drawing.rayDrawer).addVRML(drawing.extraVRMLStart);
				}
			}else 
				throw new IllegalArgumentException("Unknown drawing type");
			drawing.rayDrawer.setSkipRays(drawing.skipRays);
			
			if(drawing.rayDrawer instanceof SVGRayDrawing){
				((SVGRayDrawing)drawing.rayDrawer).generateLineStyles(colorMap, allElements.getBoundarySphereRadius() / 5000);
				((SVGRayDrawing)drawing.rayDrawer).setRayColour(0);
			}
		
			drawing.rayDrawer.drawElement(allElements);
			
			drawers.add(drawing.rayDrawer);
		}

		double cost = 0;		
		int bundleNum = 0;
		for(RayBundle rayBundle : rayBundles){
			
			int colNum = bundleNum;
			
			for(RayDrawer rayDrawer : drawers){
				if(rayDrawer instanceof SVGRayDrawing){					
					((SVGRayDrawing)rayDrawer).setRayColour(colNum);
			
				}else if(rayDrawer instanceof VRMLDrawer){
					((VRMLDrawer)rayDrawer).setRayColour(colorMap[colNum]);
					
				}
			}
			
			bundleNum++;
			
			cost += rayBundle.traceAndEvaluate(drawers);
		}
	
		for(IterationDrawings drawing : drawings){
			if(drawing.extraVRMLEnd != null){
				((VRMLDrawer)drawing.rayDrawer).addVRML(drawing.extraVRMLEnd);
			}
			
			drawing.rayDrawer.destroy();
			drawing.rayDrawer = null;
		}

		return cost;
	}
	
	@Override
	public IDomain getDomain() { 
		int n = parameters.size();
		double min[] = new double[n];
		double max[] = new double[n];
			
		int i=0;
		for(Parameter parameter : parameters){
			min[i] = parameter.min();
			max[i] = parameter.max();
			i++;
		}
		return new RectangularDomain(min, max);
	}
	
	/** If set, an svg files are created containing the tracing for all evaluations */  
	public void addRegularDrawing(String fileNamePrefix, Class type, int nSkipRays){
		addRegularDrawing(fileNamePrefix, type, nSkipRays, null, null);
	}
	
	/** If set, an svg files are created containing the tracing for all evaluations */  
	public void addRegularDrawing(String fileNamePrefix, Class type, int nSkipRays, String extraVRMLStart, String extraVRMLEnd){
		
		IterationDrawings drawing = new IterationDrawings();
		drawing.fileNamePrefix = fileNamePrefix;
		drawing.type = type;
		drawing.skipRays = nSkipRays;
		drawing.extraVRMLStart = extraVRMLStart;
		drawing.extraVRMLEnd = extraVRMLEnd;
		drawings.add(drawing);
	}
	
	public void setOutputIterationPeriod(int outputIterationPeriod) {
		this.outputIterationPeriod = outputIterationPeriod;
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
	 * -hits: { iteration No, ray bundle No, ray No, imgX, imgY }
	 * -progress: { iterationNo, nRaysInImageLast, lastBestFWHM, lastBestTargetDist, lastBestIntensityLoss, {params}}
	 * @param hitsOutputFile		Filename, or null to not write debugging info	 */
	public void setOutputPrefix(String outputPrefix){ this.outputPrefix = outputPrefix;	}
	public BinaryMatrixWriter getHitsOutWriter() { return hitsWriteEnable ? hitsOutWriter : null; }
	public BinaryMatrixWriter getProgressWriter() { return progressWriter; }
	
	public int getIterationNo() { return iterationNo; }
	public Element getAllElements() { return allElements; }

	public void addParameter(Parameter parameter){ parameters.add(parameter);	}
	public void addParameters(Collection<Parameter> parameter){ parameters.addAll(parameter);	}
	public List<Parameter> getParameterList(){ return parameters; }
	
	public void addRayBundle(RayBundle rayBundle){ rayBundle.optimiser = this; rayBundles.add(rayBundle);	}
	public void setTracingElements(Element allElements) { this.allElements = allElements;	}
	public void addRayBundles(RayBundle rayBundleArray[]){ 
		for(RayBundle rayBundle : rayBundleArray){
			rayBundle.optimiser = this;
			rayBundles.add(rayBundle);
		}
	}

	public void dumpParams() {
		for(Parameter parameter : parameters){
			System.out.println(parameter.toString());
		}
	}
	
	public void dumpRayBundles() {
		for(RayBundle rayBundle : rayBundles){
			System.out.println(rayBundle.toString());
		}
	}

	public Object getCurrentEvalParams() { return currentEvalParams;	}

	public double eval() { return eval(getParams()); }

	public double getCurrentBestFWHM() { return lastBestFWHM; 	}
	public double getCurrentBestTargetDist() { return lastBestTargetDist; 	}
	public double getCurrentBestIntensityLoss() { return lastBestIntensityLoss; 	}

	public List<RayBundle> getRayBundles() { return rayBundles;	}
	
	public void destroy(){
		if(hitsOutWriter != null)
			hitsOutWriter.close();
		if(progressWriter != null)
			progressWriter.close();		
	}
	
}
