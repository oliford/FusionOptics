package fusionOptics.optimisationMulti;

import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.surfaces.Plane;
import fusionOptics.types.Element;
import fusionOptics.types.Optic;

public class AutoFocus {
	private Element all, thingToMove;
	private Plane imagePlane;
	
	private int nRaysPerBundle = 100;
	private double wavelen;
	private double maxMovement;
	private double movementDir[];
	private RayBundle rayBundles[];
	
	public AutoFocus(Element all, Plane imagePlane, double wavelen) {
		this.all = all;
		this.imagePlane = imagePlane;
		this.wavelen = wavelen;
	}

	public void setMoveableElement(Element thingtoMove, double dir[], double maxMovement){
		this.thingToMove = thingtoMove;
		this.maxMovement = maxMovement;
		this.movementDir = dir;
	}

	public void setMoveableElement(Plane thingtoMove, double maxMovement){
		this.thingToMove = thingtoMove;
		this.maxMovement = maxMovement;
		this.movementDir = thingtoMove.getNormal();
	}
	
	public void setParallelRays(Element target, double dir[], double maxAngle, int nAngles){
		
		rayBundles = RayBundle.initRaysParallelMultiAngle(target, imagePlane, dir, 
										maxAngle, 0.100, wavelen, nAngles, nRaysPerBundle, 0);
		
	}
	
	public void setImagingRays(Plane objectPlane, Element target, double maxDist, int nPoints){
		double pointsDelta[] = Util.mul(objectPlane.getUp(), (maxDist*2)/(nPoints-1));
		rayBundles = RayBundle.initRaysImagingAlongLine(target, imagePlane, 
				objectPlane.getCentre(), pointsDelta, nPoints, 
				wavelen, nRaysPerBundle, 0);
	}
	
	public void go(int maxInterations) {
		
		//double targetMagnification = 5.0 / 3.0;
		for(RayBundle rayBundle : rayBundles){
			rayBundle.setIntensityWeight(0.0);
			rayBundle.setTargetWeight(0.0);
			rayBundle.setSharpnessWeight(1.0);			
		}		
		
		OptimiseMulti optim = new OptimiseMulti();
		optim.addRayBundles(rayBundles);
		optim.setTracingElements(all);
		
		optim.addParameter(new MoveableElement(imagePlane, new double[]{1,0,0}, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, maxMovement/100));
		

		//optim.setOutputPrefix(outPath + "/autoFocus");
		//optim.addRegularDrawing(outPath + "/autoFocus", SVGRayDrawing.class, 10);
		//optim.addRegularDrawing(outPath + "/autoFocus", VRMLDrawer.class, 10);
		//optim.setOutputIterationPeriod(10);
		
		//optim.eval();
		optim.dumpParams();
		optim.dumpRayBundles();
		
		optim.optimise(maxInterations);
		
		optim.dumpParams();
		optim.dumpRayBundles();
	}
}
