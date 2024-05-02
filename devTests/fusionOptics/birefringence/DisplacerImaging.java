package fusionOptics.birefringence;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.OpticApprox;
import fusionOptics.Util;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoUniaxialInterface;
import fusionOptics.lenses.Custom50mmF1;
import fusionOptics.materials.AlphaBariumBorate;
import fusionOptics.materials.UniaxialFixedIndexGlass;
import fusionOptics.optics.Box;
import fusionOptics.optimisationMulti.MoveableElement;
import fusionOptics.optimisationMulti.OptimiseMulti;
import fusionOptics.optimisationMulti.RayBundle;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import net.jafama.FastMath;
import uk.co.oliford.jolu.AsciiMatrixFile;
import uk.co.oliford.jolu.OneLiners;

/** Look at the phase and contrast of a single perfect displacer plate
 *  via a real imaging system (aspheric or camera lens)
 * 
 * 
 * @author oliford */
public class DisplacerImaging {
	public static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/displacerImaging/Custom50mm";//Nikon50mm/opticAxis_up";
	public static VRMLDrawer vrmlOut;
	public static SVGRayDrawing svgOut;
	
	public static double ni = 1.0;
	/*public static double no = 1.54264;
	public static double ne = 1.5517;
	public static double L = 1.973e-3;
	public static double wavelen = 632.8e-9;
	public static Material plateMat = new UniaxialFixedIndexGlass(no, ne);	
	*/
	public static double wavelen = 653.5e-9;
	public static Material plateMat = new AlphaBariumBorate();	
	public static double no = plateMat.getRefractiveIndex(0, wavelen, 300);
	public static double ne = plateMat.getRefractiveIndex(1, wavelen, 300);
	public static double L = 3.8e-3;
	public static double theta = 45 * Math.PI/180; //angle of optic axis to surface
	public static double maxAlpha = 9*Math.PI/180;
	public static double maxXY = 0.030;
	
	public static int na = 7, np = 7;
	
	
	public static double focusScan[] = OneLiners.linSpace(-0.001, 0.001, 50);
	
	
	public static double opticAxisAngleToNormal = Math.PI/2 - theta;
	public static double opticAxisRotation = 0 * Math.PI / 180;
	
	public static double opticAxis[] = new double[]{ 
			FastMath.cos(opticAxisAngleToNormal), 
			FastMath.sin(opticAxisAngleToNormal) * FastMath.cos(opticAxisRotation),
			FastMath.sin(opticAxisAngleToNormal) * FastMath.sin(opticAxisRotation)  };
	
		
	public static Medium plateMedium = new Medium(plateMat, new double[][]{ opticAxis }, 300);
	
	public static Box plate = new Box("plate", new double[]{ L/2, 0, 0 }, L, 0.2, 0.2, plateMedium, IsoUniaxialInterface.ideal());
	
	//public static EdmundOptics50mmAspheric lens = new EdmundOptics50mmAspheric(new double[]{ 0.050, 0, 0 }, new double[]{ 1, 0, 0 });
	public static Custom50mmF1 lens = new Custom50mmF1(true) ;
		
	/*public static BK7 lensMat = new BK7();
	public static Medium lens1Medium = new Medium(lensMat);
	public static SimplePlanarAsphericLens lens1 = SimplePlanarAsphericLens.fromFocalLengthAndEdgeThickness("lens1", 
													new double[]{0.045, 0, 0}, 
													new double[]{1,0,0}, 
													0.025, //radius 
													0.105, //focal length
													0.002, //edge thickness 
													0, //conic const
													new double[]{0, 0, 0, 0, 0}, //poly coeffs
													lens1Medium, 
													IsoIsoInterface.ideal(), 
													wavelen);

	public static Medium lens2Medium = new Medium(lensMat);
	public static SimplePlanarAsphericLens lens2 = SimplePlanarAsphericLens.fromFocalLengthAndEdgeThickness("lens2",
													new double[]{0.060, 0, 0}, 
													new double[]{1,0,0}, 
													0.025, //radius 
													0.105, //focal length
													0.002, //edge thickness 
													0, //conic const
													new double[]{0, 0, 0, 0, 0}, //poly coeffs
													lens2Medium, 
													IsoIsoInterface.ideal(), 
													wavelen);
	*/
	
	//public static Nikon50mmF11 lens = new Nikon50mmF11(new double[]{ 0.100, 0, 0 });
	//public static Nikon135mmF28 lens = new Nikon135mmF28(new double[]{ 0.100, 0, 0 }, 50.0/135);
	
	//public static double imgPlanePos = 0.050 + EdmundOptics50mmAspheric.backFocalDistance -0.000310;
	public static double imgPlanePos = 0.050 + Custom50mmF1.backFocalDistance;
			
	public static double backPos[] = new double[]{ imgPlanePos, 0, 0 };
	public static double backNorm[] = new double[]{ 1, 0, 0};
	
	public static Iris iris = new Iris("iris", new double[]{
			lens.getBoundarySphereCentre()[0] - 0.001, 0, 0 }
				, new double[]{1,0,0}, 0.070, lens.getBoundarySphereRadius()*0.480, Absorber.ideal());
	
	public static Square imagePlane = new Square("back", backPos, backNorm, new double[]{ 0, 0, 1 }, 2, 2, Absorber.ideal());
	
	
	public static Optic all = new Optic("all", new Element[]{ plate, iris, lens, imagePlane });
	

	public static void main(String[] args) {
		lens.shift(new double[]{ 0.050, 0, 0 });
		
		//autoFocus();
		scanAll(true);
		/*
		double startPos = imagePlane.getCentre()[0];
		AsciiMatrixFile.mustWrite(outPath + "/focusScan.txt", focusScan);
		
		String outRoot = outPath;
		for(int i=0; i < focusScan.length; i++){
			System.out.println("Scan " + i + ": dx = " + focusScan[i]);
			imagePlane.setCentre(new double[]{ startPos + focusScan[i], 0, 0 });			
			double contrast[][] = scanAll(false);
			AsciiMatrixFile.mustWrite(outPath + "/contrast_"+i+".txt", contrast, false);
		}
		//*/
	}
	
	public static double[][] scanAll(boolean output) {
		double fw = Util.calcWaveplateFullWaveThickness(new UniaxialFixedIndexGlass(no, ne), wavelen);
		System.out.println(fw/4);
			
		//double alpha[] = OneLiners.linSpace(0, 10*Math.PI/180, 100);
		//double delta[] = OneLiners.linSpace(0, 360*Math.PI/180, 360);
		
		double phsCalc[][] = new double[na*na][np*np];
		double phsTarg[][] = new double[na][na];
		double contrast[][] = new double[na][na];
		
		int nWave = 500;
		
		boolean fail = false;
		
		for(int iAY=0; iAY < na; iAY++){
			double ay = (na==1) ? 0 : ((double)iAY / (0.5*(na-1))) - 1;
			for(int iAX=0; iAX < na; iAX++){
				double ax = (na==1) ? 0 : ((double)iAX / (0.5*(na-1))) - 1;

				double wave[] = new double[nWave];
				int waveSum = 0;
				
				for(int iPY=0; iPY < np; iPY++){
					double py = (np==1) ? 0 : ((double)iPY / (0.5*(np-1))) - 1;
					for(int iPX=0; iPX < np; iPX++){
						double px = (np==1) ? 0 : ((double)iPX / (0.5*(np-1))) - 1;
																	
						double alpha = FastMath.sqrt(ax*ax + ay*ay) * maxAlpha;
						double delta = FastMath.atan2(ay,ax);
						//if(alpha <= maxAlpha){						
							phsCalc[iAY*na+iAX][iPY*np+iPX] = phaseDiffTracer(alpha, delta, px*maxXY, py*maxXY);
							phsTarg[iAY][iAX] = 2*Math.PI*OpticApprox.waveplateOPD(ni, no, ne, theta, delta, alpha, L)/wavelen;
							
							
						//}else{
						//	phsCalc[iAY*na+iAX][iPY*np+iPX] = Double.NaN;
						//	phsTarg[iAY][iAX] = Double.NaN;
						//}
						
						if(!Double.isNaN(phsCalc[iAY*na+iAX][iPY*np+iPX])){
							for(int i=0; i < nWave; i++){
								wave[i] += FastMath.sin(i*2*Math.PI/nWave + phsCalc[iAY*na+iAX][iPY*np+iPX]);
							}
							waveSum++;
						}
						
					}
				}
				
				double max=Double.NEGATIVE_INFINITY, min=Double.POSITIVE_INFINITY;
				for(int i=0; i < nWave; i++){
					wave[i] /= waveSum;
					if(wave[i] < min) min = wave[i];
					if(wave[i] > max) max = wave[i];
				}
				contrast[iAY][iAX] = (max - min) / 2;
				
				System.out.print(".");
			}
			System.out.println(iAY + " / " + na);
				
		}
		
		if(vrmlOut != null)
			vrmlOut.destroy();
		if(svgOut != null)
			svgOut.destroy();
		
		if(output){
			AsciiMatrixFile.mustWrite(outPath + "/VeirasFig9aCalc.txt", phsCalc, false);		
			AsciiMatrixFile.mustWrite(outPath + "/VeirasFig9aTarg.txt", phsTarg, false);
			AsciiMatrixFile.mustWrite(outPath + "/contrast.txt", contrast, false);
		}
	
		return contrast;
	}
	
	private static double phaseDiffTracer(double incidenceAngle, double opticAxisAngToIncidencePlane, double lensX, double lensY){
		
		double incident[] = Util.reNorm(new double[]{
				FastMath.cos(incidenceAngle), 
				FastMath.sin(incidenceAngle) * FastMath.cos(opticAxisAngToIncidencePlane), 
				FastMath.sin(incidenceAngle) * FastMath.sin(opticAxisAngToIncidencePlane) });
		
		//back.setNormal(incident); //align back plane perp to ray as Veiras appears to give difference only in dist // to ray
		
		//System.out.println(incidenceAngle +" " + opticAxisAngToIncidencePlane + " " + lensX + " " + lensY);
		RaySegment ray = new RaySegment();
		
		ray.dir = incident;		
		ray.startPos = new double[]{ -incident[0] , -incident[1] + lensX, -incident[2] + lensY };		
		ray.length = Double.POSITIVE_INFINITY;
		ray.up = Util.reNorm(Util.cross(Util.cross(incident, new double[]{ 0, 0, 1 }), incident));
		ray.E0 = new double[][]{ { 1,0,1,0 } }; // 45° to the right
		ray.wavelength = wavelen;
		
		Tracer.trace(all, ray, 1000, 0.01, true);

		if(vrmlOut == null){
			vrmlOut = new VRMLDrawer(outPath + "/phaseMaps.vrml", 0.002);
			vrmlOut.setDrawPolarisationFrames(true);
			
			vrmlOut.drawOptic(all);			
		}
		if(svgOut == null){
			svgOut = new SVGRayDrawing(outPath + "/phaseMaps.svg", new double[]{-1,-1,-1,1,1,1}, true);
			svgOut.drawElement(all);
		}

		
		
		//check the ray direction
		RaySegment oRay = ray.endHit.transmittedOrdinary;
		RaySegment eRay = ray.endHit.transmittedExtraordinary;
		
		RaySegment oEnd = oRay, eEnd = eRay;

		//find the ends, and strip all reflections
		while(oEnd != null && oEnd.endHit != null && oEnd.endHit.transmittedOrdinary != null){
			oEnd.endHit.reflectedOrdinary = null;
			oEnd.endHit.reflectedExtraordinary = null;
			oEnd = oEnd.endHit.transmittedOrdinary;
			
		}
		while(eEnd != null && eEnd.endHit != null && eEnd.endHit.transmittedOrdinary != null){
			eEnd.endHit.reflectedOrdinary = null;
			eEnd.endHit.reflectedExtraordinary = null;
			eEnd = eEnd.endHit.transmittedOrdinary;			
		}
		
		if(oEnd == null || oEnd.endHit == null || oEnd.endHit.surface != imagePlane ||
			eEnd == null || eEnd.endHit == null || eEnd.endHit.surface != imagePlane){
			
			//
			//VRMLDrawer.dumpRay(outPath + "/displMaps.vrml", all, ray, 0.005);
			return 0;
		}
		
		if(na*np <= 100){
			vrmlOut.drawRay(ray);
			svgOut.drawRay(ray);
		}
		
		double phaseO = FastMath.atan2(oEnd.E1[0][1], oEnd.E1[0][0]);
		double phaseE = FastMath.atan2(eEnd.E1[0][3], eEnd.E1[0][2]);
		
		long fullWaveDiff = oRay.nWaves + oEnd.nWaves - eRay.nWaves - eEnd.nWaves;
		
		double phaseDiffTracer = fullWaveDiff*2*Math.PI + phaseO - phaseE;
		
		Pol.recoverAll();
		return phaseDiffTracer;
	}
	

	private static void autoFocus() {
		int nAngles = 50;
		int nRaysPerAngle = 100;
		
		Optic sys = new Optic("imagingSys", new Element[]{ iris, lens,  imagePlane }); //without plate
		
		RayBundle rayBundles[] = RayBundle.initRaysParallelMultiAngle(lens.getSurfaces().get(0), imagePlane, new double[]{1,0,0 }, 
										maxAlpha, 0.100, wavelen, nAngles, nRaysPerAngle, 0);
		
		//double targetMagnification = 5.0 / 3.0;
		for(RayBundle rayBundle : rayBundles){
			rayBundle.setIntensityWeight(0.0);
			rayBundle.setTargetWeight(0.0);
			rayBundle.setSharpnessWeight(1.0);
			
			//double objPos = Util.dot(rayBundle.getSourcePos(), Util.reNorm(sourceSeparation));
			//double targDist = objPos * targetMagnification;
			//System.out.println(objPos + "\t" + targDist);						
			//rayBundle.setTargetPos(new double[]{ targDist, 0.0 });
		}		
		
		OptimiseMulti optim = new OptimiseMulti();
		optim.addRayBundles(rayBundles);
		optim.setTracingElements(sys);
		
		optim.addParameter(new MoveableElement(imagePlane, new double[]{1,0,0}, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0.001));
		//optim.addParameters(ShapeableAshperic.allAsphericParams(lens.getAsphericSurface()));
		/*optim.addParameter(new BendableDish(lens1.getAsphericSurface(), -10, 10, 100));
		optim.addParameter(new ShapeableAshperic(lens1.getAsphericSurface(), -1, -100, 100, 1e0));
		optim.addParameter(new ShapeableAshperic(lens1.getAsphericSurface(), 1, -100, 100, 1e3));
		optim.addParameter(new ShapeableAshperic(lens1.getAsphericSurface(), 2, -100, 100, 1e5));
		optim.addParameter(new ShapeableAshperic(lens1.getAsphericSurface(), 3, -100, 100, 1e8));
		optim.addParameter(new ShapeableAshperic(lens1.getAsphericSurface(), 4, -100, 100, 1e10));
		

		optim.addParameter(new BendableDish(lens2.getAsphericSurface(), -10, 10, 100));
		optim.addParameter(new ShapeableAshperic(lens2.getAsphericSurface(), -1, -100, 100, 1e0));
		optim.addParameter(new ShapeableAshperic(lens2.getAsphericSurface(), 1, -100, 100, 1e3));
		optim.addParameter(new ShapeableAshperic(lens2.getAsphericSurface(), 2, -100, 100, 1e5));
		optim.addParameter(new ShapeableAshperic(lens2.getAsphericSurface(), 3, -100, 100, 1e8));
		optim.addParameter(new ShapeableAshperic(lens2.getAsphericSurface(), 4, -100, 100, 1e10));
		*/

		optim.setOutputPrefix(outPath + "/autoFocus");
		//optim.addRegularDrawing(outPath + "/autoFocus", SVGRayDrawing.class, 10);
		optim.addRegularDrawing(outPath + "/autoFocus", VRMLDrawer.class, 10);
		optim.setOutputIterationPeriod(10);
		
		//optim.eval();
		optim.dumpParams();
		optim.dumpRayBundles();
		
		optim.optimise(20000);
		
		optim.dumpParams();
		optim.dumpRayBundles();
	}
}
