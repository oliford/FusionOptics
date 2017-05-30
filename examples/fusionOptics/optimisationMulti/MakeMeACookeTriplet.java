package fusionOptics.optimisationMulti;

import oneLiners.OneLiners;
import otherSupport.ColorMaps;
import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.Reflector;
import fusionOptics.materials.BK7;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.materials.SF11;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.optics.SimplePlanarConvexLens;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Plane;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;

/** Simple example for optimisation of a multi-lens imaging system 
 * for overall focus and correct image size  
 * 
 * @author oliford
 */
public class MakeMeACookeTriplet {

	/** Output path for drawings and data */
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/makeMeACookeTriplet";
	
	public static void main(String[] args) {
		/** Number of start positions (object plane) */
		int nSourcePoints = 5;
		/** Number of rays from each start position */
		int nRays = 500;
		
		
		// Wavelength / m 
		double wavelengths[] = new double[]{ 500e-9, 600e-9 };
		
		// Material 1: for the 2 positive lenses (high Abbe number) 
		Material highVMat = new BK7();
		// Material 1: for the negative lens (low Abbe number)
		Material lowVMat = new SF11();
		
		//Initial positions, curvatures etc
		double objectDist = 1.300;
		
		double targetMagnification = 0.200;
		double imageDist = objectDist * targetMagnification;
				
		double focalLength = 1.0 / (1.0 / objectDist + 1.0 / imageDist); 
		double simpleCurvature = (highVMat.getRefractiveIndex(0, wavelengths[0], 300)-1)*focalLength;
				
		double l1CentreThickness = 0.020;
		double l1curvatureRadius1 = 3*simpleCurvature;
		double l1curvatureRadius2 = 3*simpleCurvature;
		double l1Diameter = 0.080;
		
		double l1l2Dist = 0.020; 
		double l2CentreThickness = 0.010;
		double l2curvatureRadius1 = 10*simpleCurvature;
		double l2curvatureRadius2 = 10*simpleCurvature;
		double l2Diameter = 0.060;
		
		double l2l3Dist = 0.020; 
		double l3CentreThickness = 0.020;
		double l3curvatureRadius1 = 3*simpleCurvature;
		double l3curvatureRadius2 = 3*simpleCurvature;
		double l3Diameter = 0.080;
		
		
		
		double l1S1Pos = objectDist - l2CentreThickness/2 - l1l2Dist - l1CentreThickness;
		double l2S1Pos = objectDist - l2CentreThickness/2;
		double l3S1Pos = objectDist + l2CentreThickness/2 + l2l3Dist;
		
		
		//Create optic objects
		double opticAxis[] = new double[]{ 1, 0, 0 };
		double opticAxisRev[] = new double[]{ -1, 0, 0 };
		
		Medium l1Medium = new Medium(highVMat);
		Dish l1Front = new Dish("l1Front", new double[]{ l1S1Pos, 0,0 }, opticAxis, l1curvatureRadius1, l1Diameter*0.55, l1Medium, null, IsoIsoStdFresnel.ideal());
		Dish l1Back = new Dish("l1Back", new double[]{ l1S1Pos+l1CentreThickness, 0,0 }, opticAxisRev, l1curvatureRadius2, l1Diameter*0.55, l1Medium, null, IsoIsoStdFresnel.ideal());
		Iris l1Iris1 = new Iris("l1Iris1", new double[]{ l1S1Pos + l1CentreThickness/2, 0,0 }, new double[]{ -1.0, 0.0, 0.0}, l1Diameter*1.00, l1Diameter*0.5, Absorber.ideal());
		
		Medium l2Medium = new Medium(lowVMat);
		Dish l2Front = new Dish("l2Front", new double[]{ l2S1Pos, 0,0 }, opticAxisRev, l2curvatureRadius1, l2Diameter*0.55, null, l2Medium, IsoIsoStdFresnel.ideal());
		Dish l2Back = new Dish("l2Back", new double[]{ l2S1Pos+l2CentreThickness, 0,0 }, opticAxis, l2curvatureRadius2, l2Diameter*0.55, null, l2Medium, IsoIsoStdFresnel.ideal());
		Iris l2Iris1 = new Iris("l2Iris1", new double[]{ l2S1Pos + l2CentreThickness/2, 0,0 }, new double[]{ -1.0, 0.0, 0.0}, l2Diameter*1.00, l2Diameter*0.5, Absorber.ideal());
		
		Medium l3Medium = new Medium(highVMat);
		Dish l3Front = new Dish("l3Front", new double[]{ l3S1Pos, 0,0 }, opticAxis, l3curvatureRadius1, l3Diameter*0.55, l3Medium, null, IsoIsoStdFresnel.ideal());
		Dish l3Back = new Dish("l3Back", new double[]{ l3S1Pos+l3CentreThickness, 0,0 }, opticAxisRev, l3curvatureRadius2, l3Diameter*0.55, l3Medium, null, IsoIsoStdFresnel.ideal());
		Iris l3Iris1 = new Iris("l3Iris1", new double[]{ l3S1Pos + l3CentreThickness/2, 0,0 }, new double[]{ -1.0, 0.0, 0.0}, l3Diameter*1.00, l3Diameter*0.5, Absorber.ideal());
		
		Square imagePlane = new Square("imagePlane", new double[]{ objectDist + imageDist, 0.0, 0.0}, new double[]{ 1.0, 0.0, 0.0}, new double[]{ 0.0, 0.0, 1.0 }, 0.50, 0.50, Absorber.ideal());
		
		// the 'all' Optic just contains all the other elements (optics and surfaces) 
		Optic all = new Optic("all", new Element[]{ l1Front, l1Iris1, l1Back,  l2Front, l2Iris1, l2Back,  l3Front, l3Iris1, l3Back, imagePlane });
		
		
		double sourceCentre[] = new double[]{ 0, 0, 0 };
		double sourceSeparation[] = new double[]{ 0, 0.100, 0 };
		double rayDir[] = new double[]{1, 0, 0};
		//double maxAngle = 10.0/180.0*3.14156;
		//double rayLength = 1.0;
				
		
		RayBundle rayBundles[] = new RayBundle[nSourcePoints*wavelengths.length];
		for(int i=0; i < wavelengths.length; i++){
			RayBundle rayBundlesWL[] = RayBundle.initRaysImagingAlongLine(l1Front, imagePlane, sourceCentre, sourceSeparation, nSourcePoints, wavelengths[i], nRays, 0);
			//RayBundle rayBundles[] = RayBundle.initRaysParallelMultiAngle(l1Front,imagePlane, rayDir, maxAngle, rayLength, wavelen, nSourcePoints, nRays, 0);
			//RayBundle rayBundles[] = RayBundle.initRaysParallelMultiAngle(l1Front, imagePlane, rayDir, maxAngle, rayLength, wavelen, nSourcePoints, nRays, 0);
			
			System.arraycopy(rayBundlesWL, 0, rayBundles, i*nSourcePoints, nSourcePoints);
		}
		
		//double targetMagnification = 5.0 / 3.0; // No magnigication for rays from infinity, set targetdist to 40mm
		for(RayBundle rayBundle : rayBundles){
			rayBundle.setIntensityWeight(0.0);
			rayBundle.setTargetWeight(0.2);
			rayBundle.setSharpnessWeight(1.0);
			
			double objPos = Util.dot(rayBundle.getSourcePos(), Util.reNorm(sourceSeparation));
			double targDist = objPos * targetMagnification;
			//System.out.println(objPos + "\t" + targDist);
						
			rayBundle.setTargetPos(new double[]{ targDist, 0.0 });		
		}		
		
		OptimiseMulti optim = new OptimiseMulti();
		optim.addRayBundles(rayBundles);
		optim.setTracingElements(all);
		
		optim.addParameter(new MoveableElement(imagePlane, new double[]{ 1.0, 0, 0 }, -2.5, 2.5));		

		optim.addParameter(new BendableDish(l1Front, 0.8 * l1Front.getRadiusOfCurv(), 1.2 * l1Front.getRadiusOfCurv()));
		optim.addParameter(new BendableDish(l1Back, 0.8 * l1Front.getRadiusOfCurv(), 1.2 * l1Back.getRadiusOfCurv()));
		optim.addParameter(new MoveableElement(l1Front, opticAxis, -0.010, 0.010, 0.0001));
		optim.addParameter(new MoveableElement(l1Back, opticAxis, -0.010, 0.010, 0.0001));
		
		optim.addParameter(new BendableDish(l2Front, 0.8 * l2Front.getRadiusOfCurv(), 1.2 * l2Front.getRadiusOfCurv()));
		optim.addParameter(new BendableDish(l2Back, 0.8 * l2Front.getRadiusOfCurv(), 1.2 * l2Back.getRadiusOfCurv()));
		optim.addParameter(new MoveableElement(l2Front, opticAxis, -0.010, 0.010, 0.0001));
		optim.addParameter(new MoveableElement(l2Back, opticAxis, -0.010, 0.010, 0.0001));
		
		optim.addParameter(new BendableDish(l3Front, 0.8 * l3Front.getRadiusOfCurv(), 1.2 * l3Front.getRadiusOfCurv()));
		optim.addParameter(new BendableDish(l3Back, 0.8 * l3Front.getRadiusOfCurv(), 1.2 * l3Back.getRadiusOfCurv()));
		optim.addParameter(new MoveableElement(l3Front, opticAxis, -0.010, 0.010, 0.0001));
		optim.addParameter(new MoveableElement(l3Back, opticAxis, -0.010, 0.010, 0.0001));
		//*/
		
		optim.setOutputPrefix(outPath + "/optim");
		//optim.addRegularDrawing(outPath + "/optim", SVGRayDrawing.class, 10);
		optim.addRegularDrawing(outPath + "/optim", VRMLDrawer.class, 10);
		optim.setOutputIterationPeriod(10);
		
		//optim.eval();
		optim.dumpParams();
		optim.dumpRayBundles();
		
		optim.optimise(500000);
		
		optim.dumpParams();
		optim.dumpRayBundles();
		
	}
}
