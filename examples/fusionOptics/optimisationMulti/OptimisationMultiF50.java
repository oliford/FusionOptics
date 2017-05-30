package fusionOptics.optimisationMulti;

import oneLiners.OneLiners;
import otherSupport.ColorMaps;
import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.Reflector;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Plane;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;
import fusionOptics.lenses.Custom50mmF1;

/** Simple example for optimisation of a multi-lens imaging system 
 * for overall focus, correct image size and good intensity. 
 * 
 * @author oliford
 */
public class OptimisationMultiF50 {

	/** Output path for drawings and data */
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/optimMultiTest";
	
	public static void main(String[] args) {
		/** Number of start positions (object plane) */
		int nSourcePoints = 2;
		/** Number of rays from each start position */
		int nRays = 5000;
		
		// Wavelength / m 
		double wavelen = 653e-9;
		
		// Create Materials for arbitary fixed-index glasses, and create a medium for each 
//		Medium SLAM7 = new Medium(new IsotropicFixedIndexGlass(1.74352613));
//		Medium ZNSE= new Medium(new IsotropicFixedIndexGlass(2.57887416));
//		Medium NSF57 = new Medium(new IsotropicFixedIndexGlass(1.83690172));
//		
//		
//		
//		// Have to reverse all coefficients for negative surfaces?
//
//		
//		double s1Dist = 0.5;
//		double s1Pos[] =  new double[]{s1Dist, 0.0, 0.0};
//		double s1Curv = 3748.74e-3;
//		double s1Norm[] = new double[]{-1.0, 0.0, 0.0};
//		double s1Thick = 14.52579e-3;
//		double s1Radius = 27.0e-3;
//		Medium s1Medium1 = null;
//		Medium s1Medium2 = NSF57;
//		double s1ConicConst = 0.0; // As order 06/02/15
//		double s1PolyCoeffs[] = new double[]{ 0, 3.4662607e-006, -9.5847605e-011, -4.6348623e-014 };
//		s1PolyCoeffs = Aspheric.rescaleCoeffs(s1PolyCoeffs, 1e-3);
//		
//		double s2Dist = s1Dist + s1Thick;
//		double s2Pos[] =  new double[]{s2Dist, 0.0, 0.0};
//		double s2Curv = 192.7492e-3;
//		double s2Norm[] = new double[]{1.0, 0.0, 0.0};
//		double s2Thick = 26.50819e-3;
//		double s2Radius = 31.0e-3;
//		Medium s2Medium1 = null;
//		Medium s2Medium2 = NSF57;
//		double s2ConicConst = 0.0; // As order 06/02/15
//		double s2PolyCoeffs[] = new double[]{ 0, -2.8773777e-006, 7.5466219e-010, -5.3131379e-014 };
//		s2PolyCoeffs = Aspheric.rescaleCoeffs(s2PolyCoeffs, 1e-3);
//
//		double s3Dist = s2Dist + s2Thick;
//		double s3Pos[] =  new double[]{s3Dist, 0.0, 0.0};
//		double s3Curv = 170.8094e-3;
//		double s3Norm[] = new double[]{-1.0, 0.0, 0.0};
//		double s3Thick = 0.010;
//		double s3Radius = 0.040;
//		Medium s3Medium1 = null;
//		Medium s3Medium2 = NSF57;
//		
//		double s4Dist = s3Dist + s3Thick;
//		double s4Pos[] =  new double[]{s4Dist, 0.0, 0.0};
//		double s4Curv = 91.92364e-3;
//		double s4Norm[] = new double[]{-1.0, 0.0, 0.0};
//		double s4Thick = 10.96722e-3;
//		double s4Radius = 0.041;
//		Medium s4Medium1 = NSF57;
//		Medium s4Medium2 = null;
//		
//		double s5Dist = s4Dist + s4Thick;
//		double s5Pos[] =  new double[]{s5Dist, 0.0, 0.0};
//		double s5Curv = 138.1855e-3;
//		double s5Norm[] = new double[]{1.0, 0.0, 0.0};
//		double s5Thick = 14.92255e-3;
//		double s5Radius = 46.0e-3;
//		Medium s5Medium1 = ZNSE;
//		Medium s5Medium2 = null;
//		
//		double s6Dist = s5Dist + s5Thick;
//		double s6Pos[] =  new double[]{s6Dist, 0.0, 0.0};
//		double s6Curv = 1046.628e-3;
//		double s6Norm[] = new double[]{-1.0, 0.0, 0.0};
//		double s6Thick = 2.0e-3;
//		double s6Radius = 46.0e-3;
//		Medium s6Medium1 = ZNSE;
//		Medium s6Medium2 = null;
//		
//		double s7Dist = s6Dist + s6Thick;
//		double s7Pos[] =  new double[]{s7Dist, 0.0, 0.0};
//		double s7Curv = 41.44526e-3;
//		double s7Norm[] = new double[]{1.0, 0.0, 0.0};
//		double s7Thick = 14.9948e-3;
//		double s7Radius = 35.0e-3;
//		Medium s7Medium1 = ZNSE;
//		Medium s7Medium2 = null;
//		double s7ConicConst = 0.0; // As order 06/02/15
//		double s7PolyCoeffs[] = new double[]{ 0, -1.4887119e-008, -1.5803008e-011, -1.4029144e-015 };
//		s7PolyCoeffs = Aspheric.rescaleCoeffs(s7PolyCoeffs, 1e-3);
//
//		double s8Dist = s7Dist + s7Thick;
//		double s8Pos[] =  new double[]{s8Dist, 0.0, 0.0};
//		double s8Curv = 31.02156e-3;
//		double s8Norm[] = new double[]{1.0, 0.0, 0.0};
//		double s8Thick =  14.10915e-3;
//		double s8Radius = 26e-3;
//		Medium s8Medium1 = null;
//		Medium s8Medium2 = ZNSE;
//		
//		double s9Dist = s8Dist + s8Thick;
//		double s9Pos[] =  new double[]{s9Dist, 0.0, 0.0};
//		double s9Curv = 146.0394e-3;
//		double s9Norm[] = new double[]{1.0, 0.0, 0.0};
//		double s9Thick = 9.782205e-3;
//		double s9Radius = 24.0e-3;
//		Medium s9Medium1 = SLAM7;
//		Medium s9Medium2 = null;
//		
//		double s10Dist = s9Dist + s9Thick;
//		double s10Pos[] =  new double[]{s10Dist, 0.0, 0.0};
//		double s10Curv = 242.5449e-3;
//		double s10Norm[] = new double[]{1.0, 0.0, 0.0};
//		double s10Thick = 0.0;
//		double s10Radius = 20.0e-3;
//		Medium s10Medium1 = null;
//		Medium s10Medium2 = SLAM7;
//		double s10ConicConst = 0.0; // As order 06/02/15
//		double s10PolyCoeffs[] = new double[]{ 0, 2.1353608e-006, 2.2377491e-010, 2.0783215e-012 };
//		s10PolyCoeffs = Aspheric.rescaleCoeffs(s10PolyCoeffs, 1e-3);
//	
//		
//		
//		
//		
//		Iris iris1 = new Iris("iris1", new double[] {s1Dist - 10.078e-3, 0.0, 0.0}, new double[]{ -1.0, 0.0, 0.0},  30.0e-3, 10.0e-3, Absorber.ideal());
//		Iris iris2 = new Iris("iris2", s1Pos, new double[]{ -1.0, 0.0, 0.0}, s1Radius - 0.005, 0.001, Absorber.ideal());
//		Iris iris3 = new Iris("iris3", new double[] {s8Dist + 0.012, 0.0, 0.0}, new double[]{ -1.0, 0.0, 0.0}, s8Radius + 0.05, 50.30929e-3/2.0, Absorber.ideal());
//		
//		//new Aspheric(name + "-curved", curvedSurfaceCentre, reverseAxis, -radiusOfCurvatureFront, 1.0001 * clearAperture / 2, conicConst, polyCoeffs, null, medium, IsoIsoInterface.ideal());
//		Aspheric lens1Front = new Aspheric("L1-front", s1Pos, s1Norm, s1Curv, s1Radius, s1ConicConst, s1PolyCoeffs, s1Medium1, s1Medium2, IsoIsoStdFresnel.ideal());
//		//Dish lens1Front = new Dish("L1-front", s1Pos, s1Norm, s1Curv, s1Radius, s1Medium1, s1Medium2, IsoIsoStdFresnel.ideal());
//		//Dish lens1Back = new Dish("L1-back", s2Pos, s2Norm, s2Curv, s2Radius, s2Medium1, s2Medium2, IsoIsoStdFresnel.ideal());
//		Aspheric lens1Back = new Aspheric("L1-back", s2Pos, s2Norm, s2Curv, s2Radius, s2ConicConst, s2PolyCoeffs, s2Medium1, s2Medium2, IsoIsoStdFresnel.ideal());
//		Dish lens2Front = new Dish("L2-front", s3Pos, s3Norm, s3Curv, s3Radius, s3Medium1, s3Medium2, IsoIsoStdFresnel.ideal());
//		Dish lens2Back = new Dish("L2-back", s4Pos, s4Norm, s4Curv, s4Radius, s4Medium1, s4Medium2, IsoIsoStdFresnel.ideal());
//		Dish lens3Front = new Dish("L3-front", s5Pos, s5Norm, s5Curv, s5Radius, s5Medium1, s5Medium2, IsoIsoStdFresnel.ideal());
//
//
//		Dish lens3Back = new Dish("L3-back", s6Pos, s6Norm, s6Curv, s6Radius, s6Medium1, s6Medium2, IsoIsoStdFresnel.ideal());
//		//Dish lens4Front = new Dish("L3-middle", s7Pos, s7Norm, s7Curv, s7Radius, s7Medium1, s7Medium2, IsoIsoStdFresnel.ideal());
//		Aspheric lens4Front = new Aspheric("L4-front", s7Pos, s7Norm, s7Curv, s7Radius, s7ConicConst, s7PolyCoeffs, s7Medium1, s7Medium2, IsoIsoStdFresnel.ideal());
//		Dish lens4Back = new Dish("L4-back", s8Pos, s8Norm, s8Curv, s8Radius, s8Medium1, s8Medium2, IsoIsoStdFresnel.ideal());
//		Dish lens5Front = new Dish("L5-front", s9Pos, s9Norm, s9Curv, s9Radius, s9Medium1, s9Medium2, IsoIsoStdFresnel.ideal());
//		//Dish lens5Back = new Dish("L4-back", s10Pos, s10Norm, s10Curv, s10Radius, s10Medium1, s10Medium2, IsoIsoStdFresnel.ideal());
//		Aspheric lens5Back = new Aspheric("L5-back", s10Pos, s10Norm, s10Curv, s10Radius, s10ConicConst, s10PolyCoeffs, s10Medium1, s10Medium2, IsoIsoStdFresnel.ideal());
//
//		
//		
//
//		
//		Square imagePlane = new Square("imagePlane", new double[]{ s10Dist + 17.5304e-3, 0.0, 0.0}, new double[]{ 1.0, 0.0, 0.0}, new double[]{ 0.0, 0.0, 1.0 }, 0.50, 0.50, Absorber.ideal());
//		
		// the 'all' Optic just contains all the other elements (optics and surfaces) 
//		Optic all = new Optic("all", new Element[]{ iris1,  iris3, lens1Front,  lens1Back, lens2Front, lens2Back,lens3Front, lens3Back,lens4Front, lens4Back,lens5Front, lens5Back, imagePlane });
		//Optic all = new Optic("all", new Element[]{ lens4, lens5, iris3, imagePlane });
		Custom50mmF1 lens = new Custom50mmF1() ;
		Square imagePlane = new Square("imagePlane", new double[]{ Custom50mmF1.backFocalDistance, 0.0, 0.0}, new double[]{ 1.0, 0.0, 0.0}, new double[]{ 0.0, 0.0, 1.0 }, 0.50, 0.50, Absorber.ideal());
		Optic all = new Optic("all", new Element[]{ lens, imagePlane });
		
		double sourceCentre[] = new double[]{ 0, 0, 0 };
		//double sourceSeparation[] = new double[]{ 0, 0.1, 0 };
		double rayDir[] = new double[]{1, 0, 0};
		double maxAngle = 10.0/180.0*3.14156;
		double rayLength = 1.0;
				
		//RayBundle rayBundles[] = RayBundle.initRaysImagingAlongLine(lens1, imagePlane, sourceCentre, sourceSeparation, nSourcePoints, wavelen, nRays, 0);
		//RayBundle rayBundles[] = RayBundle.initRaysParallelMultiAngle(lens1,imagePlane, rayDir, maxAngle, rayLength, wavelen, nSourcePoints, nRays, 0);
		
		Dish lens1Surf1 = (Dish)lens.getSurfaces().get(1);
		
		System.out.println("\n\nLens1Surf1:  " + lens1Surf1.getName() + "\n\n");
		
		RayBundle rayBundles[] = RayBundle.initRaysParallelMultiAngle(lens1Surf1, imagePlane, rayDir, maxAngle, rayLength, wavelen, nSourcePoints, nRays, 0);
		
		OptimiseMulti optim = new OptimiseMulti();
		optim.addRayBundles(rayBundles);
		optim.setTracingElements(all);
				
		//double targetMagnification = 5.0 / 3.0; // No magnigication for rays from infinity, set targetdist to 40mm
		for(RayBundle rayBundle : rayBundles){
			rayBundle.setIntensityWeight(0.0);
			rayBundle.setTargetWeight(0.0);
			rayBundle.findTarget();
			rayBundle.setSharpnessWeight(1.0);
			//double objPos = Util.dot(rayBundle.getSourcePos(), Util.reNorm(sourceSeparation));
			//double targDist = 0.04;//  objPos * targetMagnification;
			//System.out.println(objPos + "\t" + targDist);
						
			//rayBundle.setTargetPos(new double[]{ targDist, 0.0 });		
		}		
		

		optim.addParameter(new MoveableElement(imagePlane, new double[]{ 1.0, 0, 0 }, -2.5, 2.5));		

//		optim.addParameters(ShapeableAshperic.allAsphericParams(lens1Front));
//		optim.addParameters(ShapeableAshperic.allAsphericParams(lens1Back));
//		optim.addParameters(ShapeableAshperic.allAsphericParams(lens4Front));
//		optim.addParameters(ShapeableAshperic.allAsphericParams(lens5Back));
//
//		
//		optim.addParameter(new BendableDish(lens2Front, 0.8 * lens2Front.getRadiusOfCurv(), 1.2 * lens2Front.getRadiusOfCurv()));
//		optim.addParameter(new BendableDish(lens2Back, 0.8 * lens2Back.getRadiusOfCurv(), 1.2 * lens2Back.getRadiusOfCurv()));
//		optim.addParameter(new BendableDish(lens3Front, 0.8 * lens3Front.getRadiusOfCurv(), 1.2 * lens3Front.getRadiusOfCurv()));
//		optim.addParameter(new BendableDish(lens3Back, 0.8 * lens3Back.getRadiusOfCurv(), 1.2 * lens3Back.getRadiusOfCurv()));
//		optim.addParameter(new BendableDish(lens4Back, 0.8 * lens4Back.getRadiusOfCurv(), 1.2 * lens4Back.getRadiusOfCurv()));
//		optim.addParameter(new BendableDish(lens5Front, 0.8 * lens5Front.getRadiusOfCurv(), 1.2 * lens5Front.getRadiusOfCurv()));
//
//		
		//optim.addParameter(new BendableDish(lens2Front, 0.8 * lens2Front.getRadiusOfCurv(), 1.2 * lens2Front.getRadiusOfCurv()));
		//optim.addParameter(new BendableDish(lens2Back, 0.8 * lens2Back.getRadiusOfCurv(), 1.2 * lens2Back.getRadiusOfCurv()));
		//optim.addParameter(new BendableDish(lens3Front, 0.8 * lens3Front.getRadiusOfCurv(), 1.2 * lens3Front.getRadiusOfCurv()));
		//optim.addParameter(new BendableDish(lens4Back, 0.8 * lens4Back.getRadiusOfCurv(), 1.2 * lens4Back.getRadiusOfCurv()));
		//Dish lens1Surf2 = (Dish)lens1.getSurfaces().get(1);
		//optim.addParameter(new BendableDish(lens1Surf2, 0.8 * lens1Surf2.getRadiusOfCurv(), 1.2 * lens1Surf2.getRadiusOfCurv()));
		
		//Dish lens2Surf1 = (Dish)lens2.getSurfaces().get(0);
		//optim.addParameter(new BendableDish(lens2Surf1, 0.8 * lens2Surf1.getRadiusOfCurv(), 1.2 * lens2Surf1.getRadiusOfCurv()));
		//Dish lens2Surf2 = (Dish)lens2.getSurfaces().get(1);
		//optim.addParameter(new BendableDish(lens2Surf2, 0.8 * lens2Surf2.getRadiusOfCurv(), 1.2 * lens2Surf2.getRadiusOfCurv()));
		
		optim.setOutputPrefix(outPath + "/optim");
		optim.addRegularDrawing(outPath + "/optim", SVGRayDrawing.class, 10);
		optim.addRegularDrawing(outPath + "/optim", VRMLDrawer.class, 10);
		optim.setOutputIterationPeriod(10);
		
		//optim.eval();
		optim.dumpParams();
		optim.dumpRayBundles();
		
		optim.optimise(100);
		
		optim.dumpParams();
		optim.dumpRayBundles();
		optim.destroy();
	}
}
