package fusionOptics.examples.optimisationMulti;

import oneLiners.OneLiners;
import otherSupport.ColorMaps;
import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.Reflector;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.optimisationMulti.BendableDish;
import fusionOptics.optimisationMulti.OptimiseMulti;
import fusionOptics.optimisationMulti.RayBundle;
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

/** Simple example for optimisation of a multi-lens imaging system 
 * for overall focus, correct image size and good intensity. 
 * 
 * @author oliford
 */
public class OptimisationMultiDoubleGauss {

	/** Output path for drawings and data */
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/optimMultiTest";
	
	public static void main(String[] args) {
		/** Number of start positions (object plane) */
		int nSourcePoints = 2;
		/** Number of rays from each start position */
		int nRays = 500;
		
		// Wavelength / m 
		double wavelen = 500e-9;
		
		// Create Materials for arbitary fixed-index glasses, and create a medium for each 
		Medium BASF51 = new Medium(new IsotropicFixedIndexGlass(1.7237));
		Medium LAKN7 = new Medium(new IsotropicFixedIndexGlass(1.6516));
		Medium SF1 = new Medium(new IsotropicFixedIndexGlass(1.7174));
		Medium BASF2 = new Medium(new IsotropicFixedIndexGlass(1.6645));
		Medium LAF2 = new Medium(new IsotropicFixedIndexGlass(1.744));
		
		
		
		
		
		double s1Dist = 0.5;
		double s1Pos[] =  new double[]{s1Dist, 0.0, 0.0};
		double s1Curv = 0.073878;
		double s1Norm[] = new double[]{1.0, 0.0, 0.0};
		double s1Thick = 0.007861;
		double s1Radius = 0.0334;
		Medium s1Medium1 = BASF51;
		Medium s1Medium2 = null;
		
		double s2Dist = s1Dist + s1Thick;
		double s2Pos[] =  new double[]{s2Dist, 0.0, 0.0};
		double s2Curv = 0.404472;
		double s2Norm[] = new double[]{1.0, 0.0, 0.0};
		double s2Thick = 0.001284;
		double s2Radius = 0.0334;
		Medium s2Medium1 = null;
		Medium s2Medium2 = BASF51;
		
		double s3Dist = s2Dist + s2Thick;
		double s3Pos[] =  new double[]{s3Dist, 0.0, 0.0};
		double s3Curv = 0.042314;
		double s3Norm[] = new double[]{1.0, 0.0, 0.0};
		double s3Thick = 0.013931;
		double s3Radius = 0.027;
		Medium s3Medium1 = LAKN7;
		Medium s3Medium2 = null;
		
		double s4Dist = s3Dist + s3Thick;
		double s4Pos[] =  new double[]{s4Dist, 0.0, 0.0};
		double s4Curv = 0.271759;
		double s4Norm[] = new double[]{-1.0, 0.0, 0.0};
		double s4Thick = 0.013135;
		double s4Radius = 0.025;
		Medium s4Medium1 = LAKN7;
		Medium s4Medium2 = SF1;
		
		double s5Dist = s4Dist + s4Thick;
		double s5Pos[] =  new double[]{s5Dist, 0.0, 0.0};
		double s5Curv = 0.024961;
		double s5Norm[] = new double[]{1.0, 0.0, 0.0};
		double s5Thick = 0.023285 + 0.004976;
		double s5Radius = 0.0131;
		Medium s5Medium1 = null;
		Medium s5Medium2 = SF1;
		
		
		double s6Dist = s5Dist + s5Thick;
		double s6Pos[] =  new double[]{s6Dist, 0.0, 0.0};
		double s6Curv = 0.031056;
		double s6Norm[] = new double[]{-1.0, 0.0, 0.0};
		double s6Thick = 0.001999;
		double s6Radius = 0.0163;
		Medium s6Medium1 = null;
		Medium s6Medium2 = BASF2;
		
		double s7Dist = s6Dist + s6Thick;
		double s7Pos[] =  new double[]{s7Dist, 0.0, 0.0};
		double s7Curv = 0.077345;
		double s7Norm[] = new double[]{1.0, 0.0, 0.0};
		double s7Thick = 0.010946;
		double s7Radius = 0.0213;
		Medium s7Medium1 = LAF2;
		Medium s7Medium2 = BASF2;
		
		double s8Dist = s7Dist + s7Thick;
		double s8Pos[] =  new double[]{s8Dist, 0.0, 0.0};
		double s8Curv = 0.040486;
		double s8Norm[] = new double[]{-1.0, 0.0, 0.0};
		double s8Thick = 0.000199;
		double s8Radius = 0.0219;
		Medium s8Medium1 = LAF2;
		Medium s8Medium2 = null;
		
		double s9Dist = s8Dist + s8Thick;
		double s9Pos[] =  new double[]{s9Dist, 0.0, 0.0};
		double s9Curv = 0.103569;
		double s9Norm[] = new double[]{1.0, 0.0, 0.0};
		double s9Thick = 0.008657;
		double s9Radius = 0.0255;
		Medium s9Medium1 = LAF2;
		Medium s9Medium2 = null;
		
		double s10Dist = s9Dist + s9Thick;
		double s10Pos[] =  new double[]{s10Dist, 0.0, 0.0};
		double s10Curv = 0.192555;
		double s10Norm[] = new double[]{-1.0, 0.0, 0.0};
		double s10Thick = 0.0;
		double s10Radius = 0.0255;
		Medium s10Medium1 = LAF2;
		Medium s10Medium2 = null;
				
		
		
		
		Iris iris1 = new Iris("iris1", s1Pos, new double[]{ -1.0, 0.0, 0.0}, s1Radius + 0.02, s1Radius, Absorber.ideal());
		Iris iris2 = new Iris("iris1", s5Pos, new double[]{ -1.0, 0.0, 0.0}, s5Radius + 0.02, s5Radius, Absorber.ideal());

		
		Dish lens1Front = new Dish("L1-front", s1Pos, s1Norm, s1Curv, s1Radius, s1Medium1, s1Medium2, IsoIsoStdFresnel.ideal());
		Dish lens1Back = new Dish("L1-back", s2Pos, s2Norm, s2Curv, s2Radius, s2Medium1, s2Medium2, IsoIsoStdFresnel.ideal());
		Dish lens2Front = new Dish("L2-front", s3Pos, s3Norm, s3Curv, s3Radius, s3Medium1, s3Medium2, IsoIsoStdFresnel.ideal());
		Dish lens2Middle = new Dish("L2-middle", s4Pos, s4Norm, s4Curv, s4Radius, s4Medium1, s4Medium2, IsoIsoStdFresnel.ideal());
		Dish lens2Back = new Dish("L2-back", s5Pos, s5Norm, s5Curv, s5Radius, s5Medium1, s5Medium2, IsoIsoStdFresnel.ideal());


		Dish lens3Front = new Dish("L3-front", s6Pos, s6Norm, s6Curv, s6Radius, s6Medium1, s6Medium2, IsoIsoStdFresnel.ideal());
		Dish lens3Middle = new Dish("L3-middle", s7Pos, s7Norm, s7Curv, s7Radius, s7Medium1, s7Medium2, IsoIsoStdFresnel.ideal());
		Dish lens3Back = new Dish("L3-back", s8Pos, s8Norm, s8Curv, s8Radius, s8Medium1, s8Medium2, IsoIsoStdFresnel.ideal());
		Dish lens4Front = new Dish("L4-front", s9Pos, s9Norm, s9Curv, s9Radius, s9Medium1, s9Medium2, IsoIsoStdFresnel.ideal());
		Dish lens4Back = new Dish("L4-back", s10Pos, s10Norm, s10Curv, s10Radius, s10Medium1, s10Medium2, IsoIsoStdFresnel.ideal());
		
		
		
//		
//		
//		
//		SimpleDoubleConvexLensNoSym lens1 = SimpleDoubleConvexLensNoSym.fromRadiusOfCurvAndCentreThickness("lens1", new double[]{ lens1Dist, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.0334, 0.073878, -0.404472, lens1Thick, BASF51, IsoIsoStdFresnel.ideal());
//		Iris iris1 = new Iris("iris1", new double[]{ lens1Dist, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.05, 0.0334, Absorber.ideal());
//		double lens2Dist = lens1Dist + lens1Thick + 0.001294;
//		double lens2Thick= 0.013931;
//		SimpleDoubleConvexLensNoSym lens2 = SimpleDoubleConvexLensNoSym.fromRadiusOfCurvAndCentreThickness("lens2", new double[]{ lens2Dist, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.027, 0.042314, 0.271759, lens2Thick, LAKN7, IsoIsoStdFresnel.ideal());
//		Iris iris2 = new Iris("iris2", new double[]{ lens2Dist, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.06, 0.027, Absorber.ideal());
//		double lens3Dist = lens2Dist + lens2Thick;
//		double lens3Thick= 0.01314;
//		SimpleDoubleConvexLensNoSym lens3 = SimpleDoubleConvexLensNoSym.fromRadiusOfCurvAndCentreThickness("lens3", new double[]{ lens3Dist, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.017, -0.271759, 0.0249610, lens3Thick, SF1, IsoIsoStdFresnel.ideal());
//		Iris iris3 = new Iris("iris3", new double[]{ lens3Dist, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.06, 0.014, Absorber.ideal());
//		double lens4Dist = lens3Dist + lens3Thick + 0.023285 + 0.004976;
//		double lens4Thick= 0.001999;
//		SimpleDoubleConvexLensNoSym lens4 = SimpleDoubleConvexLensNoSym.fromRadiusOfCurvAndCentreThickness("lens4", new double[]{ lens4Dist, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.0213, -0.031056, 0.077345, lens4Thick, BASF2, IsoIsoStdFresnel.ideal());
//		double lens5Dist = lens4Dist + lens4Thick + 0.000001;
//		double lens5Thick= 0.010946;
//		SimpleDoubleConvexLensNoSym lens5 = SimpleDoubleConvexLensNoSym.fromRadiusOfCurvAndCentreThickness("lens5", new double[]{ lens5Dist, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.0213, 0.077345, -0.040486, lens5Thick, LAF2, IsoIsoStdFresnel.ideal());
//		double lens6Dist = lens5Dist + lens5Thick + 0.000199;
//		double lens6Thick= 0.008657;
//		SimpleDoubleConvexLensNoSym lens6 = SimpleDoubleConvexLensNoSym.fromRadiusOfCurvAndCentreThickness("lens6", new double[]{ lens6Dist, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.017, 0.103569, -0.192555, lens6Thick, LAF2, IsoIsoStdFresnel.ideal());
//		
		
		
//		
//		Dish l3Front = lens3.getFrontSurface();
//		Medium l3FrontFrontMedium = l3Front.getFrontMedium();
//		//Double n1 = l3FrontFrontMedium.RefrgetactiveIndex(0, 650.0);
//		//System.out.print("Refractive index 1: " + n1 + "\n");
//		
//		Medium l3FrontBackMedium = l3Front.getBackMedium();
//		Double n2 = l3FrontBackMedium.RefrgetactiveIndex(0, 650.0);
//		System.out.print("Refractive index 2: " + n2 + "\n");
//		
//		Dish l3Back = lens3.getBackSurface();
//		Medium l3BackFrontMedium = l3Back.getFrontMedium();
//		//Double n3 = l3BackFrontMedium.RefrgetactiveIndex(0, 650.0);
//		//System.out.print("Refractive index 3: " + n3 + "\n");
//		
//		Medium l3BackBackMedium = l3Back.getBackMedium();
//		Double n4 = l3BackBackMedium.RefrgetactiveIndex(0, 650.0);
//		System.out.print("Refractive index 4: " + n4 + "\n");

		
		Square imagePlane = new Square("imagePlane", new double[]{ s10Dist + 0.042, 0.0, 0.0}, new double[]{ 1.0, 0.0, 0.0}, new double[]{ 0.0, 0.0, 1.0 }, 0.50, 0.50, Absorber.ideal());
		
		// the 'all' Optic just contains all the other elements (optics and surfaces) 
		Optic all = new Optic("all", new Element[]{ iris1, iris2, lens1Front,  lens1Back, lens2Front, lens2Middle, lens2Back,lens3Front, lens3Middle, lens3Back,lens4Front, lens4Back, imagePlane });
		//Optic all = new Optic("all", new Element[]{ lens4, lens5, iris3, imagePlane });
		
		double sourceCentre[] = new double[]{ 0, 0, 0 };
		//double sourceSeparation[] = new double[]{ 0, 0.1, 0 };
		double rayDir[] = new double[]{1, 0, 0};
		double maxAngle = 10.0/180.0*3.14156;
		double rayLength = 1.0;
				
		//RayBundle rayBundles[] = RayBundle.initRaysImagingAlongLine(lens1, imagePlane, sourceCentre, sourceSeparation, nSourcePoints, wavelen, nRays, 0);
		//RayBundle rayBundles[] = RayBundle.initRaysParallelMultiAngle(lens1,imagePlane, rayDir, maxAngle, rayLength, wavelen, nSourcePoints, nRays, 0);
		RayBundle rayBundles[] = RayBundle.initRaysParallelMultiAngle(lens1Front,imagePlane, rayDir, maxAngle, rayLength, wavelen, nSourcePoints, nRays, 0);

				
		//double targetMagnification = 5.0 / 3.0; // No magnigication for rays from infinity, set targetdist to 40mm
		for(RayBundle rayBundle : rayBundles){
			rayBundle.setIntensityWeight(0.0);
			rayBundle.setTargetWeight(6.0);
			rayBundle.setSharpnessWeight(1.0);
			//double objPos = Util.dot(rayBundle.getSourcePos(), Util.reNorm(sourceSeparation));
			double targDist = 0.04;//  objPos * targetMagnification;
			//System.out.println(objPos + "\t" + targDist);
						
			rayBundle.setTargetPos(new double[]{ targDist, 0.0 });		
		}		
		
		OptimiseMulti optim = new OptimiseMulti();
		optim.addRayBundles(rayBundles);
		optim.setTracingElements(all);
		//optim.addParameter(new MoveableElement(imagePlane, new double[]{ 1.0, 0, 0 }, -2.5, 2.5));		

		optim.addParameter(new BendableDish(lens4Back, 0.8 * lens4Back.getRadiusOfCurv(), 1.2 * lens4Back.getRadiusOfCurv()));
		//Dish lens1Surf2 = (Dish)lens1.getSurfaces().get(1);
		//optim.addParameter(new BendableDish(lens1Surf2, 0.8 * lens1Surf2.getRadiusOfCurv(), 1.2 * lens1Surf2.getRadiusOfCurv()));
		
		//Dish lens2Surf1 = (Dish)lens2.getSurfaces().get(0);
		//optim.addParameter(new BendableDish(lens2Surf1, 0.8 * lens2Surf1.getRadiusOfCurv(), 1.2 * lens2Surf1.getRadiusOfCurv()));
		//Dish lens2Surf2 = (Dish)lens2.getSurfaces().get(1);
		//optim.addParameter(new BendableDish(lens2Surf2, 0.8 * lens2Surf2.getRadiusOfCurv(), 1.2 * lens2Surf2.getRadiusOfCurv()));
		
		optim.setOutputPrefix(outPath + "/optim");
		//optim.addRegularDrawing(outPath + "/optim", SVGRayDrawing.class, 10);
		optim.addRegularDrawing(outPath + "/optim", VRMLDrawer.class, 10);
		optim.setOutputIterationPeriod(10);
		
		//optim.eval();
		optim.dumpParams();
		optim.dumpRayBundles();
		
		optim.optimise(100);
		
		optim.dumpParams();
		optim.dumpRayBundles();
		
	}
}
