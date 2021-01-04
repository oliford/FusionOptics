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
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
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
public class OptimisationMultiTest {

	/** Output path for drawings and data */
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/optimMultiTest";
	
	public static void main(String[] args) {
		/** Number of start positions (object plane) */
		int nSourcePoints = 8;
		/** Number of rays from each start position */
		int nRays = 200;
		
		// Wavelength / m 
		double wavelen = 500e-9;
		
		// Create Materials for arbitary fixed-index glasses, and create a medium for each 
		Medium glass120 = new Medium(new IsotropicFixedIndexGlass(1.2));
		Medium glass130 = new Medium(new IsotropicFixedIndexGlass(1.3));
		Medium glass105 = new Medium(new IsotropicFixedIndexGlass(1.05));
		
		SimpleDoubleConvexLens lens1 = SimpleDoubleConvexLens.fromFocalLengthAndEdgeThickness("lens1", new double[]{ 2.0, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.4, 2.0, 0.005, glass130, IsoIsoStdFresnel.ideal(), wavelen);
		Iris iris1 = new Iris("iris1", new double[]{ 2.0, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.8, 0.4, Absorber.ideal());
		SimpleDoubleConvexLens lens2 = SimpleDoubleConvexLens.fromFocalLengthAndEdgeThickness("lens2", new double[]{ 3.0, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.5, 5.0, 0.005, glass120, IsoIsoStdFresnel.ideal(), wavelen);
		Iris iris2 = new Iris("iris2", new double[]{ 3.0, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.8, 0.5, Absorber.ideal());
		
		Square imagePlane = new Square("imagePlane", new double[]{ 8.5, 0.0, 0.0}, new double[]{ 1.0, 0.0, 0.0}, new double[]{ 0.0, 0.0, 1.0 }, 5.0, 5.0, Absorber.ideal());
		
		// the 'all' Optic just contains all the other elements (optics and surfaces) 
		Optic all = new Optic("all", new Element[]{ lens1, iris1, lens2, iris2, imagePlane });
		
		double sourceCentre[] = new double[]{ 0, 0, 0 };
		double sourceSeparation[] = new double[]{ 0, 0.1, 0 };
				
		RayBundle rayBundles[] = RayBundle.initRaysImagingAlongLine(lens1, imagePlane, sourceCentre, sourceSeparation, nSourcePoints, wavelen, nRays, 0);
		
		double targetMagnification = 5.0 / 3.0;
		for(RayBundle rayBundle : rayBundles){
			rayBundle.setIntensityWeight(0.0);
			rayBundle.setTargetWeight(6.0);
			rayBundle.setSharpnessWeight(1.0);
			double objPos = Util.dot(rayBundle.getSourcePos(), Util.reNorm(sourceSeparation));
			double targDist = objPos * targetMagnification;
			System.out.println(objPos + "\t" + targDist);
						
			rayBundle.setTargetPos(new double[]{ targDist, 0.0 });		
		}		
		
		OptimiseMulti optim = new OptimiseMulti();
		optim.addRayBundles(rayBundles);
		optim.setTracingElements(all);
		optim.addParameter(new MoveableElement(imagePlane, new double[]{ 1.0, 0, 0 }, -2.5, 2.5));		

		Dish lens1Surf1 = (Dish)lens1.getSurfaces().get(0);
		optim.addParameter(new BendableDish(lens1Surf1, 0.8 * lens1Surf1.getRadiusOfCurv(), 1.2 * lens1Surf1.getRadiusOfCurv()));
		Dish lens1Surf2 = (Dish)lens1.getSurfaces().get(1);
		optim.addParameter(new BendableDish(lens1Surf2, 0.8 * lens1Surf2.getRadiusOfCurv(), 1.2 * lens1Surf2.getRadiusOfCurv()));
		
		Dish lens2Surf1 = (Dish)lens2.getSurfaces().get(0);
		optim.addParameter(new BendableDish(lens2Surf1, 0.8 * lens2Surf1.getRadiusOfCurv(), 1.2 * lens2Surf1.getRadiusOfCurv()));
		Dish lens2Surf2 = (Dish)lens2.getSurfaces().get(1);
		optim.addParameter(new BendableDish(lens2Surf2, 0.8 * lens2Surf2.getRadiusOfCurv(), 1.2 * lens2Surf2.getRadiusOfCurv()));
		
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
