package fusionOptics.examples.tracer;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.interfaces.Reflector;
import fusionOptics.lenses.Nikon50mmF11;
import fusionOptics.lenses.SchneiderXenon25mmF095;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.DoubleGaussLens;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.pointSpread.DualGaussianPSF;
import fusionOptics.pointSpread.GaussianPSF;
import fusionOptics.pointSpread.MiniImagePSF;
import fusionOptics.pointSpread.PSFGrid;
import fusionOptics.pointSpread.PointSpreadBuilder;
import fusionOptics.pointSpread.PointSpreadFunction;
import fusionOptics.pointSpread.PointsPSF;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Paraboloid;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;
import otherSupport.ColorMaps;
import binaryMatrixFile.BinaryMatrixFile;


import oneLiners.OneLiners;


/** Imaging by two parabolic mirrors
 * 
 * For VRML output only. 
 * 
 * @author oliford
 *
 */
public class ParabolicImagingExample {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/parabolaImaging";
	
	final static double maxTheta = 30 * Math.PI / 180;
	final static int nRaysPerSource = 10000;
	final static double rt2 = Math.sqrt(2);
	
	final static double startY[] = OneLiners.linSpace(-0.5, 0.5, 4);
	final static double startZ[] = OneLiners.linSpace(-0.5, 0.5, 3);
	
	final static double wavelen = 593e-9;
	
	public static void main(String[] args) {
		
		Square backPlane = new Square("backPlane", new double[]{ -0.3, 0, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 10.100, 10.300, Absorber.ideal());
		
		
		double objectToMirrors = 4.0;
		double focusDist = 0.5;
		double sepDist = 0.1;	
		double imageDist = 0.3;
		double aperatureBeforeMirror = 0.050;
		
		/* M1 focus at source/dest */
		double pos[] = { objectToMirrors, 0, 0};
		double focus[] = { objectToMirrors - focusDist, 0.0, 0 };
		double normal[] = {0, 1, 0};
		//*/
		
		/* M2 focus at source/dest */
		double pos2[] = { objectToMirrors, sepDist, 0};
		double focus2[] = { objectToMirrors - focusDist, sepDist, 0 };
		double normal2[] = {0, -1, 0};
		//*/
		
		/* M1 focus behind other mirror */
		/*double pos[] = { 3.0, 0, 0 };
		double focus[] = { 3.0, 1.0, 0 };
		double normal[] = {-1, 0, 0};
		//*/
		
		/* M2 focus behind other mirror */
		/*double pos2[] = { 3.0, sepDist, 0};
		double focus2[] = { 3.0, sepDist - 1.0, 0 };
		double normal2[] = {-1, 0, 0};
		//*/
		
		double apertureDiameter = 0.030;
		double[] aperaturePos = { objectToMirrors - aperatureBeforeMirror, 0, 0};
		Iris iris1 = new Iris("iris1", aperaturePos, new double[] {1,0,0}, apertureDiameter/2*1.5, apertureDiameter/2, Absorber.ideal());
		Disc target1 = new Disc("target1", aperaturePos, new double[] {-1,0,0}, apertureDiameter/2, NullInterface.ideal());
		
		Paraboloid paraboloid1 = new Paraboloid("Paraboloid1", pos, focus, normal, 0.030, null, null, Reflector.ideal());
		
		Paraboloid paraboloid2 = new Paraboloid("Paraboloid2", pos2, focus2, normal2, 0.030, null, null, Reflector.ideal());
				
		Square imgPlane = new Square("imgPlane", new double[]{ objectToMirrors - imageDist, sepDist, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.100, 0.100, Absorber.ideal());
		
		Optic all = new Optic("all", new Element[]{ backPlane, target1, iris1, paraboloid1, paraboloid2, imgPlane });
		
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/imgTest.vrml", 0.05);
		if(vrmlOut != null) {
			vrmlOut.setDrawPolarisationFrames(false);
			vrmlOut.setSkipRays(47);
		}
				
		//double col[][] = ColorMaps.jet(nRaysPerSource);
		double col[][] = ColorMaps.alternating2D2x2(startY.length, startZ.length);
		
		
		for(int iY=0; iY < startY.length; iY++) {
			for(int iZ=0; iZ < startZ.length; iZ++) {
				double startPos[] = new double[] { 0, startY[iY], startZ[iZ] };
	
				for(int i=0; i < nRaysPerSource; i++) {
					
					RaySegment ray = new RaySegment();
					ray.startPos = startPos;
					ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, target1);
					
					ray.length = Double.POSITIVE_INFINITY;
					ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, new double[]{0,0,1})), ray.dir);
					
					
					//ray.E0 = PointSpreadFunction.getInputStatesForMuellerCalc();
					ray.E0 = new double[][]{{1,0,0,0}}; 
					ray.wavelength = wavelen;
					
					Tracer.trace(all, ray, 30, 0.01, true);
					
					if(vrmlOut != null)
						vrmlOut.drawRay(ray, col[iZ*startY.length+iY]);

					ray.processIntersections(imgPlane);
					
					Pol.recoverAll();
				}
				
				System.out.print(".");
			}
		}
		
		//Util.throwAwayIrises(all, 100);
		
		if(vrmlOut != null) {
			vrmlOut.drawOptic(all);
			vrmlOut.destroy();
		}
		
		VRMLDrawer.dumpRay(outPath + "/opticsOnly.vrml", all, null);
		
	}
		
}
