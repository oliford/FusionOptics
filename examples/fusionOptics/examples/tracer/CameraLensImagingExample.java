package fusionOptics.examples.tracer;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoStdFresnel;
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
import fusionOptics.surfaces.Iris;
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


/** Simple imaging by a single lens.
 * 
 * For VRML output only. 
 * See GeneratePSFInterpolation for the full fast imaging exercise. 
 * 
 * @author oliford
 *
 */
public class CameraLensImagingExample {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/cameraLensImaging";
	
	final static double maxTheta = 30 * Math.PI / 180;
	final static int nRaysPerSource = 1000;
	final static double rt2 = Math.sqrt(2);
	
	final static double gridDef[][] = {
			{ -0.300, 0.300, 5 }, 
			{ -0.200, 0.200, 5 }, 
			{ -0.200, 0.200, 5 }
		};
	
	final static double imageX[] = OneLiners.linSpace(-0.3, 0.3, 500);
	final static double imageY[] = OneLiners.linSpace(-0.3, 0.3, 500);
	
	final static double wavelen = 593e-9;
	
	public static void main(String[] args) {
		
		//SchneiderXenon25mmF095 lens = new SchneiderXenon25mmF095(); double f = 0.025;	//lens focal length
		Nikon50mmF11 lens = new Nikon50mmF11(); double f = 0.050;	//lens focal length		
		
		double u = 1.000;	//object distance
		double v = 1.0 / (1.0/f - 1.0/u);	//ideal image distance		
		Square backPlane = new Square("backPlane", new double[]{ -0.3, 0, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.100, 0.300, Absorber.ideal());
		
		lens.shift(new double[]{ u, 0, 0 });
		
		Square imgPlane = new Square("imgPlane", new double[]{ u+v, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.040, 0.040, Absorber.ideal());
		
		PointSpreadBuilder psfBuild = new PointSpreadBuilder(imgPlane, outPath + "/psfData.bin");
		PSFGrid psfGrid = new PSFGrid("cameraLensImgTestS25mm", gridDef, DualGaussianPSF.class);
		
		Optic all = new Optic("all", new Element[]{ backPlane, lens, /*lensIris,*/ imgPlane });
		
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/imgTest.vrml", 0.0001);
		if(vrmlOut != null) {
			vrmlOut.setDrawPolarisationFrames(false);
			vrmlOut.setSkipRays(47);
		}
				
		//double col[][] = ColorMaps.jet(nRaysPerSource);
		double col[][] = ColorMaps.alternating2D2x2(psfGrid.getNY(), psfGrid.getNX());
		
		
		for(int iZ=0; iZ < psfGrid.getNZ(); iZ++) {
			for(int iY=0; iY < psfGrid.getNY(); iY++) {
				for(int iX=0; iX < psfGrid.getNX(); iX++) {
					double startPos[] = psfGrid.gridPos(iX, iY, iZ);
					
					psfBuild.startNewPSF(startPos,
							//new GaussianPSF()
							new DualGaussianPSF(20)
							//new PointsPSF()
							//new MiniImagePSF(20, 20)
							);
					
					for(int i=0; i < nRaysPerSource; i++) {
						
						RaySegment ray = new RaySegment();
						ray.startPos = startPos;
						ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, lens.getSurfaces().get(0));
						
						ray.length = Double.POSITIVE_INFINITY;
						ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, new double[]{0,0,1})), ray.dir);
						
						
						ray.E0 = PointSpreadFunction.getInputStatesForMuellerCalc();
						ray.wavelength = wavelen;
						
						Tracer.trace(all, ray, 30, 0.01, true);
						
						if(vrmlOut != null)
							vrmlOut.drawRay(ray, col[iY*psfGrid.getNX()+iX]);

						ray.processIntersections(imgPlane, psfBuild);
						psfBuild.nextCoherentSet();
						
						Pol.recoverAll();
					}
					
					PointSpreadFunction psf = psfBuild.psfDone(true);
					psfGrid.put(iX, iY, iZ, psf);
				
					System.out.print(".");
				}
			}
			System.out.println("\n" + iZ + ", ");
		}
		
		BinaryMatrixFile.mustWrite(outPath + "/psfAllData.bin", psfGrid.getAllData(), false);
		
		Util.throwAwayIrises(all, 100);
		
		if(vrmlOut != null) {
			vrmlOut.drawOptic(all);
			vrmlOut.destroy();
		}
		
		VRMLDrawer.dumpRay(outPath + "/opticsOnly.vrml", all, null);
		
	}
		
}
