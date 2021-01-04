package fusionOptics.examples.pointSpread;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.DoubleGaussLens;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.optimisation.AutoFocusOld;
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
import otherSupport.ColorMaps;
import binaryMatrixFile.BinaryMatrixFile;
import oneLiners.OneLiners;
import seed.optimization.HookeAndJeeves;


/** Simple imaging by a single lens. Also produces the PSF interpolation data for the imaging system.
 * Tests the basic optics and PSF collection and evetually the basic polarisation transfer.
 * 
 * After running this, PSFImaging can be used to draw images faster.
 * 
 * @author oliford
 *
 */
public class GeneratePSFInterpolationExample {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/imgTest";
	
	final static double maxTheta = 30 * Math.PI / 180;
	final static int nRaysPerSource = 1000;
	final static double rt2 = Math.sqrt(2);
	
	/** Grid definition for interpolation region { {x0, x1, nX}, {y0, y1, nY}, { z0, z1, nZ } } */
	final static double gridDef[][] = { 
			{ -0.300, 0.300, 10 }, 
			{ -0.300, 0.300, 10 }, 
			{ -0.300, 0.300, 10 }
		};
	
	/** Region on image plane */ 
	final static double imageX[] = OneLiners.linSpace(-0.3, 0.3, 500);
	final static double imageY[] = OneLiners.linSpace(-0.3, 0.3, 500);
	
	final static double wavelen = 593e-9;
	
	public static void main(String[] args) {
		
		double f = 0.050;	//lens focal length
		double u = 0.600;	//object distance
		double v = 1.0 / (1.0/f - 1.0/u);	//ideal image distance
		double autoFocus = -0.00125885009765625; //set to 0 to run autofocus, then copy the best value in here
		
		Material lensMat = new IsotropicFixedIndexGlass(1.5);
		Medium lensMed = new Medium(lensMat);
		
		Square backPlane = new Square("backPlane", new double[]{ -0.3, 0, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.100, 0.300, Absorber.ideal());
		//Iris lensIris = new Iris("lensIris", new double[]{ u, 0, 0 }, new double[]{ -1, 0, 0 }, 0.500, 0.030, null, null, Absorber.ideal());
				
		DoubleGaussLens lens = new DoubleGaussLens(); //focal length = 0.1
		lens.shift(new double[]{ u + autoFocus, 0, 0 });
		lens.scale(0.5);
		
		Square imgPlane = new Square("imgPlane", new double[]{ u+v, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 1.100, 1.300, Absorber.ideal());
		
		PointSpreadBuilder psfBuild = new PointSpreadBuilder(imgPlane, outPath + "/psfData.bin");
		PSFGrid psfGrid = new PSFGrid("imgTest", gridDef, DualGaussianPSF.class);
		
		Optic all = new Optic("all", new Element[]{ backPlane, lens, /*lensIris,*/ imgPlane });
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/imgTest.vrml", 0.01);
		if(vrmlOut != null) {
			vrmlOut.setDrawPolarisationFrames(false);
			vrmlOut.setSkipRays(97);
		}
		
		if(autoFocus == 0){
			AutoFocusOld af = new AutoFocusOld();
			af.setTracingElements(all, imgPlane);
			af.initRaysImaging(lens.getSurfaces().get(0), new double[]{0,0,0}, new double[]{0,0,0.002}, 10, wavelen, 10000);
			//af.initRaysParallel(lens.getSurfaces().get(0), new double[]{1,0,0}, 10*Math.PI/180, 0.10, wavelen, 20, 200);
			af.setMovement(lens, new double[]{1,0,0}, -0.010, 0.010);
			//af.setSVGOut(outPath + "/autoFocus", 1);
			//af.setHitsDebugFile(outPath + "/autoFocusHits.bin");
			af.optimise(new HookeAndJeeves(null), 100);			
		}
				
		double col[][] = ColorMaps.alternating2D2x2(psfGrid.getNY(), psfGrid.getNZ());
		
		for(int iZ=0; iZ < psfGrid.getNZ(); iZ++) {
			for(int iY=0; iY < psfGrid.getNY(); iY++) {
				for(int iX=0; iX < psfGrid.getNX(); iX++) {
					double startPos[] = psfGrid.gridPos(iX, iY, iZ);
					
					psfBuild.startNewPSF(startPos, new DualGaussianPSF(20));
					
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
		
		
		if(vrmlOut != null) {
			vrmlOut.drawOptic(all);
			vrmlOut.destroy();
		}	
	}	
}
