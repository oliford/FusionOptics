package fusionOptics.examples.polarisation;

import java.text.DecimalFormat;
import java.util.List;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.collection.IntersectionProcessor;
import fusionOptics.collection.LensDepolarisationInfoCollector;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoAntiReflective;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.IsoUniaxialInterface;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.interfaces.Reflector;
import fusionOptics.interfaces.SimplePolariser;
import fusionOptics.materials.BK7;
import fusionOptics.materials.Calcite;
import fusionOptics.materials.CrystalQuartz;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.materials.LithiumNiobate;
import fusionOptics.materials.SchottSFL6;
import fusionOptics.optics.Box;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.pointSpread.DualGaussianPSF;
import fusionOptics.pointSpread.GaussianPSF;
import fusionOptics.pointSpread.PointSpreadBuilder;
import fusionOptics.pointSpread.PointSpreadFunction;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;

import binaryMatrixFile.BinaryMatrixFile;
import binaryMatrixFile.BinaryMatrixWriter;

import otherSupport.ColorMaps;
import otherSupport.RandomManager;
import otherSupport.StatusOutput;
import oneLiners.OneLiners;
import net.jafama.FastMath;

/**
 * Investigation of depolarision / polarisation modification effect
 * due to a lens.
 * 
 * @author oliford
 */
public class LensDepolarisation {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/lensDepol";
	
	final static int maxHits = 10000;
	final static int maxRays = 100000;
	final static double maxTheta = 30 * Math.PI / 180;
	final static double rt2 = Math.sqrt(2);
	
	final static double wavelen = 593e-9;
	
	public static void main(String[] args) {
		
		double r = 0.5;
		double f = 1.0;
		double u = 1.0;
		double v = (u == f) ? 1.0 : 1 / (1/f - 1/u);
		
		//Lens lens = new Lens("lens", new double[]{ u, 0, 0 }, new double[]{ -1, 0, 0 }, r, new Medium(lensMat), IsoIsoStdFresnel.ideal(), f, wavelen);
		Medium lensGlass = new Medium(new SchottSFL6());
		//Medium lensGlass = new Medium(new Quartz(), new double[][]{ { rt2/2, 0, rt2/2 } }, 300);
		//IsotropicFixedIndexGlass lensGlass = new IsotropicFixedIndexGlass(1.805);
		Interface lensIFace = IsoIsoInterface.ideal();
		//Interface lensIFace = IsoIsoStdFresnel.ideal();
		//Interface lensIFace = new IsoIsoAntiReflective(1.4, wavelen/4/1.4, 0);
		//Interface lensIFace = IsoUniaxialInterface.ideal();
		
		SimpleDoubleConvexLens lens = SimpleDoubleConvexLens.fromFocalLengthAndEdgeThickness("lens", new double[]{ u, 0, 0 }, new double[]{ -1, 0, 0 }, r, f, 0.005, lensGlass, lensIFace, wavelen);
		Iris iris = new Iris("iris", new double[]{ u, 0, 0 }, new double[]{ -1, 0, 0 }, 1.4*r, 0.99*r, Absorber.ideal());
		Square imgPlane = new Square("fwdPlane", new double[]{ (u+v), 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 1, Absorber.ideal());		
		
		Disc midLens = new Disc("midLens", new double[]{ u, 0, 0 }, new double[]{ -1, 0, 0 }, r, NullInterface.ideal());

		/*
		 * //L1 from __WRONG__ AUG tube optics
		Iris iris = new Iris("lens1Iris", new double[]{ 0.2595-0.2595+u, 0,0 }, new double[]{ -1, 0, 0 }, 0.062, 0.025, Absorber.ideal());
		Dish lens1Front = new Dish("lens1Front", new double[]{ 0.2595-0.2595+u, 0,0 }, new double[]{ 1, 0, 0 }, 0.311, 0.025,  lensGlass, null, lensIFace);
		Disc lens1Back = new Disc("lens1Back", new double[]{ 0.264-0.2595+u, 0,0 }, new double[]{ 1, 0, 0}, 0.025, null, lensGlass, lensIFace);
		Optic lens = new Optic("lens", new Element[]{ lens1Front, lens1Back });
		*/
		
		Optic all = new Optic("all", new Element[]{ lens, iris, midLens, imgPlane });
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/lensDepol.vrml", 0.01);
		vrmlOut.setDrawPolarisationFrames(true);
		vrmlOut.setSkipRays(maxHits <= 1 ? 0 : (maxHits/99));
		
		double col[][] = ColorMaps.jet(maxRays);
		
		LensDepolarisationInfoCollector polCheck = new LensDepolarisationInfoCollector(midLens, imgPlane, 50, 0, outPath + "/polInf");
		
		PointSpreadBuilder psfBuild = new PointSpreadBuilder(imgPlane);
		psfBuild.startNewPSF(new double[]{ 0, 0, 0 }, new DualGaussianPSF());
		
		psfBuild.setMaxCoherentIntegrationRadius(0.001);
	
		StatusOutput stat = new StatusOutput(LensDepolarisation.class, maxHits);
		for(int i=0; i < maxRays; i++) {
			RaySegment ray = new RaySegment();
			
			ray.startPos = new double[]{ 0, 0, 0 };
			ray.length = Double.POSITIVE_INFINITY;
			ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, lens);
			
			ray.up = Tracer.generateRayFanConsistentPolarisationDefinition(ray.startPos, ray.dir, lens.getBoundarySphereCentre(), new double[]{0,0,1});
			ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, new double[]{0,0,1})), ray.dir);
			
			ray.E0 = PointSpreadFunction.getInputStatesForMuellerCalc();
			
			ray.wavelength = wavelen;
			
			Tracer.trace(all, ray, 1000, 0.1, false);
			
			vrmlOut.drawRay(ray, col[i]);
			
			ray.processIntersections(imgPlane, psfBuild, polCheck);
			psfBuild.nextCoherentSet();
			
			List<Intersection> hits = ray.getIntersections(imgPlane);
			
			if(hits.size() > 0){
				Intersection hit = hits.get(0); //there should only be 1
				
				double theta = FastMath.acos(ray.dir[0]);
				double phi = FastMath.atan2(ray.dir[1], ray.dir[2]);
				
				//put final ray back in sense of 'z'
				hit.incidentRay.rotatePolRefFrame(new double[]{ 0, 0, 1 });
			}
			
			
			Pol.recoverAll();
			
			int n = psfBuild.getNPointsCollected();
			stat.doStatus(n);
			if(n >= maxHits)
				break;
		}
		stat.done();

		DualGaussianPSF psf = (DualGaussianPSF) psfBuild.psfDone(false);
		System.out.println("Muller: ");
		psf.dumpNormalisedMueller();

		vrmlOut.drawOptic(all);
		vrmlOut.destroy();
		
		polCheck.write();
	}
	
	
}
