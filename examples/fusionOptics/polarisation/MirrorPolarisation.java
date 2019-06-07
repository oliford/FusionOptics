package fusionOptics.polarisation;

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
import fusionOptics.polarisation.MullerTest;
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
 * Investigation of polarisation modification due to mirrors.
 * 
 * Basically, there isn't any. It took me a long time to be confident in this conclusion and to
 * finally fix the Relfection interface so that it works. See the note in Reflector.pureReflection().
 * 
 * Two days wasted... grumble.
 * 
 * @author oliford
 */
public class MirrorPolarisation {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/mirrorPol";
	
	final static int maxHits = 10000;
	final static int maxRays = 10000;
	final static double maxTheta = 30 * Math.PI / 180;
	final static double rt2 = Math.sqrt(2);
	
	final static double wavelen = 593e-9;
	
	public static void main(String[] args) {
		Square mirror = new Square("mirror", new double[]{ 2, 0, 0 }, 
				new double[]{ -1, 0, 0 }, 
				new double[]{ 0, 0, 1 }, 1, 1, Reflector.ideal());
		
		Util.rotateOnZ(mirror, mirror.getCentre(), -45*Math.PI/180);
		
		//Two interrogation planes, one before the mirror and one after.
		Square polPlane1 = new Square("polPlane1", new double[]{ 1, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1.3, 1.3, NullInterface.ideal());
		Square polPlane2 = new Square("polPlane2", new double[]{ 2, 2, 0 }, new double[]{ 0, -1, 0 }, new double[]{ 0, 0, 1 }, 5, 5, Absorber.ideal());
		
		// An ideal collimating lens, which should have no effect if the projection is working properly.
		SimpleDoubleConvexLens collLens = SimpleDoubleConvexLens.fromFocalLengthAndEdgeThickness("collLens", new double[]{ 2, 1, 0 }, new double[]{ 0, -1, 0 }, 
				1.3, //radius
				//IsoIsoStdFresnel.ideal(),
				//new IsoIsoAntiReflective(1.25, wavelen/4/1.25, 0), 
				3.0, //focal length
				0.005, // edge thickness
				new Medium(new IsotropicFixedIndexGlass(1.5)), IsoIsoInterface.ideal(), wavelen);

		Disc shield = new Disc("shield", new double[]{ 0.2, 0.3, 0 }, new double[]{ 0, 1, 0 }, 1.3, Absorber.ideal());
		
		Util.rotateOnX(polPlane1, polPlane1.getCentre(), 17*Math.PI/180);
		Util.rotateOnY(polPlane2, polPlane2.getCentre(), 20*Math.PI/180);
		//Util.rotateOnZ(imgPlane, imgPlane.getCentre(), 30*Math.PI/180);
		
		Optic all = new Optic("all", new Element[]{
				collLens,
				shield,
				polPlane1,
				mirror, polPlane2  });
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/mirrorPol.vrml", 0.02);
		vrmlOut.setDrawPolarisationFrames(true);
		vrmlOut.setSkipRays(maxHits <= 1 ? 0 : (maxHits/99));
		
		double col[][] = ColorMaps.jet(maxRays);
		
		int nPol = PointSpreadFunction.inputStokesForMuellerCalc.length;
		IntersectionProcessor polCheck[] = new LensDepolarisationInfoCollector[nPol*2]; 
		for(int i=0; i < nPol; i++){
			polCheck[i] = new LensDepolarisationInfoCollector(polPlane1, polPlane1, 50, i, outPath + "/polInf1"+i);
			polCheck[nPol+i] = new LensDepolarisationInfoCollector(polPlane2, polPlane2, 50, i, outPath + "/polInf2"+i);
		}
		
		PointSpreadBuilder psfBuild = new PointSpreadBuilder(polPlane2);
		psfBuild.startNewPSF(new double[]{ 0, 0, 0 }, new DualGaussianPSF());
		psfBuild.setMaxCoherentIntegrationRadius(0.001);
		
		//double commonDef[] = new double[]{ 0, 0, 1 };
		//double commonDef[] = new double[]{ 0, 1, 0 };
		double commonDef[] = new double[]{ 0, rt2/2, rt2/2 };
		
		StatusOutput stat = new StatusOutput(MullerTest.class, maxHits);
		for(int i=0; i < maxRays; i++) {
			RaySegment ray = new RaySegment();
			
			ray.startPos = new double[]{ 0, 0, 0 };
			
			ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, mirror);
			ray.up = Tracer.generateRayFanConsistentPolarisationDefinition(ray.startPos,  ray.dir, mirror.getBoundarySphereCentre(), commonDef);
			//ray.dir = Util.reNorm(new double[]{ 1,-0.5,0.2 });
/*			// For drawing the nice picture:
			double theta0 = (-30 + (i*30)) * Math.PI/ 180;
			ray.dir = new double[]{ Math.cos(theta0), 0, Math.sin(theta0) };			
			if(i==3){
				ray.startPos = new double[]{ 2-rt2/2, rt2/2, 0.2 };
				ray.dir = new double[]{ rt2/2, -rt2/2, 0 };
			}//*/
			
			ray.length = Double.POSITIVE_INFINITY;
			
			//the simple common sense of up
			//ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, commonDef)), ray.dir);
			
			ray.E0 = PointSpreadFunction.getInputStatesForMuellerCalc();
			
			ray.wavelength = wavelen;
			
			Tracer.trace(all, ray, 1000, 0.1, false);
			
			vrmlOut.drawRay(ray, col[i]);
			
			ray.processIntersections(polPlane2, psfBuild);
			ray.processIntersections(polPlane2, polCheck);
			psfBuild.nextCoherentSet();
			
			List<Intersection> hits = ray.getIntersections(polPlane1);
			
			if(hits.size() > 0){
				Intersection hit = hits.get(0); //there should only be 1
				
				RaySegment virtRay = new RaySegment();
				
				double EInPlane[][] = Pol.projectToPlanesView(hit, true);
				virtRay.dir = polPlane1.getNormal().clone();
				virtRay.up = polPlane1.getUp().clone();
				//vrmlOut.drawPolarisationFrame(virtRay, hit.pos, EInPlane, 0.5);
				
			}
			
			Pol.recoverAll();
			
			int nDone = psfBuild.getNPointsCollected();
			stat.doStatus(nDone);
			if(nDone >= maxHits)
				break;
		}
		stat.done();
		
		
		System.out.println("nPoints = " + psfBuild.getNPointsCollected());
		DualGaussianPSF psf = (DualGaussianPSF) psfBuild.psfDone(false);
		System.out.println("Normed Muller: ");
		psf.dumpNormalisedMueller();
		System.out.println("Muller: ");
		psf.dumpMueller();
		
		double M[] = psf.getMeanMueller();
		System.out.println("Linear depolarisation(quad add) = " + FastMath.sqrt(M[5]*M[5] + M[6]*M[6])/M[0]); 
		System.out.println("Linear depolarisation(add) = " + (M[5] + M[6])/M[0]); 
		
		vrmlOut.drawOptic(all);
		vrmlOut.destroy();
		
		
		for(int i=0; i < nPol*2; i++)
			((LensDepolarisationInfoCollector)polCheck[i]).write();

	}	
}
