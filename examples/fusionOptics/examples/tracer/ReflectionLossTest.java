package fusionOptics.tracer;

import java.util.List;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.Reflector;
import fusionOptics.materials.BK7;
import fusionOptics.materials.Calcite;
import fusionOptics.materials.CrystalQuartz;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.Box;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;

import otherSupport.ColorMaps;
import otherSupport.RandomManager;

import binaryMatrixFile.BinaryMatrixFile;

import oneLiners.OneLiners;
import net.jafama.FastMath;

/** Examine transmittance and reflectance vs angle and polarisation for a single interface */
public class ReflectionLossTest {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/reflLoss";
	
	final static double maxTheta = 30 * Math.PI / 180;
	final static int nRays = 1;
	final static double rtHalf = Math.sqrt(0.5);
	final static double rt2 = Math.sqrt(2);
	
	final static double wavelen = 593e-9;
	
	public static void main(String[] args) {
		
		Material glassMat = new IsotropicFixedIndexGlass(1.5);
		Medium glassMed = new Medium(glassMat);
		
		Square backPlane = new Square("backPlane", new double[]{ 0, 0, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 1, 0.2, Reflector.ideal());
		Square fwdPlane = new Square("fwdPlane", new double[]{ 1, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 1, 0.2, Reflector.ideal());
		
		Square surf1 = new Square("surf1", new double[]{ 0.51, 0, 0 }, new double[]{ rtHalf, -rtHalf, 0 }, new double[]{ rtHalf, rtHalf, 0 }, rt2, 0.2, null, glassMed, IsoIsoStdFresnel.ideal());
		Square surf2 = new Square("surf2", new double[]{ 0.49, 0, 0 }, new double[]{ -rtHalf, rtHalf, 0 }, new double[]{ rtHalf, rtHalf, 0 }, rt2, 0.2, null, glassMed, IsoIsoStdFresnel.ideal());
		
		//a block with no refraction and no interface loss, but very low bulk transmission
		Medium blockMed = new Medium(new IsotropicFixedIndexGlass(1.0, 0.1));
		Box block = new Box("block", new double[]{ 0.2, 0.2, -0.1}, new double[]{ 0.5, 0.45, 0.1 }, blockMed, IsoIsoInterface.ideal());
			
		Optic all = new Optic("all", new Element[]{ backPlane, surf1, surf2, block, fwdPlane });
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/reflLoss.vrml", 0.01);
		vrmlOut.setDrawPolarisationFrames(true);
		vrmlOut.setDrawOnlyStrongest(false);
		
		double col[][] = new double[][]{ { 1, 1, 0  } };//ColorMaps.jet(nRays);
		
		//double out[][] = new double[nRays][3];
		
		for(int i=0; i < nRays; i++) {
			
			RaySegment ray = new RaySegment();
			ray.startPos = new double[]{ -0.2, -0.62, 0 };
			ray.dir = Util.reNorm(new double[]{ 0.2, 0.03, 0 });//*/
			
			ray.length = Double.POSITIVE_INFINITY;
			ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, new double[]{0,0,1})), ray.dir);
			
			ray.E0 = new double[][]{ {rtHalf/2,0,rtHalf/2,0 }, {rtHalf/2,0,-rtHalf/2,0} };
			
			ray.rotatePolRefFrame(new double[]{ //0,Math.sqrt(0.5),Math.sqrt(0.5) 
					RandomManager.instance().nextUniform(0.1, 1.0),
					RandomManager.instance().nextUniform(0.1, 1.0),
					RandomManager.instance().nextUniform(0.1, 1.0),					
			});
			
			
			//ray.rotatePolRefFrame(new double[]{ 0,0,1 });
			
			ray.wavelength = wavelen;
			
			Tracer.trace(all, ray, 500, 0.001, true);
			
			//ray.endHit.transmittedOrdinary.rotatePolRefFrame(new double[]{ 0, 0, 1});
			
			
			//in general, ray.up might be any direction, but when there is a surface intersection, it will be left as the 's' direction
			//I know this, and I'm feeling lazy
			
			/*List<Intersection> frontHits = ray.getIntersections(fwdPlane);
			if(frontHits.size() > 0){
				RaySegment tRay = frontHits.get(0).incidentRay;
				out[i][1] = (tRay.stokes[0] + tRay.stokes[1]) / 2; //E0_up^2 (Ts)
				out[i][2] = (tRay.stokes[0] - tRay.stokes[1]) / 2; //E0_right^2 (Tp)
			}*/
			
			//What to expect:
			//Angles here are clockwise from 'up' looking in the ray travel direction
			// Starting with +45' ( s = [ 1 0 1 0 ] )
			// 
			// Transmitted should ALWAYS be the same - none of the transmission ampltiude coefficients are -ve 
			// 
			// For reflected:
			//   ni < nt (glass entry), ang < brewster: 			Rs=-ve, Rp=+ve ==> Flipped -45'
			//   ni < nt (glass entry), ang > brewster: 			Rs=-ve, Rp=-ve ==> Same +45'
			//   ni > nt (glass exit), ang < brewster: 				Rs=+ve, Rp=-ve ==> Flipped -45'
			//   ni > nt (glass exit), brewster < ang < critical: 	Rs=+ve, Rp=+ve ==> Same +45'
			//
			
			
			vrmlOut.drawRay(ray, col[i]);
		}
		vrmlOut.drawOptic(all);
		vrmlOut.destroy();
		
		//BinaryMatrixFile.mustWrite(outPath + "/stokes.bin", out, false);
	}
	
	
}
