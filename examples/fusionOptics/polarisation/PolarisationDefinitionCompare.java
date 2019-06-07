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

import algorithmrepository.Algorithms;
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
public class PolarisationDefinitionCompare {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/polDefs";
	
	final static double rt2 = Math.sqrt(2);	
	final static double wavelen = 593e-9;
	final static int nPhi = 20;
	final static double theta = 65 * Math.PI / 180;
	
	final static double globalUp[] = new double[]{0,0,1};
	
	public static void main(String[] args) {
		
		double r = 1.3;
		double u = 1.0;
		double f = 1.0;
		double v = (u == f) ? 1.0 : 1 / (1/f - 1/u);
		
		Medium lensGlass = new Medium(new SchottSFL6());
		Interface lensIFace = IsoIsoInterface.ideal();
		
		u = u - 0.2;
		double d = 0.4;
		SimpleDoubleConvexLens lens1 = SimpleDoubleConvexLens.fromFocalLengthAndEdgeThickness("lens1", new double[]{ u - d, 0, 0 }, new double[]{ -1, 0, 0 }, r, 3*f, 0.005, lensGlass, lensIFace, wavelen);
		Iris iris1 = new Iris("iris1", new double[]{ u - d, 0, 0 }, new double[]{ -1, 0, 0 }, 1.4*r, 0.99*r, Absorber.ideal());
		
		SimpleDoubleConvexLens lens2 = SimpleDoubleConvexLens.fromFocalLengthAndEdgeThickness("lens2", new double[]{ u, 0, 0 }, new double[]{ -1, 0, 0 }, r, 3*f, 0.005, lensGlass, lensIFace, wavelen);
		Iris iris2 = new Iris("iris2", new double[]{ u, 0, 0 }, new double[]{ -1, 0, 0 }, 1.4*r, 0.99*r, Absorber.ideal());
		
		SimpleDoubleConvexLens lens3 = SimpleDoubleConvexLens.fromFocalLengthAndEdgeThickness("lens3", new double[]{ u + d, 0, 0 }, new double[]{ -1, 0, 0 }, r, 3*f, 0.005, lensGlass, lensIFace, wavelen);
		Iris iris3 = new Iris("iris3", new double[]{ u + d, 0, 0 }, new double[]{ -1, 0, 0 }, 1.4*r, 0.99*r, Absorber.ideal());
				
		Square imgPlane = new Square("fwdPlane", new double[]{ (u+v), 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 3, 3, Absorber.ideal());		
		
		Disc midLens = new Disc("midLens", new double[]{ u, 0, 0 }, new double[]{ -1, 0, 0 }, r, NullInterface.ideal());

		Optic all = new Optic("all", new Element[]{ lens1, iris1, lens2, iris2, lens3, iris3, midLens, imgPlane });
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/polDefs.vrml", 0.01);
		vrmlOut.setDrawPolarisationFrames(true);
		
		double col[][] = ColorMaps.jet(nPhi);
		
		for(int i=0; i < nPhi; i++) {
			RaySegment ray = new RaySegment();
			
			double phi = i * FastMath.PI / nPhi; 
						
			ray.startPos = new double[]{ 0, 0, 0 };
			ray.length = Double.POSITIVE_INFINITY;
			//ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, lens);
			ray.dir = new double[]{ 
					FastMath.cos(theta), 
					FastMath.sin(theta) * FastMath.cos(phi),
					FastMath.sin(theta) * FastMath.sin(phi),
				};
			
			double M[][] = Algorithms.rotationMatrix(Util.reNorm(new double[]{0, 0.1, 0.2}), 5 * Math.PI / 180);
			ray.dir = Algorithms.rotateVector(M, ray.dir);
			
			//if( (i % 2) == 0){
				
			//}else{
				///ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, globalUp)), ray.dir);
			//}
				
			ray.up = Tracer.generateRayFanConsistentPolarisationDefinition(ray.startPos, ray.dir, lens1.getBoundarySphereCentre(), globalUp);			
			ray.E0 = new double[][]{{1,0,0,0}, {0,0,0,0}}; //PointSpreadFunction.getInputStatesForMuellerCalc();
			
			ray.rotatePolRefFrame(globalUp);
			ray.E0[1] = new double[]{1,0,0,0};
			
			ray.wavelength = wavelen;
			
			Tracer.trace(all, ray, 1000, 0.1, false);
			
			vrmlOut.drawRay(ray, col[i % 2]);
			
			
		}

		vrmlOut.drawOptic(all);
		vrmlOut.destroy();
		
	}
	
	
}
