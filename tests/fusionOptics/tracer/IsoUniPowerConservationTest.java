package fusionOptics.tracer;

import java.util.List;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.IsoUniaxialInterface;
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
import junit.framework.TestCase;

/** Check that total intensity is preserved for light on iso/iso and iso/uni interfaces */
public class IsoUniPowerConservationTest extends TestCase {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/powerCheck";
	
	final static double maxTheta = 30 * Math.PI / 180;
	final static int nRays = 1;
	final static double rtHalf = Math.sqrt(0.5);
	final static double rt2 = Math.sqrt(2);
	
	final static double wavelen = 593e-9;
	
	public void testPowerConservations(){
		
		double tol = 1e-4;
		
		Material glassMatIso = new IsotropicFixedIndexGlass(1.5);
		Medium glassMedIso = new Medium(glassMatIso);
		
		Material glassMatUni = new CrystalQuartz();
		Medium glassMedUni = new Medium(glassMatUni);
		
		Square backPlane = new Square("backPlane", new double[]{ 0, 0, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 1, 0.2, Absorber.ideal());
		Square fwdPlane = new Square("fwdPlane", new double[]{ 2, 0, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 1, 0.2,  Absorber.ideal());
		
		Square surfA1 = new Square("surfA1", new double[]{ 0.8, -0.25, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.4, 0.2, null, glassMedIso, IsoIsoStdFresnel.ideal());
		Square surfA2 = new Square("surfA2", new double[]{ 1.2, -0.25, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.4, 0.2, null, glassMedIso, IsoIsoStdFresnel.ideal());
		
		Square surfB1 = new Square("surfB1", new double[]{ 0.8, 0.25, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.4, 0.2, null, glassMedUni, IsoUniaxialInterface.ideal());
		Square surfB2 = new Square("surfB2", new double[]{ 1.2, 0.25, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.4, 0.2, null, glassMedUni, IsoUniaxialInterface.ideal());
		
		Optic all = new Optic("all", new Element[]{ backPlane, surfA1, surfA2, surfB1, surfB2, fwdPlane });
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/powerCheck.vrml", 0.01);
		vrmlOut.setDrawPolarisationFrames(true);
		vrmlOut.setDrawOnlyStrongest(false);
		
		double col[][] = new double[][]{ { 1, 1, 0  } };//ColorMaps.jet(nRays);
		//double out[][] = new double[nRays][3];
		
		//for(int i=0; i < nRays; i++) {
		
		/* Isotropic at normal incidence */
		RaySegment ray = new RaySegment();
		ray.startPos = new double[]{ 0.2, -0.25, 0 };
		ray.dir = new double[]{ 1, 0, 0 };			
		ray.length = Double.POSITIVE_INFINITY;
		ray.up = new double[]{0,0,1};			
		ray.E0 = new double[][]{ {rtHalf,0,rtHalf,0 } };
		
		ray.wavelength = wavelen;
					
		Tracer.trace(all, ray, 500, 0.001, true);
		vrmlOut.drawRay(ray, col[0]);
		ray.dumpPath();
		
		double I0 = Pol.intensity(ray.E1);
		double Ir1 = Pol.intensity(ray.endHit.reflectedOrdinary.E0);
		double It1 = Pol.intensity(ray.endHit.transmittedOrdinary.E0);
		double Ir2 = Pol.intensity(ray.endHit.transmittedOrdinary.endHit.reflectedOrdinary.E0);
		double It2 = Pol.intensity(ray.endHit.transmittedOrdinary.endHit.transmittedOrdinary.E0);

		System.out.println(I0);
		System.out.println(It1 + "\t" + Ir1 + "\t" + (It1 + Ir1));
		System.out.println(It2 + "\t" + Ir2 + "\t" + (It2 + Ir2));
		
		//refl + transmitted should match at both boundaries
		assertEquals(It1 + Ir1, I0, tol);
		assertEquals(It2 + Ir2, It1, tol);

		
		/* Isotropic at 28' incidence */
		double ang = 28 * Math.PI / 180;
		RaySegment ray2 = new RaySegment();
		ray2.startPos = new double[]{ 0.6, -0.45, 0 };
		ray2.dir = new double[]{ Math.cos(ang), Math.sin(ang), 0 };			
		ray2.length = Double.POSITIVE_INFINITY;
		ray2.up = new double[]{0,0,1};			
		ray2.E0 = new double[][]{ {rtHalf,0,rtHalf,0 } };
		ray2.wavelength = wavelen;
					
		Tracer.trace(all, ray2, 500, 0.001, true);
		vrmlOut.drawRay(ray2, col[0]);
		ray2.dumpPath();
		
		I0 = Pol.intensity(ray2.E1);
		Ir1 = Pol.intensity(ray2.endHit.reflectedOrdinary.E0);
		It1 = Pol.intensity(ray2.endHit.transmittedOrdinary.E0);
		Ir2 = Pol.intensity(ray2.endHit.transmittedOrdinary.endHit.reflectedOrdinary.E0);
		It2 = Pol.intensity(ray2.endHit.transmittedOrdinary.endHit.transmittedOrdinary.E0);

		System.out.println(I0);
		System.out.println(It1 + "\t" + Ir1 + "\t" + (It1 + Ir1));
		System.out.println(It2 + "\t" + Ir2 + "\t" + (It2 + Ir2));
		
		//refl + transmitted should match at both boundaries
		assertEquals(It1 + Ir1, I0, tol);
		assertEquals(It2 + Ir2, It1, tol);
		
		
		/* Birefringent at normal incidence */
		RaySegment ray3 = new RaySegment();
		ray3.startPos = new double[]{ 0.2, 0.25, 0 };
		ray3.dir = new double[]{ 1, 0, 0 };
		ray3.length = Double.POSITIVE_INFINITY;
		ray3.up = new double[]{ 0, 0, 1 };
		ray3.E0 = new double[][]{ {rtHalf,0,rtHalf,0 } };
		ray3.wavelength = wavelen;

		Tracer.trace(all, ray3, 500, 0.001, true);
		vrmlOut.drawRay(ray3, col[0]);
		ray3.dumpPath();			
		
		I0 = Pol.intensity(ray3.E1);
		Ir1 = Pol.intensity(ray3.endHit.reflectedOrdinary.E0);
		double It1o = Pol.intensity(ray3.endHit.transmittedOrdinary.E0);
		double It1e = Pol.intensity(ray3.endHit.transmittedExtraordinary.E0);
		double It2o = Pol.intensity(ray3.endHit.transmittedOrdinary.endHit.transmittedOrdinary.E0);
		double It2e = Pol.intensity(ray3.endHit.transmittedExtraordinary.endHit.transmittedOrdinary.E0);

		System.out.println(I0);
		System.out.println(It1o + "\t" + It1e + "\t" + Ir1 + "\t" + (Ir1 + It1o + It1e));
		System.out.println(It2o + "\t" + It2e + "\t" + (It2o + It2e) + "\t" + (It1e + It1o));

		assertEquals(Ir1 + It1o + It1e, I0, tol); //refl + transmittedE+O should match at first boundary
		assertEquals(It2o + It2e, It1o + It1e, tol);
		

		/* Birefringent at angled incidence */
		ang = 28 * Math.PI / 180;
		RaySegment ray4 = new RaySegment();
		ray4.startPos = new double[]{ 0.6, 0.05, 0 };
		ray4.dir = new double[]{ Math.cos(ang), Math.sin(ang), 0 };
		ray4.length = Double.POSITIVE_INFINITY;
		ray4.up = new double[]{ 0, 0, 1 };
		ray4.E0 = new double[][]{ {rtHalf,0,rtHalf,0 } };
		ray4.wavelength = wavelen;

		Tracer.trace(all, ray4, 500, 0.001, true);
		vrmlOut.drawRay(ray4, col[0]);
		ray4.dumpPath();			
		
		I0 = Pol.intensity(ray4.E1);
		Ir1 = Pol.intensity(ray4.endHit.reflectedOrdinary.E0);
		It1o = Pol.intensity(ray4.endHit.transmittedOrdinary.E0);
		It1e = Pol.intensity(ray4.endHit.transmittedExtraordinary.E0);
		It2o = Pol.intensity(ray4.endHit.transmittedOrdinary.endHit.transmittedOrdinary.E0);
		It2e = Pol.intensity(ray4.endHit.transmittedExtraordinary.endHit.transmittedOrdinary.E0);

		System.out.println(I0);
		System.out.println(It1o + "\t" + It1e + "\t" + Ir1 + "\t" + (Ir1 + It1o + It1e));
		System.out.println(It2o + "\t" + It2e + "\t" + (It2o + It2e) + "\t" + (It1e + It1o));

		assertEquals(Ir1 + It1o + It1e, I0, tol); //refl + transmittedE+O should match at first boundary
		assertEquals(It2o + It2e, It1o + It1e, tol);
		
		//}
			
		vrmlOut.drawOptic(all);
		vrmlOut.destroy();
		
	}
	
	
}
