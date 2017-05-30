package fusionOptics.polarisation;

import static org.junit.Assert.*;

import org.junit.Test;

import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.interfaces.Reflector;
import fusionOptics.materials.BK7;
import fusionOptics.materials.FusedSilica;
import fusionOptics.media.MediumInMagneticField;
import fusionOptics.optics.Box;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Material;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;


public class TestFaradayRotation {

	@Test
	public void test() {
		
		double length = 0.10; //m
		double Bmag = 1; //Tesla 
		
		double B[] = new double[]{ Bmag, 0, 0 };
		//Material mat = new FusedSilica();
		Material mat = new BK7();
		MediumInMagneticField med = new MediumInMagneticField(mat, null, 300, B);
		
		Box box = new Box("faradayGlass", new double[]{ 1.5, 0, 0 }, length, 0.2, 0.2, med, IsoIsoInterface.ideal());
		
		Square target = new Square("target", new double[]{3,0,0}, new double[]{1,0,0}, new double[]{0,0,1}, 0.2, 0.2, null, null, Reflector.ideal());
		Square target2 = new Square("target2", new double[]{0,0,0}, new double[]{1,0,0}, new double[]{0,0,1}, 0.2, 0.2, null, null, Absorber.ideal());
		
		Optic all = new Optic("all", new Element[]{ box, target, target2 });
	
		RaySegment ray = new RaySegment();
		ray.startPos = new double[]{ 0.1, -0.02, 0 };
		ray.dir = Util.reNorm(new double[]{ 1, 0.01, 0 });
		ray.up = new double[]{ 0, 0, 1 };
		ray.E0 = new double[][]{{ 1, 0, 0 ,0 }};
		ray.length = Double.POSITIVE_INFINITY;
		ray.wavelength = 600e-9;
		
		Tracer.trace(all, ray, 100, 0.001, false);
		
		Intersection hit = ray.getIntersections(target2).get(0);		
		
		VRMLDrawer.dumpRay("/tmp/faradayRay.vrml", all, ray, 0.01);
		
		//This is read with target2's normal in the +ve X direction, not in the ray
		//direction of the ray that hits it, so it's polariastion angle is clockwise looking along +ve x
		double rotation = Pol.psi( Pol.projectToPlanesView(hit, false)[0]);		
		
		//For BK7, the graph in
		// [ E.Munin "Faraday Effect And Energy Gap In Optical Matierals", J.Phys D (1992) ]
		// shows n*V = +ve 7.5 rad / T / m at 600nm		
		double n600 = mat.getRefractiveIndex(0, 600e-9, 300);
		double Vtest = 7.5 / n600;
		
		//we expected it to double, because we pass the glass twice in opposite directions.
		// From the ray's point of view, the field direction and rotation direction
		// both change, but in real 3D space, the rotation should stay clockwise looking along B
		// for any +ve V.
		double expectedRotation = Vtest * Bmag * length * 2; 
		
		System.out.println("Verdet(BK7) = " + mat.getVerdetConstant(0, ray.wavelength, 300)*180/Math.PI + " ° /T/m");
		System.out.println("Vtest(BK7) = " + Vtest *180/Math.PI + " °/T/m");
		System.out.println("Rotation = " + rotation*180/Math.PI + "°, expected = " + expectedRotation*180/Math.PI );
		
		 //5% is reasonable considering it's old experimental data and the variation in different glasses
		assertEquals(expectedRotation, rotation, 0.05 * expectedRotation);
		
		FusedSilica fusedSilica = new FusedSilica();
		n600 = fusedSilica.getRefractiveIndex(0, 600e-9, 300);
		Vtest = 6.0 / n600;
		System.out.println("Verdet(Fused Silica) = " + fusedSilica.getVerdetConstant(0, 600e-9, 300)*180/Math.PI + " ° /T/m");
		System.out.println("Vtest(Fused Silica)  = " + Vtest * 180/Math.PI + " °/T/m");
		
		ray.dumpPath();
		
	}

}
