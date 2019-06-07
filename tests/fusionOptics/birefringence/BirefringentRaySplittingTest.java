package fusionOptics.birefringence;

import java.text.DecimalFormat;
import java.util.List;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.IsoUniaxialInterface;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.interfaces.SimplePolariser;
import fusionOptics.materials.BK7;
import fusionOptics.materials.Calcite;
import fusionOptics.materials.CrystalQuartz;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.materials.LithiumNiobate;
import fusionOptics.materials.UniaxialFixedIndexGlass;
import fusionOptics.optics.Box;
import fusionOptics.optics.SimpleDoubleConvexLens;
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

import otherSupport.ColorMaps;
import otherSupport.RandomManager;

import binaryMatrixFile.BinaryMatrixFile;
import binaryMatrixFile.BinaryMatrixWriter;
import oneLiners.OneLiners;
import net.jafama.FastMath;
import junit.framework.TestCase;

/** Tests for splitting of rays by refraction in birefringent materials
 * with non-normal optic axis.
 * 
 * Fires a single ray at each incidence angle ang[] to a glass block
 * with optic axis given in the papers.
 * 
 * Checks:
 * 	1) Ray directions for ordinary and extraordinary rays at 30°, 45° and 60°
 *  2) Energy conservation of fresnel coefficients
 *  
 *  */
public class BirefringentRaySplittingTest extends TestCase {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/uniAxialTest";
	
	final static double maxTheta = 30 * Math.PI / 180;
	final static double rt2 = Math.sqrt(2);
	
	
	//final static double ang[] = new double[]{ 0, 30 * Math.PI/180, 45 * Math.PI/180, 60 * Math.PI/180 };
	final static double ang[] = OneLiners.linSpace(0, 85*Math.PI/180, 5.0*Math.PI/180);
	final static int nRays = ang.length;
	
	/** NB: The test vector in M.Avendano is wrong. They write 0.75, 0.5, 0.433 which 
	 * is copied from G.Beyerle (ref 2 in M.Avendano). However, that paper had
	 * a frame where X was the surface normal, whereas Avendano claims Z is the normal.
	 * 
	 * From reading Beyerle, 0.75 must be the optic axis component normal to the surface.
	 * 0.5 must be the component parallel to the surface and the incidence plane.
	 * 0.433 is the component parallel to the surface, but perp to the incident ray.
	 * 
	 * Our incidence plane here is (x,y) and surface is (y,z). */
	final static double opticAxis[] = new double[]{ 0.75, 0.5, 0.433 }; //to match directions 
	//final static double opticAxis[] = new double[]{  0.433, 0.75, 0.5, }; //to match Fig3 of A.Weidlich for fresnel coeffs
	
	/**
		from Avendano and Beyerle, the results for opticAxis = { 0.75, 0.5, 0.433 } should be:
		 	theta_I		ord. x,y,z				extraord. x,y,z 
	*/
	final static double dirChecks[][][] = new double[][][]{
			{ { 30*Math.PI/180 }, { 0.946133, 0.32378, 0 }, { 0.945516, 0.325546, 0.0044153 } },
			{ { 45*Math.PI/180 }, { 0.889007, 0.457894, 0 }, { 0.888783, 0.458307, 0.004536 } },
			{ { 60*Math.PI/180 }, { 0.827949, 0.560803, 0 }, { 0.828391, 0.560131, 0.00456446 } }
	};
	
	public void testAvendanoFixedCases() {
		final double wavelen = 593e-9;
		
		Material plateMat = new UniaxialFixedIndexGlass(1.54426, 1.55335);
		
		Medium plateMedium = new Medium(plateMat, new double[][]{ opticAxis }, 300);
		//Medium plateMedium = new Medium(plateMat, new double[][]{ { 1, 0, 0 } }, 300);
		
		Square backPlane = new Square("backPlane", new double[]{ 0, 0, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 15, 2, Absorber.ideal());
		Square imgPlane = new Square("fwdPlane", new double[]{ 2, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 15, 2, Absorber.ideal());

		Box plate = new Box("plate", new double[]{ 1, 0, 0 }, 0.1, 2, 1, plateMedium, IsoUniaxialInterface.ideal());

		Optic all = new Optic("all", new Element[]{ backPlane, plate, imgPlane });
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/uniAxialTest.vrml", 0.01);
		vrmlOut.setDrawPolarisationFrames(true);
		
		double col[][] = ColorMaps.jet(nRays);
		
		DecimalFormat fmt = new DecimalFormat("0.000000");
		
		BinaryMatrixWriter fresnelOut = new BinaryMatrixWriter(outPath + "/fresnel.bin", 9); 
		
		boolean failed = false;
		for(int i=0; i < nRays; i++) {
			
			RaySegment ray = new RaySegment();
			ray.dir = Util.reNorm(new double[]{
					Math.cos(ang[i]),
					Math.sin(ang[i]),
					0,
			});//*/
			
			ray.startPos = new double[]{
				plate.getBoundarySphereCentre()[0] - ray.dir[0],
				0 - ray.dir[1],
				0 - ray.dir[2] + ((-nRays/2 + i) * 0.02)
			};
			
			ray.length = Double.POSITIVE_INFINITY;
			ray.up = new double[]{0,0,1};
			
			ray.E0 = new double[][]{ { 1,0,0,0 }, {0,0,1,0}  };
			
			ray.wavelength = wavelen;
			
			Tracer.trace(all, ray, 10, 0.01, true);
			
			List<Intersection> hits = ray.getIntersections(plate.getSurfaces().get(1));
			for(Intersection hit : hits){
				
				int dirChkIdx = -1;
				for(int j=0; j < dirChecks.length; j++){
					if(FastMath.abs(dirChecks[j][0][0] - ang[i]) < 1e-6){
						dirChkIdx = j;
						break;
					}
				}
					
				double Rss=0, Rpp=0, Rsp=0, Rps=0, Tso=0, Tse=0, Tpo=0, Tpe=0;
				
				if(hit.transmittedOrdinary != null){
					double dir[] = hit.transmittedOrdinary.dir;
					Tso = Pol.intensity(hit.transmittedOrdinary.E0[0]);
					Tpo = Pol.intensity(hit.transmittedOrdinary.E0[1]);
					System.out.println(
							"O: " + 	fmt.format(dir[0]) + " " + fmt.format(dir[1]) + " " + fmt.format(dir[2]) + 
							"\tTso=" + fmt.format(Tso) + 
							", Tpo=" + fmt.format(Tpo)	 +
							", d x u = " + fmt.format(Util.dot(dir, hit.transmittedOrdinary.up))
							);
					
					if(dirChkIdx >= 0){
						double d = Util.length(Util.minus(dir, dirChecks[dirChkIdx][1]));
						if(d > 1e-5){
							System.err.println("Direction checks failed for ordinary ray at angle = " + (dirChecks[dirChkIdx][0][0]*180/Math.PI) + "°");
							failed = true;
						}
					}
				}
				if(hit.transmittedExtraordinary != null){
					double dir[] = hit.transmittedExtraordinary.dir;
					Tse=Pol.intensity(hit.transmittedExtraordinary.E0[0]);
					Tpe=Pol.intensity(hit.transmittedExtraordinary.E0[1]);
					System.out.println(
							"E: " + 	fmt.format(dir[0]) + " " + fmt.format(dir[1]) + " " + fmt.format(dir[2]) + 
							"\tTse=" + fmt.format(Tse) + 
							", Tpe=" + fmt.format(Tpe) +
							", d x u = " + fmt.format(Util.dot(dir, hit.transmittedExtraordinary.up))
							);
					
					if(dirChkIdx >= 0){
						double d = Util.length(Util.minus(dir, dirChecks[dirChkIdx][2]));
						if(d > 1e-5){
							System.err.println("Direction checks failed for extraordinary ray at angle = " + (dirChecks[dirChkIdx][0][0]*180/Math.PI) + "°");
							failed = true;
						}
					}
				}
				if(hit.reflectedOrdinary != null){
					double dir[] = hit.reflectedOrdinary.dir;
					double Es[] = hit.reflectedOrdinary.E0[0]; //what was originally s
					double Ep[] = hit.reflectedOrdinary.E0[1]; //what was originally p
					Rss = Es[Pol.uRe]*Es[Pol.uRe] + Es[Pol.uIm]*Es[Pol.uIm];
					Rpp = Ep[Pol.rRe]*Ep[Pol.rRe] + Ep[Pol.rIm]*Ep[Pol.rIm];
					Rsp = Es[Pol.rRe]*Es[Pol.rRe] + Es[Pol.rIm]*Es[Pol.rIm];
					Rps = Ep[Pol.uRe]*Ep[Pol.uRe] + Ep[Pol.uIm]*Ep[Pol.uIm];
					System.out.println(
							"R: " + 	fmt.format(dir[0]) + " " + fmt.format(dir[1]) + " " + fmt.format(dir[2]) + 
							", Rss=" + fmt.format(Rss) +
							", Rpp=" + fmt.format(Rpp) +
							", Rsp=" + fmt.format(Rsp) +
							", Rps=" + fmt.format(Rps) +
							", d x u = " + fmt.format(Util.dot(dir, hit.reflectedOrdinary.up))
							);
				}
				
				if(1.0 - (Rss+Rsp+Tso+Tse) > 1e-4){
					System.err.println("Energy conservation broken for s component in fresnel transfer coefficients");
					failed = true;
				}
				if(1.0 - (Rpp+Rps+Tpo+Tpe) > 1e-4){
					System.err.println("Energy conservation broken for p component in fresnel transfer coefficients");
					failed = true;
				}
					
				fresnelOut.writeRow(ang[i], Rss, Rpp, Rsp, Rps, Tso, Tse, Tpo, Tpe);
			}
			
			vrmlOut.drawRay(ray, col[i]);
		}
		
		vrmlOut.drawOptic(all);
		vrmlOut.destroy();
		
		fresnelOut.close();
		
		//we can't actually fail until the end, because 
		//we need to finish writing the files etc
		assertTrue("Failed. Check stderr messages.", !failed);
	}
	
	/** Check that polarisations don't magically change direction(phase sign) on
	 * entry to anisotropic media. This was fixed by fiddling with the signs in 
	 * {@link IsoUniaxialInterface}.calcFresnelCoeffs() */
	public void testEntryExitPhases(){
		final double wavelen = 0.6;
		final double rt2r = 1 / Math.sqrt(2);
		
		double iAng = 12 * Math.PI / 180;
		
		Material mat = new UniaxialFixedIndexGlass(1.0, 1.2);
		
		Medium medium = new Medium(mat, new double[][]{ { 0, 0, 1 } }, 300); //this gets set later
		
		double fw = Util.calcWaveplateFullWaveThickness(mat, wavelen);
		System.out.println("fw = " + fw);
		fw = fw * FastMath.cos(iAng);
		
		Square squares[] = new Square[] {
			new Square("plate0", new double[]{ 1.0, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 10.0, null, medium, IsoUniaxialInterface.ideal()),
			new Square("plate1", new double[]{ 1.1, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 10.0, medium, medium, NullInterface.ideal()),
			new Square("plate2", new double[]{ 1.2, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 10.0, medium, medium, NullInterface.ideal()),
			new Square("plate3", new double[]{ 1.3, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 10.0, medium, medium, NullInterface.ideal()),
			new Square("plate4", new double[]{ 1.4, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 10.0, medium, medium, NullInterface.ideal()),
			new Square("plate5", new double[]{ 1.5, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 10.0, medium, medium, NullInterface.ideal()),
			new Square("plate6", new double[]{ 1.6, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 10.0, medium, medium, NullInterface.ideal()),
			new Square("plate7", new double[]{ 1.7, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 10.0, medium, medium, NullInterface.ideal()),
			new Square("plate8", new double[]{ 1.8, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 10.0, medium, medium, NullInterface.ideal()),
			new Square("plateEnd", new double[]{ 1.0 + fw/2, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 10.0, medium, null, IsoUniaxialInterface.ideal()),
			new Square("back", new double[]{ 1.0 + fw/2 + 1.0, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 1, 15.0, null, null, Absorber.ideal())
		};
		
		Optic all = new Optic("all", squares);
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/entryExitPhases.vrml", 0.01);
		vrmlOut.setDrawPolarisationFrames(true);
		
		
		boolean fail = false;
		fail |= doDirectionsCheck(-0.10, iAng, new double[]{ 0,0,1 }, new double[]{ 1,0,1,0 }, all, vrmlOut, medium, wavelen);
		fail |= doDirectionsCheck(-0.06, iAng, new double[]{ 0,1,0 }, new double[]{ 1,0,1,0 }, all, vrmlOut, medium, wavelen);
		fail |= doDirectionsCheck(-0.02, iAng, new double[]{ 0,0,1 }, new double[]{ 1,0,-1,0 }, all, vrmlOut, medium, wavelen);
		fail |= doDirectionsCheck(0.02, iAng, new double[]{ 0,rt2r,rt2r }, new double[]{ 1,0,0,0 }, all, vrmlOut, medium, wavelen);
		fail |= doDirectionsCheck(0.06, iAng, new double[]{ 0,-rt2r,rt2r }, new double[]{ 1,0,0,0 }, all, vrmlOut, medium, wavelen);
		fail |= doDirectionsCheck(0.10, iAng, new double[]{ 0,rt2r,rt2r }, new double[]{ 0,0,1,0 }, all, vrmlOut, medium, wavelen);
		
		vrmlOut.drawOptic(all);
		vrmlOut.destroy();
		
		if(fail)
			fail("Failed, check stderr");
		
	}
	
	private boolean doDirectionsCheck(double y, double incidenceAng, double opticAxis[], double initE[], 
										Optic all, VRMLDrawer vrmlOut, Medium medium, double wavelen){
		
		medium.setOpticAxes(new double[][]{ opticAxis });

		RaySegment ray = new RaySegment();
		ray.dir = Util.reNorm(new double[]{ FastMath.cos(incidenceAng), FastMath.sin(incidenceAng), 0 });		
		ray.startPos = Util.minus(new double[]{ 2, y, 0 }, Util.mul(ray.dir, 12*wavelen));					
		ray.length = Double.POSITIVE_INFINITY;
		ray.up = Util.reNorm(new double[]{ 0, 0, 1 });
		ray.E0 = new double[][]{ initE }; // 45° to the right
		ray.wavelength = wavelen;
		
		Tracer.trace(all, ray, 100, 1E-5, true);
		
		//ray.dumpPath();
		
		//Merge the ordinary and extraordinary rays, to 
		//make the drawing look better and give us
		//the a real sense of 
		ray.rotatePolRefFrame(new double[]{0,0,1});
		RaySegment o = ray.endHit.transmittedOrdinary;
		RaySegment e = ray.endHit.transmittedExtraordinary;
		while(true){
			if(o == null || e == null || o.E1 == null || e.E1 == null){
				System.err.println(y + ": Incomplete ray path");
				return true;
			}
				
			o.rotatePolRefFrame(new double[]{0,0,1});
			e.rotatePolRefFrame(new double[]{0,0,1});
			o.E0[0][0] += e.E0[0][0];
			o.E0[0][1] += e.E0[0][1];
			o.E0[0][2] += e.E0[0][2];
			o.E0[0][3] += e.E0[0][3];
			o.E1[0][0] += e.E1[0][0];
			o.E1[0][1] += e.E1[0][1];
			o.E1[0][2] += e.E1[0][2];
			o.E1[0][3] += e.E1[0][3];
			if(o.endHit == null || o.endHit.transmittedOrdinary == null)
				break;
			o = o.endHit.transmittedOrdinary;
			e = e.endHit.transmittedOrdinary;
		};
		
		double psiInit = Pol.psi(ray.E0[0]);
		double psiEntry = Pol.psi(ray.endHit.transmittedOrdinary.E0[0]);
		double psiExit = Pol.psi(o.startHit.incidentRay.E1[0]);
		double psiEnd = Pol.psi(o.E1[0]);
				
		ray.endHit.transmittedExtraordinary = null;
		System.out.println(y + 
				": psi0 = " + psiInit*180/Math.PI + 
				", psii = " + psiEntry*180/Math.PI + 
				", psix = " + psiExit*180/Math.PI + 
				", psi1 = " + psiEnd*180/Math.PI );
		
		boolean fail = false;

		if(FastMath.abs(psiEntry - psiInit) > 0.5*Math.PI/180){
			System.err.println(" <-- MISMATCH(entry) ");
			fail = true;
		}

		//We're expecting it to be a half wave plate and flip the direction
		if( FastMath.abs(FastMath.abs(psiEnd - psiInit) - 90*Math.PI/180) > 2*Math.PI/180){
			System.err.println(" <-- MISMATCH(exit) ");
			fail = true;
		}
		
		if(FastMath.abs(psiEnd - psiExit) > 0.5*Math.PI/180){
			System.err.println(" <-- MISMATCH(entry) ");
			fail = true;
		}
		
		vrmlOut.drawRay(ray);
		
		Pol.recoverAll();
		
		return fail;
	}
}
