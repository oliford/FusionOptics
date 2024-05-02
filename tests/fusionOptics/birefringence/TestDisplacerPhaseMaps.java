package fusionOptics.birefringence;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.OpticApprox;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoUniaxialInterface;
import fusionOptics.materials.CrystalQuartz;
import fusionOptics.materials.UniaxialFixedIndexGlass;
import fusionOptics.optics.Box;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import net.jafama.FastMath;
import uk.co.oliford.jolu.BinaryMatrixFile;
import junit.framework.TestCase;

/** Check the maps of phase difference match the formula/graph in F.E.Veiras and
 * also that the fringes match reality  (experiment) 
 * 
 * This seems to work, except that the results jumps around by ±180° phase shifts
 * 
 * @author oliford */
public class TestDisplacerPhaseMaps extends TestCase {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/uniAxialTest";
	VRMLDrawer vrmlOut;
	
	public void testVeirasFig9(){
		double ni = 1.0;
		double no = 1.54264;
		double ne = 1.5517;
		double L = 1.973e-3;
		double wavelen = 632.8e-9;
		double theta = 45 * Math.PI/180; //angle of optic axis to surface
		double maxAlpha = 20*Math.PI/180;
		
		int n = 30;
		double fw = Util.calcWaveplateFullWaveThickness(new UniaxialFixedIndexGlass(no, ne), wavelen);
		System.out.println(fw/4);
			
		//double alpha[] = OneLiners.linSpace(0, 10*Math.PI/180, 100);
		//double delta[] = OneLiners.linSpace(0, 360*Math.PI/180, 360);
		
		double phsCalc[][] = new double[n][n];
		double phsTarg[][] = new double[n][n];
		
		boolean fail = false;
		
		for(int iY=0; iY < n; iY++){
			double y = ((double)iY / (0.5*n)) - 1;
			for(int iX=0; iX < n; iX++){
				double x = ((double)iX / (0.5*n)) - 1;
				
				double alpha = FastMath.sqrt(x*x + y*y) * maxAlpha;
				double delta = FastMath.atan2(y,x);
				if(alpha <= maxAlpha){						
					phsCalc[iY][iX] = phaseDiffTracer(ni, no, ne, L, wavelen, alpha, Math.PI/2 - theta, delta);
					phsTarg[iY][iX] = 2*Math.PI*OpticApprox.waveplateOPD(ni, no, ne, theta, delta, alpha, L)/wavelen;
					
					if(phsCalc[iY][iX] != phsTarg[iY][iX])
						fail = true;
				}else{
					phsCalc[iY][iX] = Double.NaN;
					phsTarg[iY][iX] = Double.NaN;
				}
				
				System.out.print(".");
					
			}
			System.out.println();
				
		}
		
		if(vrmlOut != null)
			vrmlOut.destroy();
		
		BinaryMatrixFile.mustWrite(outPath + "/VeirasFig9aCalc.bin", phsCalc, false);
		BinaryMatrixFile.mustWrite(outPath + "/VeirasFig9aTarg.bin", phsTarg, false);
		
		
		if(fail)
			fail("Tracer phase mismatches Veiras calc");
		
	}
	
	private double phaseDiffTracer(double ni, double no, double ne, double thickness, double wavelen, 
										double incidenceAngle, double opticAxisAngleToNormal, 
										double opticAxisAngToIncidencePlane){
		
		
		double normal[] = new double[]{ -1, 0, 0 };
		double incident[] = new double[]{ FastMath.cos(incidenceAngle), FastMath.sin(incidenceAngle), 0 };		
		double opticAxis[] = new double[]{ 
				FastMath.cos(opticAxisAngleToNormal), 
				FastMath.sin(opticAxisAngleToNormal) * FastMath.cos(opticAxisAngToIncidencePlane), 
				FastMath.sin(opticAxisAngleToNormal) * FastMath.sin(opticAxisAngToIncidencePlane) };
		
		Material plateMat = new UniaxialFixedIndexGlass(no, ne);		
		Medium plateMedium = new Medium(plateMat, new double[][]{ opticAxis }, 300);
		
		Box plate = new Box("plate", new double[]{ thickness/2, 0, 0 }, thickness, 0.2, 0.2, plateMedium, IsoUniaxialInterface.ideal());
		double backPos[] = Util.mul(incident, 0.1);
		double backNorm[] = incident.clone();
		Square back = new Square("back", backPos, backNorm, new double[]{ 0, 0, 1 }, 0.2, 0.2, Absorber.ideal());
		Optic all = new Optic("all", new Element[]{ plate, back });
		
		RaySegment ray = new RaySegment();
		
		ray.dir = incident;		
		ray.startPos = new double[]{ -incident[0], -incident[1], -incident[2] };		
		ray.length = Double.POSITIVE_INFINITY;
		ray.up = Util.reNorm(new double[]{ 0, 0, 1 });
		ray.E0 = new double[][]{ { 1,0,1,0 } }; // 45° to the right
		ray.wavelength = wavelen;
		
		Tracer.trace(all, ray, 100, 1e-10, true);

		if(vrmlOut == null){
			vrmlOut = new VRMLDrawer(outPath + "/phaseMaps.vrml", 0.005);
			vrmlOut.setDrawPolarisationFrames(true);
			
			vrmlOut.drawOptic(all);
			
		}
		vrmlOut.drawRay(ray);
		
		RaySegment mergedRay = ray.createMergedRay();
		//VRMLDrawer.dumpRay(outPath + "/displMaps.vrml", all, mergedRay, 0.005);
		
		//check the ray direction
		RaySegment oRay = ray.endHit.transmittedOrdinary;
		RaySegment eRay = ray.endHit.transmittedExtraordinary;
		
		if(oRay == null || eRay == null || oRay.endHit.transmittedOrdinary == null || eRay.endHit.transmittedOrdinary == null){
			System.err.println("Ray path incorrect.");
			return 0;
		}
		
		double u[] = back.getUp();
		double r[] = back.getRight();
		oRay.endHit.transmittedOrdinary.rotatePolRefFrame(u);
		eRay.endHit.transmittedOrdinary.rotatePolRefFrame(u);
		double psiO = Pol.psi(oRay.endHit.transmittedOrdinary.E1[0]);
		double psiE = Pol.psi(eRay.endHit.transmittedOrdinary.E1[0]);
		//System.out.println("Tracer psi: O = " + psiO*180/Math.PI + "\tE = " + psiE*180/Math.PI + "\tDiff = " + (psiE - psiO)*180/Math.PI);
		
		double polODir[] = new double[]{
				FastMath.cos(psiO) * u[0] + Math.sin(psiO) * r[0],
				FastMath.cos(psiO) * u[1] + Math.sin(psiO) * r[1],
				FastMath.cos(psiO) * u[2] + Math.sin(psiO) * r[2],
		};
		oRay.endHit.transmittedOrdinary.rotatePolRefFrame(polODir);
		eRay.endHit.transmittedOrdinary.rotatePolRefFrame(polODir);
		psiO = Pol.psi(oRay.endHit.transmittedOrdinary.E1[0]);
		psiE = Pol.psi(eRay.endHit.transmittedOrdinary.E1[0]);
		double phaseO = FastMath.atan2(oRay.endHit.transmittedOrdinary.E1[0][1], oRay.endHit.transmittedOrdinary.E1[0][0]);
		double phaseE = FastMath.atan2(eRay.endHit.transmittedOrdinary.E1[0][3], eRay.endHit.transmittedOrdinary.E1[0][2]);
		
		long fullWaveDiff = oRay.nWaves + oRay.endHit.transmittedOrdinary.nWaves - eRay.nWaves - eRay.endHit.transmittedOrdinary.nWaves;
		
		double phaseDiffTracer = fullWaveDiff*2*Math.PI + phaseO - phaseE;
		
		Pol.recoverAll();
		
		return phaseDiffTracer;
	}
}
