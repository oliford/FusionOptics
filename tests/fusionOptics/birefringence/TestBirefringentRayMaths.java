package fusionOptics.birefringence;

import java.text.DecimalFormat;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.OpticApprox;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoUniaxialInterface;
import fusionOptics.materials.LithiumNiobate;
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

import otherSupport.ScientificNumberFormat;
import binaryMatrixFile.BinaryMatrixFile;
import algorithmrepository.Algorithms;
import oneLiners.OneLiners;
import net.jafama.FastMath;
import junit.framework.TestCase;

/** Gradually going through the maths of several papers to get a consistent calculation
 * of the ray directions and OPD through a displacer plate or tilted wave plate.
 * 
 *  [ F.E.Veiras "Phase Shift Formulas in Uniaxial Media: An application to waveplates"
 *      J Opt. Soc. Am. 2010  ]
 *  [ M.Avendano "Optical Path Difference in a Plane Parallel uniaxial Plate" 
 *  	J. Opt. Soc. Am. A 2006 ]
 *  [ M.C.Simon "Ray tracing formulas for monoaxial optical components: vectorial formulation"
 *  	1985 ]
 *  [ M. Avendano "Huygen's principle and rays in uniaxial anisotropic media. I Crystal axis normal to refracting surface"
 *  	J. Opt. Soc. Am. A 19 1668 (2002) ]
 *  [ M. Avendano "Huygen's principle and rays in uniaxial anisotropic media. II Crystal axis orientation arbitrary"
 *  	J. Opt. Soc. Am. A 19 1674 (2002) ]
 *  [ M.C.Simon "Waves and rays in uniaxial birefringent crystals" 
 * 		Optik 118 (2007) http://dx.doi.org/10.1016/j.ijleo.2006.03.032 ]
 * 	
 * The first and last papers are the best to use.
 *  
 *  */
public class TestBirefringentRayMaths extends TestCase {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/uniAxialTest";
		
	public static double[][] MAvendano(double ni, double no, double ne, double opticAxis[], double normal[], double incident[]){
				
		// setup our working frame, which has the incidence plane in (x,z) and surface in (x,y) 
		//Z points in the same direction to ray
		double z[] = new double[]{ -normal[0], -normal[1], -normal[2] };
		double y[] = Util.cross(z, incident);
		y = Util.dot(y,y) == 0 ? Util.createPerp(z) : Util.reNorm(y);
		double x[] = Util.cross(y, z);
		
		//Rotate the optic axis to that frame
		double A[] = new double[]{
				Util.dot(opticAxis, x),
				Util.dot(opticAxis, y),
				Util.dot(opticAxis, z),
		};
		
		double N = no*no - ne*ne;

		// Incident angles 
		double cosThetaI = Util.dot(incident, z);
		double sinThetaI = Util.dot(incident, x);
	
		double rIndexRatio = ni / no;
	
		if( sinThetaI*sinThetaI >= 1/(rIndexRatio*rIndexRatio)){
			//total internal refraction
			//Reflector.pureReflection(hit, minIntensity, nonAbsorbedAmplitudeCoeff);
			throw new RuntimeException("Total internal reflection on entry to anisotropic media is not supported.");
		}
		
		//reflection direction
		double reflectedDir[] = new double[]{ sinThetaI, 0, -cosThetaI };
		
		//snells law for ordinary wave
		double sinThetaO = sinThetaI * ni / no;
		double ordinaryDir[] = (new double[]{
				sinThetaO, 
				0, 
				Math.sqrt(1 - sinThetaO*sinThetaO), 
			});
		System.out.println("MA: |No| = " + Util.length(ordinaryDir)); //should be normalised
		
		
		//extra ordinary dir and fresnel coeffs from formulae in papers
		double ret[][] = IsoUniaxialInterface.calcTransmittedExtraordinaryDir(ni, no, ne, ordinaryDir, A);
		double extraordinaryDir[] = ret[0];
		double exWaveNormal[] = ret[1];
		
		double cosThetaN = Util.dot(A, extraordinaryDir);
		double cosThetaA = Util.dot(A, exWaveNormal);
		double sinSqThetaN = 1-cosThetaN*cosThetaN;
		double sinSqThetaA = 1-cosThetaA*cosThetaA;
		
		double nNE = no*ne / Math.sqrt(no*no*sinSqThetaA + ne*ne*cosThetaA*cosThetaA); //equ 31a
		double nRE = Math.sqrt(ne*ne*sinSqThetaN + no*no*cosThetaN*cosThetaN); //equ 31b
		
		
		return new double[][]{ 
				IsoUniaxialInterface.reFrame(x,y,z, ordinaryDir), 
				IsoUniaxialInterface.reFrame(x,y,z, exWaveNormal), 
				IsoUniaxialInterface.reFrame(x,y,z, extraordinaryDir), 
				new double[]{ nNE, nRE} };
	}
	
	
	public static double[][] SimonEcharri(double ni, double no, double ne, double opticAxis[], double normal[], double incident[]){
		
		double S[] = incident;
		double z3[] = opticAxis;
		double n[] = new double[]{ -normal[0], -normal[1], -normal[2] };
		double x[] = n;
		double y[] = Util.reNorm(Util.cross(x, z3));
		double z[] = Util.cross(x, y);
		
		double u = 1.0 / ni; // speed of light = 1 here, we don't need it
		double uo = 1.0 / no;
		double ue = 1.0 / ne;
		
		double sn = Util.dot(S, n);
		double sn2 = sn*sn;
				
		double ao = FastMath.sqrt((u/uo)*(u/uo) - 1 + sn2) - sn;
		double rnNo = FastMath.sqrt(1 + ao*ao + 2*ao*sn);
		double No[] = new double[]{
				(S[0] + ao * n[0]) / rnNo,
				(S[1] + ao * n[1]) / rnNo,
				(S[2] + ao * n[2]) / rnNo,
		};
		System.out.println("SE: |No| = " + Util.length(No)); //should be normalised
		
		
		double b = (uo*uo - ue*ue) / (u*u);
		double nXz3[] = Util.cross(n, z3);
		double snz32 = FastMath.pow2(Util.dot(S, nXz3));
		double z3n2 = FastMath.pow2(Util.dot(z3, n));
		double ueu2 = (ue/u)*(ue/u);
		
		double A = FastMath.pow2(1 + b*(1 - sn2 - snz32)) 
						-4*b*((1-sn2)*(1-z3n2) - snz32);
		
		double B = 2*(1 + b*(1 - sn2 - snz32)) * (b*z3n2 + ueu2) 
					- 4*b*ueu2*((1-sn2)*(1-z3n2) - snz32);
		
		double C = FastMath.pow2(b*z3n2 + ueu2);
		
		double B24AC = B*B - 4*A*C;
		B24AC = (B24AC < 0) ? 0 : FastMath.sqrt(B24AC); 
				
		double wp = FastMath.sqrt((B + B24AC)/(2*A));
		double wm = FastMath.sqrt((B - B24AC)/(2*A));
		
		double sel = (uo-ue)*Util.dot(z3,z)*Util.dot(z3,n)*Util.dot(S,z);
		double w = (sel > 0) ? wp : wm;
		
		double ae = FastMath.sqrt((1/w)*(1/w) - 1 + sn2) - sn;
		double rnNe = FastMath.sqrt(1 + ae*ae + 2*ae*sn);
		double Ne[] = new double[]{
				(S[0] + ae*n[0]) / rnNe,
				(S[1] + ae*n[1]) / rnNe,
				(S[2] + ae*n[2]) / rnNe,
		};
		System.out.println("SE: |Ne| = " + Util.length(Ne)); //should be normalised
		
		double Nez32 = FastMath.pow2(Util.dot(Ne, z3));
		double ux = FastMath.sqrt((1 - Nez32)*ue*ue + Nez32 * uo*uo);

		double N3 = Util.dot(Ne,z3);
		double rnRe = FastMath.sqrt((1-N3*N3)*ue*ue*ue*ue + N3*N3*uo*uo*uo*uo);
		double Re[] = new double[]{
				(Ne[0]*ue*ue + (uo*uo - ue*ue)*N3*z3[0]) / rnRe,
				(Ne[1]*ue*ue + (uo*uo - ue*ue)*N3*z3[1]) / rnRe,
				(Ne[2]*ue*ue + (uo*uo - ue*ue)*N3*z3[2]) / rnRe,				
			};
		System.out.println("SE: |Re| = " + Util.length(Re)); //should be normalised
		
		System.out.println("SE: No.y:\t" + Util.dot(S, y)/u + "\t" + Util.dot(No,y)/uo);
		System.out.println("SE: No.z:\t" + Util.dot(S, z)/u + "\t" + Util.dot(No,z)/uo);

		System.out.println("SE: Ne.y:\t" + Util.dot(S, y)/u + "\t" + Util.dot(Ne,y)/ux);
		System.out.println("SE: Ne.z:\t" + Util.dot(S, z)/u + "\t" + Util.dot(Ne,z)/ux);
		
		System.out.println("SE: wp = " + wp + "\twm = "  + wm + "\tchosen=" + w + "\tux/u=" + ux/u);
				
		return new double[][]{ No, Ne, Re, new double[]{ 1.0/ux } };
	}
	
	public void testAll(){
		double incidenceAngle = 5 * Math.PI/180;
		double opticAxisAngleToNormal = 66 * Math.PI/180;
		double opticAxisAngToIncidencePlane = 78 * Math.PI/180;
		double ni = 1.0;
		//double no = 1.1;
		//double ne = 1.2;
		double no = 1.54264;
		double ne = 1.5517;
		double thickness = 1.973e-3;
		double wavelen = 632.8e-9; 
		
		//singleRayTest(wavelen, ni, no, ne, thickness, 0.16748108983784, opticAxisAngleToNormal, opticAxisAngToIncidencePlane);
		//if(true)return;
		
		double ang[] = OneLiners.linSpace(-10*Math.PI/180, 10*Math.PI/180, 100);
		double ret[][] = new double[ang.length][];
		
		for(int i=0; i < ang.length; i++){
			incidenceAngle = ang[i];
			
			ret[i] = singleRayTest(wavelen, ni, no, ne, thickness, incidenceAngle, opticAxisAngleToNormal, opticAxisAngToIncidencePlane);
			
		}
		
		BinaryMatrixFile.mustWrite(outPath + "/mathsPhaseDiffTracer.bin", ret, false);
		
		/*double ni = 1.0;
		double no = 1.54264;
		double ne = 1.5517;
		double L = 1.973e-3;
		double wavelen = 632.8e-9;
		*/
	}
	
	private double[] singleRayTest(double wavelen, double ni, double no, double ne, double thickness, double incidenceAngle, 
								double  opticAxisAngleToNormal, double opticAxisAngToIncidencePlane) {
		
		
		double normal[] = new double[]{ -1, 0, 0 };
		double incident[] = new double[]{ FastMath.cos(incidenceAngle), FastMath.sin(incidenceAngle), 0 };		
		double opticAxis[] = new double[]{ 
				FastMath.cos(opticAxisAngleToNormal), 
				FastMath.sin(opticAxisAngleToNormal) * FastMath.cos(opticAxisAngToIncidencePlane), 
				FastMath.sin(opticAxisAngleToNormal) * FastMath.sin(opticAxisAngToIncidencePlane) };
		
		double ret[][] = SimonEcharri(ni, no, ne, opticAxis, normal, incident);
		double No[] = ret[0], Ne[] = ret[1], Re[] = ret[2], nNe = ret[3][0];
		
		ret = MAvendano(ni, no, ne, opticAxis, normal, incident);
		double NoAv[] = ret[0], NeAv[] = ret[1], ReAv[] = ret[2], nNeAv = ret[3][0], nReAv = ret[3][1];

		
		double thetaOSnells = Math.asin(Math.sin(incidenceAngle)*1.0/no);
		double thetaOSimon = FastMath.asin(No[1]);
		double thetaOAvendano = FastMath.asin(NoAv[1]);
		
		double thetaNESnells = Math.asin(Math.sin(incidenceAngle)*1.0/ne);
		double thetaNESimon = FastMath.asin(Ne[1]);
		double thetaNEAvendano = FastMath.asin(NeAv[1]);
		
		double thetaRESimon = FastMath.asin(Re[1]);
		double thetaREAvendano = FastMath.asin(ReAv[1]);
		double thetaOutSimon = FastMath.asin(Re[2]);
		double thetaOutAvendano = FastMath.asin(Re[2]);
		
		
		//*/
		
		
		/* ******* Geometric OPD ********* */

		double exitO[] = new double[]{
				No[0] * (thickness / No[0]),
				No[1] * (thickness / No[0]),
				No[2] * (thickness / No[0]),
		};
		
		double exitE[] = new double[]{
				Re[0] * (thickness / Re[0]),
				Re[1] * (thickness / Re[0]),
				Re[2] * (thickness / Re[0]),
		};
		
		double lenO = Util.length(exitO);
		double oplO = lenO * no;
		
		double nRe = nNe * Util.dot(Re, Ne);
		double lenE = Util.length(exitE);
		double oplE = lenE * nRe;		
		
		//the distance along exited E ray from it's exit point before
		//reaching the plane perp to it that contains the O exit point
		double distEtoOplane = Util.dot(Util.minus(exitO, exitE), incident);
		double oplE2 = distEtoOplane * 1.0; //in vacuum
		
		double lto = exitO[1];
		double lte = exitE[1];
		double oplE2b = -ni * (lte - lto) * FastMath.sin(incidenceAngle);
		
		double lSigmae = exitE[2];
		
		double opdGeom = (oplO - oplE - oplE2);
		double phaseShiftGeomOPD = 2 * Math.PI * opdGeom / wavelen;
		
				
		/* ******* Veiras OPD ********* */		
		double theta = Math.PI/2 - opticAxisAngleToNormal;
		double sinSqAlpha = FastMath.pow2(FastMath.sin(incidenceAngle));
		double sinSqTheta = FastMath.pow2(FastMath.sin(theta));
		double cosSqTheta = FastMath.pow2(FastMath.cos(theta));
		double sinSqDelta = FastMath.pow2(FastMath.sin(opticAxisAngToIncidencePlane));
		double oplOVeiras = thickness * no*no/FastMath.sqrt(no*no - ni*ni*sinSqAlpha); //Verias equ7
		
		//Verias equ9
		double oplEVeiras = thickness*no*ne*ne/FastMath.sqrt(
					ne*ne*(ne*ne*sinSqTheta + no*no*cosSqTheta) 
					- (ne*ne - (ne*ne - no*no)*cosSqTheta*sinSqDelta)*ni*ni*sinSqAlpha );
		
		double oplEVeiras2 = nNe * thickness * Util.dot(Re, Ne) / Re[0];
		
		double opdVeiras = OpticApprox.waveplateOPD(ni, no, ne, theta, opticAxisAngToIncidencePlane, incidenceAngle, thickness);
		
		/* ******* MCSimon Waves abd Rays in Uniaxial ..., derivation of Veiras ********* */
		
		double hx = no*no + (ne*ne - no*no)*opticAxis[0]*opticAxis[0];
		double ht = no*no + (ne*ne - no*no)*opticAxis[1]*opticAxis[1];
		double hSigma = ne*ne - (ne*ne - no*no)*opticAxis[2]*opticAxis[2];
		double hxt = (ne*ne - no*no)*opticAxis[0]*opticAxis[1];
		double hxSigma = (ne*ne - no*no)*opticAxis[0]*opticAxis[2];
		double hSigmat = (ne*ne - no*no)*opticAxis[2]*opticAxis[1];
		
		double Phi = ne*ne*hx - hSigma*ni*ni*incident[1]*incident[1];
		
		double lte2 = thickness * (hxt/hx + no*hSigma*ni*incident[1]/hx/FastMath.sqrt(Phi));
		double lSigmae2 = thickness * (hxSigma/hx + no*hSigmat*ni*incident[1]/hx/FastMath.sqrt(Phi));

		//double oplE2b = ni * (lte - lto) * FastMath.sin(incidenceAngle);
		double opdE2Simon = -ni*(lte2 - lto)*FastMath.sin(incidenceAngle);
				
		/* ******* Avendano OPD (delta=0 only) ********* */
		//double thetaN = (opticAxisAngToIncidencePlane != 0) ? Double.NaN : FastMath.acos(Util.dot(opticAxis, Re));
		double thetaN = FastMath.acos(Util.dot(opticAxis, Re));
		double thetaA = FastMath.acos(Util.dot(opticAxis, Ne));
		double sinThetaA = FastMath.sin(thetaA);
		double cosThetaA = FastMath.cos(thetaA);
		double sinThetaN = FastMath.sin(thetaN);
		double cosThetaN = FastMath.cos(thetaN);
		
		double nwf = no*ne / Math.sqrt(no*no*sinThetaA*sinThetaA + ne*ne*cosThetaA*cosThetaA); //equ 31a
		double npv = Math.sqrt(ne*ne*sinThetaN*sinThetaN + no*no*cosThetaN*cosThetaN); //equ 31b
		
		
		double cosThetaO = No[0] / Util.length(No);
		double cosThetaE = Re[0] / Util.length(Re);
		double tanThetaO = FastMath.tan(FastMath.acos(cosThetaO));
		double tanThetaE = FastMath.tan(FastMath.acos(cosThetaE));
		double sinThetaI = FastMath.sin(incidenceAngle);
		
		double oplEAvendano1 = npv * thickness / cosThetaE; //equ 32
		
		double opdAvendano = thickness * ((no/cosThetaO) - (npv / cosThetaE) + ni*(tanThetaE - tanThetaO)*sinThetaI); //equ 33
		
		
		/* ************* Ray Tracer ************ */
		
		//test the tracer with the optic axis coordinate aligned and the incident ray at an angle
		incident = new double[]{ FastMath.cos(incidenceAngle), 
				FastMath.sin(incidenceAngle) * FastMath.cos(opticAxisAngToIncidencePlane), 
				FastMath.sin(incidenceAngle) * FastMath.sin(opticAxisAngToIncidencePlane) };		
		opticAxis = new double[]{ 
				FastMath.cos(opticAxisAngleToNormal), 
				FastMath.sin(opticAxisAngleToNormal) , 
				0 };
		
		
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
		ray.up = Util.reNorm(Util.cross(Util.cross(incident, new double[]{ 0, 0, 1 }), incident));
		ray.E0 = new double[][]{ { 1,0,1,0 } }; // 45° to the right
		ray.wavelength = wavelen;
		
		Tracer.trace(all, ray, 100, 0.01, true);

		RaySegment mergedRay = ray.createMergedRay();
		VRMLDrawer.dumpRay(outPath + "/mathsTest.vrml", all, mergedRay, 0.005);
		
		//check the ray direction
		RaySegment oRay = ray.endHit.transmittedOrdinary;
		RaySegment eRay = ray.endHit.transmittedExtraordinary;
		
		double thetaOTracer = FastMath.asin(oRay.dir[1]);
		double thetaETracer = FastMath.asin(eRay.dir[1]);
		double thetaOutTracer = FastMath.asin(eRay.dir[2]);
		
		double ltoTracer = oRay.endHit.pos[1];
		double lteTracer = eRay.endHit.pos[1];
		double lSigmaeTracer = eRay.endHit.pos[2];
		
		double oplOTracer = oRay.length * no;
		double oplETracer = eRay.length * eRay.raySpecificRefractiveIndex;
	
		double opdTracer = oplOTracer + oRay.endHit.transmittedOrdinary.length - oplETracer - eRay.endHit.transmittedOrdinary.length;	
		
		double opdTracer2 = 0;
		double u[] = back.getUp();
		double r[] = back.getRight();
		oRay.endHit.transmittedOrdinary.rotatePolRefFrame(u);
		eRay.endHit.transmittedOrdinary.rotatePolRefFrame(u);
		double psiO = Pol.psi(oRay.endHit.transmittedOrdinary.E1[0]);
		double psiE = Pol.psi(eRay.endHit.transmittedOrdinary.E1[0]);
		System.out.println("Tracer psi: O = " + psiO*180/Math.PI + "\tE = " + psiE*180/Math.PI + "\tDiff = " + (psiE - psiO)*180/Math.PI);
		
		double polODir[] = new double[]{
				FastMath.cos(psiO) * u[0] + Math.sin(psiO) * r[0],
				FastMath.cos(psiO) * u[1] + Math.sin(psiO) * r[1],
				FastMath.cos(psiO) * u[2] + Math.sin(psiO) * r[2],
		};
		ray.rotatePolRefFrame(new double[]{ 0, 0, 1 });
		oRay.rotatePolRefFrame(new double[]{ 0, 0, 1 });
		eRay.rotatePolRefFrame(new double[]{ 0, 0, 1 });
		oRay.endHit.transmittedOrdinary.rotatePolRefFrame(new double[]{ 0, 0, 1 });
		eRay.endHit.transmittedOrdinary.rotatePolRefFrame(new double[]{ 0, 0, 1 });
		
		psiO = Pol.psi(oRay.endHit.transmittedOrdinary.E1[0]);
		psiE = Pol.psi(eRay.endHit.transmittedOrdinary.E1[0]);
		System.out.println("Tracer psi(O as up): O = " + psiO*180/Math.PI + "\tE = " + psiE*180/Math.PI + "\tDiff = " + (psiE - psiO)*180/Math.PI);
				
		//simply (geom) calculated phase rotations of each section (using tracers indecies and ray lengths)
		double fullPhaseSeg0Geom = 2 * Math.PI * ray.length / wavelen;
		int nWavesSeg0Geom = (int)(fullPhaseSeg0Geom / (2*Math.PI));
		double partPhaseSeg0Geom = fullPhaseSeg0Geom - nWavesSeg0Geom * 2 * Math.PI;
		
		double fullPhaseSeg1GeomO = fullPhaseSeg0Geom + 2 * Math.PI * oRay.length * no / wavelen;
		int nWavesSeg1GeomO = (int)(fullPhaseSeg1GeomO / (2*Math.PI));
		double partPhaseSeg1GeomO = fullPhaseSeg1GeomO - nWavesSeg1GeomO * 2 * Math.PI;
		
		double fullPhaseSeg1GeomE =  fullPhaseSeg0Geom + 2 * Math.PI * eRay.length * eRay.raySpecificRefractiveIndex / wavelen;
		int nWavesSeg1GeomE = (int)(fullPhaseSeg1GeomE / (2*Math.PI));
		double partPhaseSeg1GeomE = fullPhaseSeg1GeomE - nWavesSeg1GeomE * 2 * Math.PI;
		
		double fullPhaseSeg2GeomO = fullPhaseSeg1GeomO + 2 * Math.PI * oRay.endHit.transmittedOrdinary.length / wavelen;
		int nWavesSeg2GeomO = (int)(fullPhaseSeg2GeomO / (2*Math.PI));
		double partPhaseSeg2GeomO = fullPhaseSeg2GeomO - nWavesSeg2GeomO * 2 * Math.PI;
		
		double fullPhaseSeg2GeomE =  fullPhaseSeg1GeomE + 2 * Math.PI * eRay.endHit.transmittedOrdinary.length / wavelen;
		int nWavesSeg2GeomE = (int)(fullPhaseSeg2GeomE / (2*Math.PI));
		double partPhaseSeg2GeomE = fullPhaseSeg2GeomE - nWavesSeg2GeomE * 2 * Math.PI;
		
		double phaseDiffGeom = fullPhaseSeg2GeomO - fullPhaseSeg2GeomE;
		
		//phase and wave count for each ray as calced by Tracer
		int nWavesSeg0Tracer = ray.nWaves;
		double partPhaseSeg0TracerO = FastMath.atan2(ray.E1[0][1], ray.E1[0][0]); 
		double partPhaseSeg0TracerE = FastMath.atan2(ray.E1[0][3], ray.E1[0][2]); 
		if(partPhaseSeg0TracerO < 0)partPhaseSeg0TracerO += 2*Math.PI; 
		if(partPhaseSeg0TracerE < 0)partPhaseSeg0TracerE += 2*Math.PI;
		double fullPhaseSeg0TracerO =  nWavesSeg0Tracer*2*Math.PI + partPhaseSeg0TracerO;
		double fullPhaseSeg0TracerE =  nWavesSeg0Tracer*2*Math.PI + partPhaseSeg0TracerE;
				
		int nWavesSeg1sTracer = ray.nWaves; //start of seg in crystal (to catch interface phase changes)
		double partPhaseSeg1sTracerO = FastMath.atan2(oRay.E0[0][1], oRay.E0[0][0]);  
		double partPhaseSeg1sTracerE = FastMath.atan2(eRay.E0[0][3], eRay.E0[0][2]);  
		if(partPhaseSeg1sTracerO < 0)partPhaseSeg1sTracerO += 2*Math.PI; 
		if(partPhaseSeg1sTracerE < 0)partPhaseSeg1sTracerE += 2*Math.PI; 
		double fullPhaseSeg1sTracerO =  nWavesSeg1sTracer*2*Math.PI + partPhaseSeg1sTracerO;
		double fullPhaseSeg1sTracerE =  nWavesSeg1sTracer*2*Math.PI + partPhaseSeg1sTracerE;
		
		int nWavesSeg1TracerO = ray.nWaves + oRay.nWaves; //rays in crystal, phase change during propagation
		int nWavesSeg1TracerE = ray.nWaves + eRay.nWaves;
		double partPhaseSeg1TracerO = FastMath.atan2(oRay.E1[0][1], oRay.E1[0][0]);  
		double partPhaseSeg1TracerE = FastMath.atan2(eRay.E1[0][3], eRay.E1[0][2]);  
		if(partPhaseSeg1TracerO < 0)partPhaseSeg1TracerO += 2*Math.PI; 
		if(partPhaseSeg1TracerE < 0)partPhaseSeg1TracerE += 2*Math.PI; 
		double fullPhaseSeg1TracerO =  nWavesSeg1TracerO*2*Math.PI + partPhaseSeg1TracerO;
		double fullPhaseSeg1TracerE =  nWavesSeg1TracerE*2*Math.PI + partPhaseSeg1TracerE;
		
		int nWavesSeg2TracerO = ray.nWaves + oRay.nWaves + oRay.endHit.transmittedOrdinary.nWaves; //ray after leaving crystal
		int nWavesSeg2TracerE = ray.nWaves + eRay.nWaves + eRay.endHit.transmittedOrdinary.nWaves;
		double partPhaseSeg2TracerO = FastMath.atan2(oRay.endHit.transmittedOrdinary.E1[0][1], oRay.endHit.transmittedOrdinary.E1[0][0]);  
		double partPhaseSeg2TracerE = FastMath.atan2(eRay.endHit.transmittedOrdinary.E1[0][3], eRay.endHit.transmittedOrdinary.E1[0][2]);  
		if(partPhaseSeg2TracerO < 0)partPhaseSeg2TracerO += 2*Math.PI; 
		if(partPhaseSeg2TracerE < 0)partPhaseSeg2TracerE += 2*Math.PI; 
		double fullPhaseSeg2TracerO =  nWavesSeg2TracerO*2*Math.PI + partPhaseSeg2TracerO;
		double fullPhaseSeg2TracerE =  nWavesSeg2TracerE*2*Math.PI + partPhaseSeg2TracerE;
		
		double phaseDiffTracer = fullPhaseSeg2TracerO - fullPhaseSeg2TracerE;
		
		System.out.println("Seg0: Geom=["+nWavesSeg0Geom+"w + "+partPhaseSeg0Geom*180/Math.PI+"° = "+fullPhaseSeg0Geom*180/Math.PI+"°]\n" +
				 "\t TracerO=["+nWavesSeg0Tracer+"w + "+partPhaseSeg0TracerO*180/Math.PI+"° = "+fullPhaseSeg0TracerO*180/Math.PI+"°] Δ="+(fullPhaseSeg0TracerO-fullPhaseSeg0Geom)*180/Math.PI+"\n" + 
					"\t TracerE=["+nWavesSeg0Tracer+"w + "+partPhaseSeg0TracerE*180/Math.PI+"° = "+fullPhaseSeg0TracerE*180/Math.PI+"°] Δ="+(fullPhaseSeg0TracerE-fullPhaseSeg0Geom)*180/Math.PI+"\n" +
				 "\t TracerO1s=["+nWavesSeg1sTracer+"w + "+partPhaseSeg1sTracerO*180/Math.PI+"° = "+fullPhaseSeg1sTracerO*180/Math.PI+"°] Δ="+(fullPhaseSeg1sTracerO-fullPhaseSeg0Geom)*180/Math.PI+"\n" + 
				"\t TracerE1s=["+nWavesSeg1sTracer+"w + "+partPhaseSeg1sTracerE*180/Math.PI+"° = "+fullPhaseSeg1sTracerE*180/Math.PI+"°] Δ="+(fullPhaseSeg1sTracerE-fullPhaseSeg0Geom)*180/Math.PI);

		System.out.println("Seg1: GeomO=["+nWavesSeg1GeomO+"w + "+partPhaseSeg1GeomO*180/Math.PI+"° = "+fullPhaseSeg1GeomO*180/Math.PI+"°]\n" +
				"\t GeomE=["+nWavesSeg1GeomE+"w + "+partPhaseSeg1GeomE*180/Math.PI+"° = "+fullPhaseSeg1GeomE*180/Math.PI+"°]\n" +
			    "\t TracerO=["+nWavesSeg1TracerO+"w + "+partPhaseSeg1TracerO*180/Math.PI+"° = "+fullPhaseSeg1TracerO*180/Math.PI+"°] Δ="+(fullPhaseSeg1TracerO-fullPhaseSeg1GeomO)*180/Math.PI+"\n" + 
				"\t TracerE=["+nWavesSeg1TracerE+"w + "+partPhaseSeg1TracerE*180/Math.PI+"° = "+fullPhaseSeg1TracerE*180/Math.PI+"°] Δ="+(fullPhaseSeg1TracerE-fullPhaseSeg1GeomE)*180/Math.PI);

		System.out.println("Seg1: GeomO=["+nWavesSeg2GeomO+"w + "+partPhaseSeg2GeomO*180/Math.PI+"° = "+fullPhaseSeg2GeomO*180/Math.PI+"°]\n" +
				"\t GeomE=["+nWavesSeg2GeomE+"w + "+partPhaseSeg2GeomE*180/Math.PI+"° = "+fullPhaseSeg2GeomE*180/Math.PI+"°]\n" +
			    "\t TracerO=["+nWavesSeg2TracerO+"w + "+partPhaseSeg2TracerO*180/Math.PI+"° = "+fullPhaseSeg2TracerO*180/Math.PI+"°] Δ="+(fullPhaseSeg2TracerO-fullPhaseSeg2GeomO)*180/Math.PI+"\n" + 
				"\t TracerE=["+nWavesSeg2TracerE+"w + "+partPhaseSeg2TracerE*180/Math.PI+"° = "+fullPhaseSeg2TracerE*180/Math.PI+"°] Δ="+(fullPhaseSeg2TracerE-fullPhaseSeg2GeomE)*180/Math.PI);

		/* ************* Output *************** */

		System.out.println();
		System.out.println("thetaO: Snells = " + thetaOSnells*180/Math.PI + "°\tSimon = " + thetaOSimon*180/FastMath.PI + "°\tAvendano = " + thetaOAvendano*180/FastMath.PI + "\tTracer = " + thetaOTracer*180/Math.PI);
		System.out.println("thetaNE: Snells = [" + thetaNESnells*180/Math.PI + "°]\tSimon = " + thetaNESimon*180/FastMath.PI + "°\tAvendano = " + thetaNEAvendano*180/FastMath.PI);
		System.out.println("thetaRE: Simon = " + thetaRESimon*180/FastMath.PI + "°\tAvendano = " + thetaREAvendano*180/FastMath.PI + "\tTracer = " + thetaETracer*180/Math.PI);
		System.out.println("thetaREOut: Simon = " + thetaOutSimon*180/FastMath.PI + "°\tAvendano = " + thetaOutAvendano*180/FastMath.PI + "\tTracer = " + thetaOutTracer*180/Math.PI);
		
		System.out.println("nwf: Avendano(ret) = " + nNeAv + "\tAvendano =" + nwf + "\tnx = " + nNe);
		System.out.println("npv: Avendano(ret) = " + nReAv + "\tAvendano =" + npv + "\tnRE = " + nRe + "\tTracer: " + eRay.raySpecificRefractiveIndex);

		System.out.println("lto: geom =" + lto + "\tTracer = " + ltoTracer);
		System.out.println("lte: geom =" + lte + "\tSimon = " + lte2 + "\tTracer = " + lteTracer);
		System.out.println("lσe: geom =" + lSigmae + "\tSimon = " + lSigmae2 + "\tTracer = " + lSigmaeTracer);
		
		System.out.println("oplO: geom = " + oplO + "\tVeiras = " + oplOVeiras + "\tTracer = " + oplOTracer);
		System.out.println("oplE: geom = " + oplE + "\tVeiras = " + oplEVeiras + "\tVeiras2 = " + oplEVeiras2 + "\tAvendano1 = " + oplEAvendano1+ "\tTracer = " + oplETracer);
		System.out.println("oplE2: geom = " + oplE2 + "\tpartVeiras = " + oplE2b + "\tVeiras(opd - oplO - oplE) = " + (oplO - oplE - opdVeiras) + "\tSimonVerias = " + opdE2Simon);
		System.out.println("opd: geom = " + opdGeom + "\tVerias = " + opdVeiras + "\tAvendano = [" + opdAvendano + "]\tTracer = " + opdTracer + ", " + opdTracer2);
		
		
		/*System.out.println("Tracer waves: Diff = " + fullWaveDiff);
		System.out.println("Tracer psi: O = " + psiO*180/Math.PI + "\tE = " + psiE*180/Math.PI + "\tDiff = " + (psiO - psiE)*180/Math.PI);
		System.out.println("Tracer phs: O = " + phaseO*180/Math.PI + "\tE = " + phaseE*180/Math.PI + "\tDiff = " + phaseDiffTracer*180/Math.PI);
		*/
		System.out.println("phaseShift: Veiras = " + phaseShiftGeomOPD*180/Math.PI + "\tGeom = " + phaseDiffGeom*180/Math.PI 
								+ "\tTracer = " + phaseDiffTracer*180/Math.PI + "\t(Tracer-Veiras)=" + (phaseDiffTracer - phaseDiffGeom)*180/Math.PI);
		System.out.println("phaseShift%360: Veiras = " + (phaseShiftGeomOPD*180/Math.PI)%360 + "\tGeom = " + (phaseDiffGeom*180/Math.PI)%360
								+ "\tTracer = " + (phaseDiffTracer*180/Math.PI)%360);
		
		double allowedDelta = 1e-6;
		
		assertEquals(thetaOSnells, thetaOSimon, allowedDelta);
		assertEquals(thetaOSnells, thetaOAvendano, allowedDelta);
		assertEquals(thetaOSnells, thetaOTracer, allowedDelta);
		
		assertEquals(thetaNESimon, thetaNEAvendano, allowedDelta);
		assertEquals(thetaRESimon, thetaREAvendano, allowedDelta);
		assertEquals(thetaOutSimon, thetaOutAvendano, allowedDelta);
		assertEquals(thetaRESimon, thetaETracer, allowedDelta);
		assertEquals(thetaOutSimon, thetaOutTracer, allowedDelta);
		
		assertEquals(nwf, nNe, allowedDelta);
		assertEquals(npv, nRe, allowedDelta);
		assertEquals(npv, nReAv, allowedDelta);
		assertEquals(npv, eRay.raySpecificRefractiveIndex, allowedDelta);

		assertEquals(lte, lte2, allowedDelta);
		assertEquals(lSigmae, lSigmae2, allowedDelta);
		assertEquals(lto, ltoTracer, allowedDelta);
		assertEquals(lte, lteTracer, allowedDelta);
		assertEquals(lSigmae, lSigmaeTracer, allowedDelta);
		
		assertEquals(oplO, oplOVeiras, allowedDelta);
		assertEquals(oplO, oplOTracer, allowedDelta);
		assertEquals(oplE, oplEVeiras, allowedDelta);
		assertEquals(oplE, oplEVeiras2, allowedDelta);
		assertEquals(oplE, oplEAvendano1, allowedDelta);
		assertEquals(oplE, oplETracer, allowedDelta);
		
		assertEquals(oplE2, oplE2b, allowedDelta);
		assertEquals(oplE2, opdE2Simon, allowedDelta);
		
		assertEquals(opdGeom, opdVeiras, allowedDelta);
		assertEquals(opdGeom, opdTracer, allowedDelta);

		assertEquals(phaseShiftGeomOPD, phaseDiffGeom, allowedDelta);
		assertEquals(phaseShiftGeomOPD%(2*Math.PI), phaseDiffTracer%(2*Math.PI), allowedDelta);
		
		//we don't expect this to work if optic axis isn't in incidence plane
		//since the last Avendano paper isn't general enough
		if(opticAxisAngToIncidencePlane == 0)
			assertEquals(opdGeom, opdAvendano, allowedDelta);

		return new double[]{ incidenceAngle, phaseDiffGeom, phaseShiftGeomOPD, phaseDiffTracer };
		
		
	}
	
}
