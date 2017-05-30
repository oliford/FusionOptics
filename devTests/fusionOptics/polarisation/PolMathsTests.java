package fusionOptics.polarisation;

import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import binaryMatrixFile.BinaryMatrixFile;

/** Low level test/dev platform for polarisation maths in Pol class */
public class PolMathsTests {

	// Some basic tests for polarisation
	public static void main(String[] args) {
		//testPhaseTangents();
		testPols();
		//modulusTest();
	
	}
	
	// Even a signed int (31-bit) can cope with wave counting up to 1.3km 
	public static void modulusTest() {
		double approxDist = 30; //30m = wave long 
		double wavelength = 653e-9; //653nm
		double waveFrac = 0.01;
		
		int nWaveLengths = (int)(approxDist / wavelength);
		double fullWaveDist = nWaveLengths * wavelength;
		
		fullWaveDist += wavelength * waveFrac;
		
		int nWaves = (int)(fullWaveDist / wavelength);
			
		double partWave = (fullWaveDist % wavelength) / wavelength;
		double partWave2 = (fullWaveDist - (nWaves * wavelength)) / wavelength;
				
		System.out.println("approxDist = " + approxDist);
		System.out.println("wavelength = " + wavelength);
		System.out.println("nWaveLengths = " + nWaveLengths);
		System.out.println("fullWaveDist = " + fullWaveDist);
		System.out.println("nWaves = " + nWaves);
		System.out.println("partWave = " + partWave +"\t" + (partWave - waveFrac));
		System.out.println("partWave2 = " + partWave2 +"\t" + (partWave2 - waveFrac));
	}
	
	public static void testPols() {
	
		double E[] = new double[]{ 3,0,3,0 };
		
		System.out.println(
		"\nintensity = " + Pol.intensity(E) + 
		"\nmag2Eu = " + Pol.mag2Eu(E) + 
		"\nmag2Er = " + Pol.mag2Er(E) + 
		"\ntanPhaseEu = " + Pol.tanPhaseEu(E) + 
		"\ntanPhaseEr = " + Pol.tanPhaseEr(E) + 
		"\ncosPhaseEu = " + Pol.cosPhaseEu(E) + 
		"\nsinPhaseEu = " + Pol.sinPhaseEu(E) + 
		"\ncosPhaseEr = " + Pol.cosPhaseEr(E) + 
		"\nsinPhaseEr = " + Pol.sinPhaseEr(E) + 
		"\ntanPhaseDiff = " + Pol.tanPhaseDiff(E) + 
		"\ncosPhaseDiff = " + Pol.cosPhaseDiff(E) + 
		"\nsinPhaseDiff = " + Pol.sinPhaseDiff(E) + 
		"\npolarisedIntensityFrac = " + Pol.polarisedIntensityFrac(E) +
		"\ns0 = " + Pol.s0(E) +
		"\ns1 = " + Pol.s1(E) +
		"\ns2 = " + Pol.s2(E) +
		"\ns3 = " + Pol.s3(E) + 
		"\nsqrt(s1^2+s2^2+s3^2) / s0 = " + Math.sqrt(
					Pol.s1(E)*Pol.s1(E) + 
					Pol.s2(E)*Pol.s2(E) + 
					Pol.s3(E)*Pol.s3(E)
				) / Pol.s0(E) + 
		//"\nsin2Psi = " + ray.sin2Psi() +
		//"\ncos2Psi = " + ray.cos2Psi() +
		"\nsin2Chi = " + Pol.sin2Chi(E) +
		"\ncos2Chi = " + Pol.cos2Chi(E)
		);
	}
	
	public static void testPhaseTangents(){
		double E[] = new double[4];
		
		int nAngs = 4;
		double Eu0 = 10;
		double Er0 = 1;
		double d[][] = new double[nAngs*nAngs][4];
		int k=0;
		for(int i=0; i < nAngs; i++){
			double phsU = i * 2 * Math.PI / nAngs; 
			E[Pol.uRe] = Eu0 * Math.cos(phsU);
			E[Pol.uIm] = Eu0 * Math.sin(phsU);
			
			for(int j=0; j < nAngs; j++){
				double phsR = j * 2 * Math.PI / nAngs;
				E[Pol.rRe] = Er0 * Math.cos(phsR);
				E[Pol.rIm] = Er0 * Math.sin(phsR);
				
				System.out.println("phsU = "+phsU + ", phsR = "+phsR+", phsR-phsU="+(phsR-phsU)+", tan(phsDiff)="+Math.tan((phsR-phsU)) + " ?= " + Pol.tanPhaseDiff(E));
				d[k][0] = phsU;
				d[k][1] = phsR;
				d[k][2] = Math.tan((phsR-phsU));
				d[k][3] = Pol.tanPhaseDiff(E);
				k++;
			}
		}
		BinaryMatrixFile.mustWrite("/tmp/tanTest.bin", d, false);
	}
}
