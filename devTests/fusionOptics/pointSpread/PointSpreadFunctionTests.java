package fusionOptics.pointSpread;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.pointSpread.GaussianPSF;
import fusionOptics.pointSpread.PointSpreadFunction;
import uk.co.oliford.jolu.BinaryMatrixFile;
import uk.co.oliford.jolu.OneLiners;
import uk.co.oliford.jolu.RandomManager;

/** Low level test dev platform for Point spread functions */ 
public class PointSpreadFunctionTests {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/psfTest";
	
	public static void main(String[] args) {
		int n = 10000;
		double cosQ = Math.cos(30*Math.PI/180);
		double sinQ = Math.sqrt(1.0 - cosQ*cosQ);
		double pos[][] = new double[n][2];
		double E[][][] = new double[n][1][];
		for(int i=0; i < n ; i++){
			double a = 1.5 * RandomManager.instance().nextNormal(0, 1);
			double b = 0.3 * RandomManager.instance().nextNormal(0, 1);
			pos[i][0] = a * cosQ - b * sinQ;
			pos[i][1] = a * sinQ + b * cosQ;
			E[i][0] = new double[]{ 1,0,0,0 };
		}
		BinaryMatrixFile.mustWrite(outPath + "/pts1.bin", pos, false);
		
		//GaussianPSF psf1 = new GaussianPSF();
		//PointsPSF psf1 = new PointsPSF();
		//psf1.setPoints(pos, E);
		
		GaussianPSF psf1 = new GaussianPSF(2, -1, 2, 0.4, 0.0, 0.005);
//		GaussianPSF psf1 = new GaussianPSF(2, -1, 2, 0.005, 0.063, 0.8);
		GaussianPSF psf2 = new GaussianPSF(2, 5, -4, 0.005, 0.0*63, 0.8);
		GaussianPSF psf3 = new GaussianPSF();
		
		System.out.println(psf1.toString());
		
		double x[] = OneLiners.linSpace(-2, 7, 101);
		double y[] = OneLiners.linSpace(-6, 4, 99);
		
		double P[][] = new double[x.length][y.length];
		psf1.addToGridMonteCarlo(x, y, P, 1, 1000);
		psf2.addToGridMath(x, y, P, 1);
		
		for(int i=0; i < 6; i++){
			psf3.combine(new PointSpreadFunction[]{ psf1, psf2 }, new double[]{ (i/5.0), 1.0-(i/5.0) });
			//spsf3.addToGridMonteCarlo(x, y, P, 1);			
		}
		
		System.out.println("nX = " + P.length + " ?= " + x.length);
		System.out.println("nY = " + P[0].length + " ?= " + y.length);
		
		BinaryMatrixFile.mustWrite(outPath + "/gaussGrid.bin", y, x, P, false);
		

		BinaryMatrixFile.mustWrite(outPath + "/pts2.bin", psf1.generatePoints(1000), true);
		
	}
}
