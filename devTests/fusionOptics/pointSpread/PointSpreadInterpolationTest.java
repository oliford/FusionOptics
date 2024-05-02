package fusionOptics.pointSpread;

import java.text.DecimalFormat;
import java.util.HashMap;

import fusionOptics.pointSpread.GaussianPSF;
import fusionOptics.pointSpread.PSFGrid;
import fusionOptics.pointSpread.PSFStatsSourceInterpolation;
import uk.co.oliford.cache.randomAccessCache.RACache;
import uk.co.oliford.cache.randomAccessCache.RACacheService;
import uk.co.oliford.jolu.OneLiners;


/** Low level test dev platform for PSF interpolation */
public class PointSpreadInterpolationTest {
	public static void main(String[] args) {
		RACache c = (RACache)RACacheService.instance().getCache("optics.PSF");
		c.emptyAllSets("interpTest");

		GaussianPSF gm00 = new GaussianPSF(0, 9, 9, 9, 9, 9); c.put("interpTest", new double[]{ 0, 0, -1 }, gm00);
		GaussianPSF gm01 = new GaussianPSF(0, 9, 9, 9, 9, 9); c.put("interpTest", new double[]{ 1, 0, -1 }, gm01);
		GaussianPSF gm10 = new GaussianPSF(0, 9, 9, 9, 9, 9); c.put("interpTest", new double[]{ 0, 1, -1 }, gm10);
		GaussianPSF gm11 = new GaussianPSF(0, 9, 9, 9, 9, 9); c.put("interpTest", new double[]{ 1, 1, -1 }, gm11);
		
		GaussianPSF g000 = new GaussianPSF(1, 0, 7, 1, 0, 1); c.put("interpTest", new double[]{ 0, 0, 0 }, g000);
		GaussianPSF g001 = new GaussianPSF(1, 3, 3, 1, 0, 1); c.put("interpTest", new double[]{ 1, 0, 0 }, g001);
		GaussianPSF g010 = new GaussianPSF(1, 6, 2, 1, 0, 1); c.put("interpTest", new double[]{ 0, 1, 0 }, g010);
		GaussianPSF g011 = new GaussianPSF(1, 9, 3, 1, 0, 1); c.put("interpTest", new double[]{ 1, 1, 0 }, g011);
		
		GaussianPSF g100 = new GaussianPSF(2, 2, 5, 4, 0, 1); c.put("interpTest", new double[]{ 0, 0, 1 }, g100);
		GaussianPSF g101 = new GaussianPSF(2, 3, 8, 4, 0, 1); c.put("interpTest", new double[]{ 1, 0, 1 }, g101);
		GaussianPSF g110 = new GaussianPSF(2, 4, 4, 4, 0, 1); c.put("interpTest", new double[]{ 0, 1, 1 }, g110);
		GaussianPSF g111 = new GaussianPSF(2, 5, 1, 4, 0, 1); c.put("interpTest", new double[]{ 1, 1, 1 }, g111);
		
		GaussianPSF gp00 = new GaussianPSF(0.0001, 9, 9, 9, 9, 9); c.put("interpTest", new double[]{ 0, 0, 2 }, gp00);
		GaussianPSF gp01 = new GaussianPSF(0, 9, 9, 9, 9, 9); c.put("interpTest", new double[]{ 1, 0, 2 }, gp01);
		GaussianPSF gp10 = new GaussianPSF(0, 9, 9, 9, 9, 9); c.put("interpTest", new double[]{ 0, 1, 2 }, gp10);
		GaussianPSF gp11 = new GaussianPSF(0, 9, 9, 9, 9, 9); c.put("interpTest", new double[]{ 1, 1, 2 }, gp11);
		
		double gridDef[][] = new double[][]{ { 0, 1, 2 }, { 0, 1, 2 }, { -1, 2, 4 } };
		
		PSFGrid psfGrid = new PSFGrid(c, "interpTest", gridDef, GaussianPSF.class);
				
		PSFStatsSourceInterpolation psfInterp = new PSFStatsSourceInterpolation(psfGrid);
		
		double z[] = OneLiners.linSpace(-0.5, 1.5, 0.1);
		DecimalFormat fmt = new DecimalFormat("##.###");
		for(int i=0; i < z.length; i++){
			GaussianPSF gi = (GaussianPSF) psfInterp.getInterpolatedPSF(0.5, 0.5, z[i]);
			
			System.out.println(fmt.format(z[i]) + "\t"+ (gi == null ? "null" : (gi.I0 + "\t" + gi.cXX)) );
		}
		
	}
}
