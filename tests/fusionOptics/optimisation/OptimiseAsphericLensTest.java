package fusionOptics.optimisation;

import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.materials.BK7;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.materials.SchottSFL6;
import fusionOptics.optics.DoubleGaussLens;
import fusionOptics.optimisation.OptimiseAsphericLens;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.types.Element;
import fusionOptics.types.Interface;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;
import seed.optimization.HookeAndJeeves;
import junit.framework.TestCase;

/** Test if the aspherical surface optimisation
 * can bring the a single aspherical surface lens into sharp focus at high f number.
 * 
 * 
 * There's som
 * 
 * @author oliford
 */
public class OptimiseAsphericLensTest extends TestCase {
	private static String outPath = System.getProperty("java.io.tmpdir") + "/rayTracing/tests";
	
	public void testAsphericalOptimisation(){
		
		double f = 1.0;	//approx focal length
		double wavelen = 550e-9; //greenish
		
		Material lensMat1 = new IsotropicFixedIndexGlass(1.8);
		Medium lensMed1 = new Medium(lensMat1);
		Material lensMat2 = new IsotropicFixedIndexGlass(1.8);
		Medium lensMed2 = new Medium(lensMat2);
		
		Aspheric lensAspheric1 = new Aspheric("lensAspheric1",
				new double[]{ 0.0, 0, 0 }, //centre of front surface
				new double[]{ 1, 0, 0 }, //centre normal
				0.7, //rad curv
				0.6, //rim radius
				0.0, //conic const (NOT SUPPORTED, must be 0)
				new double[]{ 0,0,0,0,0,0 }, //polynomial coeffs
				lensMed1, null, IsoIsoInterface.ideal());
		
		/*Surface s2 = new Aspheric("s2",
				new double[]{ 0.7, 0, 0 }, //centre of front surface
				new double[]{ 1, 0, 0 }, //centre normal
				4.0, //rad curv
				0.6, //rim radius
				0.0, //conic const (NOT SUPPORTED, must be 0)
				new double[]{ 0,0 }, //polynomial coeffs
				null, lensMed1, IsoIsoInterface.ideal());
				*/
		
		Surface s2 = new Disc("s2",
				new double[]{ 0.6, 0, 0 }, //centre of front surface
				new double[]{ 1, 0, 0 }, //centre normal
				4.0, //rad curv
				null, lensMed1, IsoIsoInterface.ideal());
		
		Iris apature = new Iris("aparture",
				new double[]{ 0,0,0}, 
				new double[]{1,0,0}, 
				0.8, 
				0.5, 
				Absorber.ideal()); 
		
		
		
		Square imgPlane = new Square("imgPlane", new double[]{ f, 0, 0 }, new double[]{ 1,0,0 }, new double[]{ 0, 0, 1}, 0.5, 1.0, null, null, Absorber.ideal()); 
		
		Optic all = new Optic("all", new Element[]{ apature, lensAspheric1, s2, imgPlane });
		
		OptimiseAsphericLens optim = new OptimiseAsphericLens();
				
		optim.setTracingElements(all, imgPlane);
		//optim.initRaysImaging(lensAspheric, new double[]{0,0,0}, new double[]{0,0,0.002}, 10, wavelen, 10000);
		optim.initRaysParallel(lensAspheric1, new double[]{1,0,0}, 10*Math.PI/180, 0.10, wavelen, 15, 1000);
		optim.setModifications(
				new Element[]{ lensAspheric1, s2 },
				0.2, 5.0, 0.9, new double[]{ 0.2, 0, 0 });
		//System.out.println(optim.eval(optim.getParams()));
		optim.setHitsOutputFile(outPath + "/optimHits.bin");
		optim.optimise(new HookeAndJeeves(null), 100);
			
		//optim.singleParamScan(0, 1.0, 4.0, 30);
		/*optim.setParams(new double[]{ 
		 	0.669844843917371,	0.1	,-0.09926072892582444,
		 	0.6897939916008341,	0.07942064264878274	,0.06590931819863885,	
		 	0.8996352764981579
			});
			//*/
		
		//optim.setHitsOutputFile(outPath + "/optimHits.bin");
		optim.setHitsOutputFile(null);
		optim.setSVGOut(outPath + "/asphericOptim", 77);
		optim.optimise(new HookeAndJeeves(null), 0);
		
		double fwhm = optim.eval(optim.getParams());
			
		// Without the aspherics (setting no polyCoeffs), this gets down to about 5.7% (of f)  
		// which isn't too bad for a f/1 lens. 
		// With 4th order aspherics, we will require < 3%. It's still not great, but it'll do
		
		//actually, it really depends so heavily on the init conditions, that I'd not like to say what
		//difference the aspherics make.
		
		//but.. It can get down to 1.5%, and never really gets better than that.
		assertEquals(0, fwhm, f * 0.03);
	}

}
