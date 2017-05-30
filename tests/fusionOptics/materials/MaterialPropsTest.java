package fusionOptics.materials;

import fusionOptics.materials.AlphaBariumBorate;
import fusionOptics.materials.BK7;
import fusionOptics.materials.BetaBariumBorate;
import fusionOptics.materials.CrystalQuartz;
import fusionOptics.materials.FusedSilica;
import fusionOptics.materials.IsotropicLinearDispersiveGlass;
import fusionOptics.materials.LithiumNiobate;
import fusionOptics.materials.SchottSFL6;
import fusionOptics.types.Material;
import junit.framework.TestCase;

/** Check that various material properties at some specific
 * wavelengths and sizes etc, match various bits of literature 
 * (not the literature that the material data came from)
 * 
 * Currently,  almost nothing actually passes :(
 * 
 * @author oliford
 *
 */
public class MaterialPropsTest extends TestCase {
	
	public void testOharaTIH6(){
		OharaS_TIH6 tih6 = new OharaS_TIH6();
		assertEquals(1.77495, tih6.getRefractiveIndex(0, 1013.98e-9, 300), 0.00001); // nt 
		assertEquals(1.79611, tih6.getRefractiveIndex(0, 656.27e-9, 300), 0.00001); // nC 
		assertEquals(1.80491, tih6.getRefractiveIndex(0, 589.29e-9, 300), 0.00001); //  nD 
		assertEquals(1.80518, tih6.getRefractiveIndex(0, 587.56e-9, 300), 0.00001); // nd 
		assertEquals(1.84729, tih6.getRefractiveIndex(0, 435.835e-9, 300), 0.00001); // ng 
		assertEquals(1.86494, tih6.getRefractiveIndex(0, 404.656e-9, 300), 0.00001); // nh
	}
	
	/** This currently doesn't pass */
	public void testFusedSilicaDispersion(){
		//refractiveIndex.info
		//assertEquals("dispersionQuartz", -0.0333e6, new Quartz().getLinearDispersion(0, 600e-9, 300)); handbook of optics, quartz?
		//assertEquals("dispersionQuartz", -0.0162e6, new Quartz().getLinearDispersion(0, 600e-9, 300)); // SOPRA, quartz?
		//from differentiating myself from 
		// [ http%3A%2F%2Fwww.lightmachinery.com%2FMaterials%2Fstandard%2520fused%2520silica%25207980.pdf&ei=y8E9T4KOFMzCmQWI-aHEBw&usg=AFQjCNECBVFAMiD4rngs4GT387lvHy68RA ]
		assertEquals("dispersionQuartz", -0.0390e6, new FusedSilica().getLinearDispersion(0, 600e-9, 300)); 
	}
	
	/** This currently doesn't pass */
	public void testFusedSilicaVerdetConst(){
		/* Original data is from Bequerel formula, an entirely alternate measurement
		 * comes from
		 * [ C.Z.Tan "Wavelength dependence of the Faraday effect in glassy SiO2"
				J Phys & Chem Solids 60 (1999) 1689–1692 ]
		
				The numbers given are radians per Amp in their glass rod, in their coil
				with H = 18484.3 I, presumably in S.I units of H, so				
				B = mu0 18484.3 I.
				length of rod = 19.97mm
				
				Err... no.
				They say their V = -e µ0/(2mc)λ dN/dλ, which is same as the Bequerel forumla but with the µ_0.
				That converts from their odd radians/Amp unit back to radians/Tesla/m
				The experiment geometry is not in the number, just the silly unit.
				
		 */
		double mu0 = 1.25663706e-6;
		double TtoA = mu0; // * 18484.3;
		double l = 0.01997; //19.97mm long bit of glass
		
		FusedSilica SiO2 = new FusedSilica();
		
		//assertEquals("verdetQuartz(l=500nm)", 6.5e-6 * TtoA / l, SiO2.getVerdetConstant(0, 500e-9, 300));
		//assertEquals("verdetQuartz(l=600nm)", 4.7e-6 * TtoA / l, SiO2.getVerdetConstant(0, 600e-9, 300));
		//assertEquals("verdetQuartz(l=750nm)", 2.5e-6 * TtoA, SiO2.getVerdetConstant(0, 750e-9, 300));		
		assertEquals("verdetQuartz(l=600nm)", 4.7e-6 / mu0, SiO2.getVerdetConstant(0, 600e-9, 300));
		
	}

	public void testABBORefractiveIndecies(){
		/* Original data from
		 * [ http://www.unitedcrystals.com/BBOProp.html United Crystals website ]
		   [ http://www.agoptics.com/Birefrigent-Crystal/alpha-BBO.htm ]
		   [ ... via John Howard ANU ]
		   
		   Test points from [ http://refractiveindex.info/?group=CRYSTALS&material=BBO ]
		 */
		AlphaBariumBorate BBO = new AlphaBariumBorate();
		
			
	}

	public void testBBBORefractiveIndecies(){
		BetaBariumBorate BBO = new BetaBariumBorate();

		//assertEquals("betaBBO no(l=220nm)", 1.8284, BBO.getRefractiveIndex(0, 229e-9, 300)); FAILS
		//assertEquals("betaBBO no(l=600nm)", 1.55994, BBO.getRefractiveIndex(0, 600e-9, 300)); //from refractiveIndex.info, FAILS
		//nope, they don't work
					
	}
	
	public void testAbbeBasedGlass(){
		BK7 bk7 = new BK7();
		
		double nd = bk7.getRefractiveIndex(0, IsotropicLinearDispersiveGlass.lambda_d, 300);
		double nF = bk7.getRefractiveIndex(0, IsotropicLinearDispersiveGlass.lambda_F, 300);
		double nC = bk7.getRefractiveIndex(0, IsotropicLinearDispersiveGlass.lambda_C, 300);
		
		double Vd = (nd - 1.0) / (nF - nC);
		System.out.println("Vd = " + Vd);
		
		IsotropicLinearDispersiveGlass g = new IsotropicLinearDispersiveGlass(nd, Vd, 1.0);

		System.out.println(nd + "\t" + g.getRefractiveIndex(0, IsotropicLinearDispersiveGlass.lambda_d, 300));
		System.out.println(nF + "\t" + g.getRefractiveIndex(0, IsotropicLinearDispersiveGlass.lambda_F, 300));
		System.out.println(nC + "\t" + g.getRefractiveIndex(0, IsotropicLinearDispersiveGlass.lambda_C, 300));

		//this is very approximate, because it's not really linear
		assertEquals(nd, g.getRefractiveIndex(0, IsotropicLinearDispersiveGlass.lambda_d, 300), 1e-8);
		assertEquals(nF, g.getRefractiveIndex(0, IsotropicLinearDispersiveGlass.lambda_F, 300), 0.001);
		assertEquals(nC, g.getRefractiveIndex(0, IsotropicLinearDispersiveGlass.lambda_C, 300), 0.001);
		
	}
	
	public void testNothingReally(){
		System.out.println("Quartz: " + (new CrystalQuartz()).getVerdetConstant(0, 589.3e-9, 300) * 180 / Math.PI * 0.1 * 1e-3 +
											" degrees per mm at 100mT");
		//System.out.println("LiNbO3: " + (new LithiumNiobate()).getVerdetConstant(0, 589.3e-9, 300) * 180 / Math.PI * 0.1 * 1e-3 +
		//									" degrees per mm at 100mT");
		System.out.println("SFL6: " + (new SchottSFL6()).getVerdetConstant(0, 589.3e-9, 300) * 180 / Math.PI * 0.1 * 1e-3 +
											" degrees per mm at 100mT");
		System.out.println("BK7: " + (new BK7()).getVerdetConstant(0, 589.3e-9, 300) * 180 / Math.PI * 0.1 * 1e-3 +
											" degrees per mm at 100mT");
		System.out.println("LiNbO3 No(589.3nm) = " + (new LithiumNiobate()).getRefractiveIndex(0, 589.3e-9, 300));
		System.out.println("LiNbO3 No(589.4nm) = " + (new LithiumNiobate()).getRefractiveIndex(0, 589.4e-9, 300));
		System.out.println("SFL6 No(589.3nm) = " + (new SchottSFL6()).getRefractiveIndex(0, 589.3e-9, 300));
		System.out.println("SFL6 No(589.4nm) = " + (new SchottSFL6()).getRefractiveIndex(0, 589.4e-9, 300));
		System.out.println("BK7 No(589.4nm) = " + (new BK7()).getRefractiveIndex(0, 589.4e-9, 300));
		
		//stuff for magneton exps
		System.out.println("LiNbO3 No(491nm) = " + (new LithiumNiobate()).getRefractiveIndex(0, 491e-9, 300));
		System.out.println("LiNbO3 Ne(491nm) = " + (new LithiumNiobate()).getRefractiveIndex(1, 491e-9, 300));
		System.out.println("LiNbO3 No(650nm) = " + (new LithiumNiobate()).getRefractiveIndex(0, 650e-9, 300));
		System.out.println("LiNbO3 Ne(650nm) = " + (new LithiumNiobate()).getRefractiveIndex(1, 650e-9, 300));
		System.out.println("LiNbO3 No(720nm) = " + (new LithiumNiobate()).getRefractiveIndex(0, 720e-9, 300));
		System.out.println("LiNbO3 Ne(720nm) = " + (new LithiumNiobate()).getRefractiveIndex(1, 720e-9, 300));
		
		//stuff for imse design
		Material mat = new BetaBariumBorate();
		System.out.println("aBBO No(653.5nm) = " + mat.getRefractiveIndex(0, 653.5e-9, 300));
		System.out.println("aBBO Ne(653.5nm) = " +mat.getRefractiveIndex(1, 653.5e-9, 300));
		System.out.println("aBBO No(654.5nm) = " + mat.getRefractiveIndex(0, 654.5e-9, 300));
		System.out.println("aBBO Ne(654.5nm) = " +mat.getRefractiveIndex(1, 654.5e-9, 300));
		System.out.println("aBBO dNo/dlambda(653.5nm) = " + mat.getLinearDispersion(0, 653.5e-9, 300));
		System.out.println("aBBO dNe/dlambda(653.5nm) = " + mat.getLinearDispersion(1, 653.5e-9, 300));
		System.out.println("aBBO B(653.5nm) = " + (mat.getRefractiveIndex(1, 653.5e-9, 300) - mat.getRefractiveIndex(0, 653.5e-9, 300)));
		System.out.println("aBBO dB/dlambda(653.5nm) = " + (mat.getLinearDispersion(1, 653.5e-9, 300) - mat.getLinearDispersion(0, 653.5e-9, 300)));
		
		System.out.println("aBBO No(464nm) = " + mat.getRefractiveIndex(0, 464e-9, 300));
		System.out.println("aBBO Ne(464nm) = " +mat.getRefractiveIndex(1, 464e-9, 300));
		
		//phs(464.742e-9);
		//phs(465.025e-9);
		//phs(465.147e-9);
		phs(653.9e-9);
		phs(653.5e-9);
		phs(653.1e-9);
		
		//(new AlphaBariumBorate()).getLinearDispersion(0, 600e-9, 300);
	}
	
	private void phs(double wavelen) {
		Material mat = new BetaBariumBorate();
		
		System.out.println("aBBO dPhs(wavelen) = " + (mat.getRefractiveIndex(1, wavelen, 300) - mat.getRefractiveIndex(0, wavelen, 300)) * 2.7e-3 / wavelen);
		
		//mat = new LithiumNiobate();
		//System.out.println("LiNi dPhs(wavelen) = " + (mat.getRefractiveIndex(1, wavelen, 300) - mat.getRefractiveIndex(0, wavelen, 300)) * 2.7e-3 / wavelen);
		
	}
	
	public void testSapphire() {
		Sapphire sapphire = new Sapphire();
		
		assertEquals(1.7866, sapphire.getRefractiveIndex(0, 400.0e-9, 300), 0.0001); // nt 
		assertEquals(1.7682, sapphire.getRefractiveIndex(0, 587.6e-9, 300), 0.0001); // nt 
		assertEquals(1.7633, sapphire.getRefractiveIndex(0, 700.0e-9, 300), 0.0001); // nt 
		
	}
}
