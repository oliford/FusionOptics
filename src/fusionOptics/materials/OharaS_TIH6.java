package fusionOptics.materials;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;
import fusionOptics.types.Material;

public class OharaS_TIH6 extends DispersionCoefficientBased {
	
	/** [ www.oharacorp.com/pdf/estih06.pdf ] 
	 * Coeffs for visible range.  */
	public static final double coeffs[][] = new double[][]{{
		1.77227611E+00,		//A1
		3.45691250E-01,		//A2
		2.40788501E+00,		//A3
		1.31182633E-02,		//B1
		6.14479619E-02,		//B2
		2.00753254E+02,		//B3
	}};	

	public OharaS_TIH6() {	super(coeffs); }
	
	
	@Override
	public double getVerdetConstant(int modeNumber, double wavelen,
			double temperature) {
		
		if(wavelen < 640e-9 || wavelen > 670e-9){
			throw new RuntimeException("Not implemented. This was only measured at 656nm");
		}
		
		return 0.4; // We measured somewhere between -0.14 and +0.4, so this is the worst case scenario
	}

}
