package fusionOptics.materials;


import fusionOptics.types.Material;

public class OharaS_LAH64 extends DispersionCoefficientBased {
	
	/** [ https://www.ohara-gmbh.com/fileadmin/user_upload/export-data/pdf/product_datasheets/S-LAH64_English_.pdf ] 
	 * Coeffs for visible range.  */
	public static final double coeffs[][] = new double[][]{{
		1.83021453E+00,		//A1
		2.91563590E-01,		//A2
		1.28544024E+00,		//A3
		9.04823290E-03,		//B1
		3.30756689E-02,		//B2
		8.93675501E+01,		//B3
	}};	

	public OharaS_LAH64() {	super(coeffs); }
	
	
	@Override
	public double getVerdetConstant(int modeNumber, double wavelen,
			double temperature) {
		
		throw new RuntimeException("Not implemented.");		
	}

}
