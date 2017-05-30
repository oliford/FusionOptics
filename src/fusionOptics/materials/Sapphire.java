package fusionOptics.materials;

public class Sapphire extends DispersionCoefficientBased {
	
	/** [ http://refractiveindex.info/?shelf=main&book=Al2O3&page=Malitson ] 
	 * Coeffs for visible range.  */
	public static final double coeffs[][] = new double[][]{{
		1.023798,		//A1
		1.058264,		//A2
		5.280792,		//A3
		0.00377588,		//B1
		0.01225442,		//B2
		321.36155343,		//B3	
	}};	

	public Sapphire() {	super(coeffs); }
	
	
	@Override
	public double getVerdetConstant(int modeNumber, double wavelen, double temperature) {
		throw new RuntimeException("Not implemented.");
		
	}

}
