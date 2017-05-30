package fusionOptics.materials;

/** YVO4 - Yttrium Orthovanadate.
 * A birefringent crystal that we might use for waveplates and displacers. 
 * 
 */ 
public class YVO4 extends SellmeierBased {

	public YVO4() {
		/* SellmeierBased is:
		 * 
		 * F = (T - T0) * (T + T0 + 546);
		 *	n^2 = A1 + (A2 + B1*F) / (L*L - (A3+B2*F)*(A3+B2*F)) + B3*F - A4*L*L;
		 *
		 * [ http://www.unitedcrystals.com/YVO4Prop.html ]
		 * Matches wikipedia.
		 *
		 */
		super(new double[][]{
				{ //ordinary						
					3.77834, 0.069736, Math.sqrt(0.04724), 0.0108133, //A1,A2,A3,A4
					0, 0, //B1, B2
					0, //B3, we sould be able to do
					27, //T0 / deg.C - Not given by ref, although it does give dn/dT
					
				},{ //extraordinary
					4.59905, 0.110534, Math.sqrt(0.04813), 0.0122676,
					0, 0, //B1, B2 
					0, //B3, we sould be able to do
					27 //T0 / deg.C					
				},
		}, 
		Double.NaN); //magneto-optic anomaly
		
	}
}
