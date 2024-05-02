package fusionOptics.materials;

import java.util.Arrays;

import fusionOptics.types.Material;

import net.jafama.FastMath;

/** Some info from John of what the Ferro-electric Liquid Crystal dispersion might look like 
 * We really don't know
 * function FLC_G_sell, lambda, n_e=n_e, n_o=n_o, dnedl=dnedl, dnodl=dnodl, dmudl=dmudl,  T = T
 * [ Wu, S.T. 1986 Phys. Rev. A, 33, 1270 ]



return, mu
end

;--------------------------------------------------------------------
function FLC_G, thickness, lambda, biref=biref, $
                 n_e=n_e, n_o=n_o, kappa=kappa
;+
; return wave delay caused by birefringent LiTaO3 crystal
;-

    biref=FLC_G_sell(lambda, n_e=n_e, n_o=n_o, dmudl=dmudl, T = T)
    kappa = lambda/biref*dmudl

return, biref*thickness/lambda

end

 * 
 * */ 
public class LiquidCrystal extends Material {
	public static final double G = 2.117e+12; // m-2 , 2.117x10-6 nm-2;
	public static final double lstar = 253e-9; // not a clue, sorry
	public static final double fixedNO = 1.49;
	@Override
	public int getNAxes() {
		return 1;
	}

	@Override
	public double getRefractiveIndex(int modeNumber, double wavelength,
			double temperature) {
		// should let n have same functional dependence as birefringence
		
		if(modeNumber == 0)
			return fixedNO; // at 593nm
		 
		// birefringence
		double mu = G*wavelength*wavelength*lstar*lstar / (wavelength*wavelength - lstar*lstar);

		// birefringence dispersion
		//double dμdλ = -2*G*FastMath.pow(lstar,4)*wavelength / FastMath.pow(wavelength*wavelength - lstar*lstar, 2);

		return fixedNO + mu;
		
		//no=1.49;
		//ne=no + 2.117e+12*l^2*(253e-9)^2 / (l^2 - (253e-9)^2);
		
		
	}

	@Override
	public double getTransmission(int modeNumber, double wavelength,
			double temperature) {
		return 1;
	}

	@Override
	public double getVerdetConstant(int modeNumber, double wavelen,
			double temperature) {
		throw new UnsupportedOperationException();
	}
	
	//Materials with no modifiable content are equal if they are the same type
	// but the static variables may have changed at compile and the comparison of hashcode might be against an old
	// code
	@Override
	public int hashCode() { 
		long t; int r = 1;
		t = Double.doubleToLongBits(G); 		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(lstar); 		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(fixedNO); 		r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	
	@Override
	public boolean equals(Object obj) { return (obj != null) && obj instanceof LiquidCrystal; }	

}
