package fusionOptics.optimisation;

import fusionOptics.materials.IsotropicLinearDispersiveGlass;
import fusionOptics.types.Element;
import fusionOptics.types.Material;
import jafama.FastMath;
import seed.optimization.Optimizer;

/** Modification of material refractive indecies to focus an image
 * on specified points (rather than just focusing anywhere as best as pos) 
 * 
 * Useful if you have the drawing of a lens structure without the details
 * of the glasses.
 * 
 * @author oliford
 *
 */
public class OptimiseIndices extends OptimiseOptic {
	private double minN[], maxN[];
	private Material mats[];
	
	private double angleStep;
	private double intendedFocalLength;
	
	public OptimiseIndices(double intendedFocalLength) {
		this.intendedFocalLength = intendedFocalLength;
	}
	
	public void setModifications(Material mats[], double minN[], double maxN[]){
		this.minN = minN;
		this.maxN = maxN;
		this.mats = mats;
		setParameterSpace(minN, maxN);
	}

	@Override
	protected void setParams(double[] p) {
		for(int i=0; i < mats.length; i++){
			//if(mats[i] instanceof IsotropicLinearDispersiveGlass){
				((IsotropicLinearDispersiveGlass)mats[i]).setRefractiveIndexD(p[i]);
			//}
		}
	}

	@Override
	protected double[] getParams() {
		double p[] = new double[mats.length];
		for(int i=0; i < mats.length; i++){
			//if(mats[i] instanceof IsotropicLinearDispersiveGlass){
				p[i] = ((IsotropicLinearDispersiveGlass)mats[i]).getRefractiveIndexD();
			//}
		}
		return p;
	}

	public void optimise(Optimizer opt, int nIters) {
		optimise(opt, getParams(), nIters);
	}
	
	@Override
	public void initRaysParallel(Element rayTarget, double[] rayDir0,
			double maxAngle, double rayLength, double wavelength, int nSets,
			int nRaysPerSet) {
		super.initRaysParallel(rayTarget, rayDir0, maxAngle, rayLength, wavelength,
				nSets, nRaysPerSet);
		this.angleStep = FastMath.tan(maxAngle) / (nSets - 1);
	}
	
	@Override
	protected double extraCostFunction() {
		return super.extraCostFunction() +  0.003*(1.0 - ((double)nRaysInImageLast)/nRaysTotal);
	}
	
	@Override
	protected double[] targetImagePos() {
		return new double[]{ setNo * angleStep * intendedFocalLength, 0 };
	}
	
}
