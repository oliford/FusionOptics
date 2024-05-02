package fusionOptics.optimisation;

import fusionOptics.Util;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Plane;
import fusionOptics.types.Element;
import seed.digeom.IFunction;
import seed.optimization.Optimizer;
import uk.co.oliford.jolu.OneLiners;


/** Optimiser to set the best curvature radii for a simple lens  
 * 
 * @author oliford
 *
 */
public class OptimiseLensCurvature extends OptimiseOptic {
	
	private SimpleDoubleConvexLens lenses[];
	/** The two surfaces of the lens */
	private Dish[] s;

	public OptimiseLensCurvature(SimpleDoubleConvexLens lenses[], Plane imagePlane, 
				Element all, double rayStart[], double wavelength, int nRays) {
		if(true)throw new RuntimeException("Not updates to newer OptimiseOptic, see AutoFocus for new architecture");
		/*super(imagePlane, 
				OneLiners.fillArray(0, lenses.length*2), 
				OneLiners.fillArray(Double.POSITIVE_INFINITY, lenses.length*2),
				all, lenses[0], rayStart, wavelength, nRays);
			*/
		this.lenses = lenses;
		s = new Dish[lenses.length*2];
		for(int i=0; i < lenses.length; i++){
			s[i*2+0] = ((Dish)lenses[i].getSurfaces().get(0));
			s[i*2+1] = ((Dish)lenses[i].getSurfaces().get(1));
		}
		
	}
	
	public OptimiseLensCurvature(SimpleDoubleConvexLens lenses[], Plane imagePlane, 
				Element all, double rayDir[], double rayLength, double wavelength, int nRays) {
		if(true)throw new RuntimeException("Not updates to newer OptimiseOptic, see AutoFocus for new architecture");
		/*super(imagePlane, 
				OneLiners.fillArray(0, lenses.length*4-1), 
				OneLiners.fillArray(Double.POSITIVE_INFINITY, lenses.length*4-1),
				all, lenses[0], rayDir, rayLength, wavelength, nRays);
			*/
		this.lenses = lenses;
		s = new Dish[lenses.length*2];
		for(int i=0; i < lenses.length; i++){
			s[i*2+0] = ((Dish)lenses[i].getSurfaces().get(0));
			s[i*2+1] = ((Dish)lenses[i].getSurfaces().get(1));
		}
		
	}
	
	public double[] optimise(Optimizer opt, int nIters) {
		return optimise(opt, getParams(), nIters);
	}

	@Override
	protected void setParams(double[] p) {
		s[0].setRadiusOfCurv(p[0]);		
		double c0[] = s[0].getCentre(); 
		for(int i=1; i < s.length; i++){
			double d[] = Util.reNorm(Util.minus(s[i].getCentre(), s[0].getCentre()));
			s[i].setRadiusOfCurv(p[i*2-1]);
			s[i].setCentre(new double[]{
					c0[0] + p[i*2-0] * d[0],
					c0[1] + p[i*2-0] * d[1],
					c0[2] + p[i*2-0] * d[2],
			});
		}
	}

	@Override
	protected double[] getParams() {
		double p[] = new double[s.length*2-1]; 
		p[0] = s[0].getRadiusOfCurv();		
		for(int i=1; i < s.length; i++){
			p[i*2-1] = s[i].getRadiusOfCurv();
			p[i*2-0] = Util.length(Util.minus(s[i].getCentre(), s[0].getCentre()));
		}
		return p;
	}

}
