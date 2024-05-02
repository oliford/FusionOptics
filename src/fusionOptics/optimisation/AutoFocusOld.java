package fusionOptics.optimisation;

import java.util.List;

import fusionOptics.Util;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Plane;
import fusionOptics.types.Element;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;
import seed.digeom.IFunction;
import seed.optimization.Optimizer;
import uk.co.oliford.jolu.OneLiners;


/** Optimiser to move given elements along their normals, or first surface normals to achieve best focus
 * 
 * @author oliford
 *
 */
public class AutoFocusOld extends OptimiseOptic {
	
	private Element elems[];
	private double pos0[][];
	private double axis[];
	
	/** Sets a single moveable element, the movement direction and how much */
	public void setMovement(Element elem, double axis[], double min, double max){
		setMovement(new Element[]{ elem }, axis, min, max);
	}
	
	/** Sets the elements that can be moved, the movement direction and how much */
	public void setMovement(Element elems[], double axis[], double min, double max){
		this.elems = elems;
		this.axis = Util.reNorm(axis);
		setParameterSpace(
				OneLiners.fillArray(min, elems.length),
				OneLiners.fillArray(max, elems.length));
		
	}
	
	private void storeInitPositions(){
		int n = elems.length;
		pos0 = new double[n][];

		for(int i=0; i < n; i++){
			pos0[i] = elems[i].getBoundarySphereCentre().clone();
		}
	}
	
	public double[] optimise(Optimizer opt, int nIters) {
		if(elems == null) throw new IllegalArgumentException("No elements to move. Call AutoFocus.setMovement() first.");
		storeInitPositions();		
		return optimise(opt, OneLiners.fillArray(0.0000, elems.length), nIters);
	}

	@Override
	protected void setParams(double[] p) {
		for(int i=0; i < elems.length; i++){
			double c[] = elems[i].getBoundarySphereCentre();
			double oldShift[] = Util.minus(c, pos0[i]);
			//we have to remove the current shift and then reshift by p[i]
			elems[i].shift(new double[]{
					-oldShift[0] + p[i] * axis[0],
					-oldShift[1] + p[i] * axis[1],
					-oldShift[2] + p[i] * axis[2],
			});
		}
	}

	@Override
	/** Not tested */
	protected double[] getParams() {
		double p[] = new double[elems.length];
		for(int i=0; i < elems.length; i++){
			double c[] = elems[i].getBoundarySphereCentre();
			double oldShift[] = Util.minus(c, pos0[i]);
			p[i] = Util.dot(axis, oldShift);
		}
		return p;
	}

}
