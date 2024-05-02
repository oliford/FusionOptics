package fusionOptics.optimisation;

import fusionOptics.Util;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Plane;
import fusionOptics.types.Element;
import seed.digeom.IFunction;
import seed.optimization.Optimizer;
import uk.co.oliford.jolu.OneLiners;


/** Optimiser to bend a single apsherical surface to get the best focus.  
 * @deprecated Superceeded by OptimiseMulti and ShapableAspheric
 * 
 * @author oliford
 */
public class OptimiseAsphericLens extends OptimiseOptic {

	public Element elem[];
	int nParams = 0;
	
	double zeroPos[];
	double shiftAxis[];

	public void setModifications(Element elem[], double minRadCurv, double maxRadCurv, double maxCoeffs, double maxShift[]){
		this.elem = elem;
		for(int i=0; i < elem.length; i++){
			if(elem[i] instanceof Aspheric)
				nParams += ((Aspheric)elem[i]).getPolyCoeffs().length + 1;
			
			if(i > 0)
				nParams++;
					
		}
		
		double maxShiftLen = Util.length(maxShift);
		shiftAxis = Util.reNorm(maxShift.clone());
		
		double pMin[] = new double[nParams];
		double pMax[] = new double[nParams];
		
		int k = 0;
		for(int i=0; i < elem.length; i++){
			if(elem[i] instanceof Aspheric){
				pMin[k] = minRadCurv;
				pMax[k] = maxRadCurv;
				k++;
				int nPoly = ((Aspheric) elem[i]).getPolyCoeffs().length;
				for(int j=0; j < nPoly; j++){
					pMin[k] = -maxCoeffs;
					pMax[k] = maxCoeffs;
					k++;				
				}
			}
			
			if(i == 0){
				zeroPos = elem[i].getBoundarySphereCentre();
			}else{
				double currentShift = Util.dot(Util.minus(elem[i].getBoundarySphereCentre(), zeroPos), shiftAxis);
				pMin[k] = currentShift - maxShiftLen;
				pMax[k] = currentShift + maxShiftLen;
				k++;
			}
		}		
		setParameterSpace(pMin, pMax);
	}
		
	public double[] optimise(Optimizer opt, int nIters) {
		return optimise(opt, getParams(), nIters);
	}
	
	public double[] getParams(){		
		double p[] = new double[nParams];
		int k = 0;
		for(int i=0; i < elem.length; i++){
			if(elem[i] instanceof Aspheric){
				p[k++] = ((Aspheric)elem[i]).getRadiusOfCurv();
				
				double c[] = ((Aspheric) elem[i]).getPolyCoeffs();
				System.arraycopy(c, 0, p, k, c.length);
				k += c.length;
			}
			if(i > 0){
				p[k++] = Util.dot(Util.minus(elem[i].getBoundarySphereCentre(), zeroPos), shiftAxis);				
			}
		}
		return p;
	}

	@Override
	protected void setParams(double[] p) {
		
		int k = 0;
		for(int i=0; i < elem.length; i++){
			if(elem[i] instanceof Aspheric){
				((Aspheric) elem[i]).setRadiusOfCurv(p[k++]);
				
				int nPoly = ((Aspheric) elem[i]).getPolyCoeffs().length;
				double c[] = new double[nPoly];
				System.arraycopy(p, k, c, 0, nPoly);
				((Aspheric) elem[i]).setPolyCoeffs(c);
				k += c.length;
			}
			
			if(i > 0){
				double currentPos[] = elem[i].getBoundarySphereCentre();
				elem[i].shift(new double[]{
						zeroPos[0] + p[k] * shiftAxis[0] - currentPos[0],
						zeroPos[1] + p[k] * shiftAxis[1] - currentPos[1],
						zeroPos[2] + p[k] * shiftAxis[2] - currentPos[2],
				});
				k++;				
			}
		}	
		
	}

}
