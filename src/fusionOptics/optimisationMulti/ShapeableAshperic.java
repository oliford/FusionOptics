package fusionOptics.optimisationMulti;

import java.util.ArrayList;

import fusionOptics.surfaces.Aspheric;

public class ShapeableAshperic extends Parameter {
	
	/** The aspheric to reshape */
	private Aspheric aspheric;
	
	/** Which polyCoeff to modify.
	 * -1 is the conic constant, 0 is the R^2, 1 is R^4 etc....
	 */
	private int order;
	
	/** Creates optimisation parameters for curvature radius, conic constant and each polynomial coefficient >= R^4 */
	public static ArrayList<Parameter> allAsphericParams(Aspheric aspheric){
		int nPoly = aspheric.getPolyCoeffs().length; 
		
		ArrayList<Parameter> params = new ArrayList<Parameter>(nPoly + 2);
		
		params.add(new BendableDish(aspheric));
		
		params.add(new ShapeableAshperic(aspheric, -1, -200.0, 200.0));
		
		double scale = 10;
		for(int i=1; i < nPoly; i++){
			ShapeableAshperic polyParam = new ShapeableAshperic(aspheric, i, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
			polyParam.setScale(scale);
			params.add(polyParam);			
			scale *= 100;
		}
		return params;
	}
	
	/**
	 * @param aspheric
	 * @param order -1 is the conic constant, 0 is the R^2, 1 is R^4 etc....
	 * @param min
	 * @param max
	 */
	public ShapeableAshperic(Aspheric aspheric, int order, double min, double max) {
		super(min, max); //parameter is relative to the initial bounding sphere centre, so starts at 0
		this.aspheric = aspheric;
		this.order = order;
		this.initVal = get();
	}
	
	public ShapeableAshperic(Aspheric aspheric, int order, double min, double max, double scale) {
		super(min, max); //parameter is relative to the initial bounding sphere centre, so starts at 0
		this.aspheric = aspheric;
		this.order = order;
		this.initVal = get();
		this.setScale(scale);
	}
	
	public ShapeableAshperic(Aspheric aspheric, int order, double fractionChange) {
		this.order = order;
		this.aspheric = aspheric; 
		this.initVal = get();
		this.min = (1.0 - fractionChange) * initVal;
		this.max = (1.0 + fractionChange) * initVal;
	}
	
	@Override
	public String toString() { 
		return "ShapeableAshperic["+aspheric.getName()+"]: order=" + order +", init = " + initVal() + ", current="+get()+", min = " + min() + ", max = " + max() + ", scale = " + scale;
	}
	
	@Override
	public double get() { //distance from init pos along the given vector
		return ((order == -1) ? aspheric.getConicConstant() : aspheric.getPolyCoeffs()[order]) / scale;
	}
	
	@Override
	public void set(double val) {
		if(order == -1){
			aspheric.setConicConstant(val*scale);
		}else{
			double polyCoeffs[] = aspheric.getPolyCoeffs();
			polyCoeffs[order] = val*scale;
			aspheric.setPolyCoeffs(polyCoeffs);
		}
	}

}
