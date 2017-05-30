package fusionOptics.optimisationMulti;

/** Base class of all parameters optimised by OptimiseMulti */
public abstract class Parameter {
	protected double min, max;
	protected double initVal;
	protected double scale = 1.0;
	
	public Parameter(double min, double max) {
		this.min = min;
		this.max = max;
		this.initVal = Double.NaN;
	}
	
	public Parameter(double min, double max, double initVal) {
		this.min = min;
		this.max = max;
		this.initVal = initVal;
	}
	
	public Parameter() { //must be done by derived class constructor
		this.min = Double.NaN;
		this.max = Double.NaN;
		this.initVal = Double.NaN;
	}

	public abstract double get();
	public abstract void set(double val);
	
	public double min(){ return min; }
	public double max(){ return max; }
	public double initVal(){ return initVal; }
	
	/** Scaling of parameters (done by derived classes) - this should be done in Digeom land, but that support was in the Minerva/Digeom gateway */
	public void setScale(double scale) {
		min *= this.scale / scale;
		max *= this.scale / scale;
		initVal *= this.scale / scale;
		this.scale = scale;
	}
	
	public double getScale() { return scale; }
}
