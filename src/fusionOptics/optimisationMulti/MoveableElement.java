package fusionOptics.optimisationMulti;

import fusionOptics.Util;
import fusionOptics.types.Element;

/** Elements we can move */
public class MoveableElement extends Parameter {
	/** The element to move */
	private Element element;
	/** Direction of shift */
	private double vector[];
	/** Position of boundingSphereCentre() corresponding to shift=0  */
	private double zeroPos[];
	
	
	public MoveableElement(Element element, double vector[], double min, double max) {
		this(element, vector, min, max, 1.0);
	}
	
	public MoveableElement(Element element, double vector[], double min, double max, double scale) {
		super(min, max, 0); //parameter is relative to the initial bounding sphere centre, so starts at 0
		
		this.element = element;
		this.zeroPos = element.getCentre().clone();
		this.vector = Util.reNorm(vector);
		this.scale = scale;
	}
	
	@Override
	public String toString() { 
		double c[] = element.getCentre();
		return "MoveableElement["+element.getName()+"]: init = " + initVal() + ", current="+get()+", min = " + min() + ", max = " + max() + ", scale = " + scale + ", pos={ "+c[0]+","+c[1]+","+c[2]+"}";
	}
	
	@Override
	/** Get position relative to zeroPos[], along vector */
	public double get() { //distance from init pos along the given vector
		return Util.dot(Util.minus(element.getCentre(), zeroPos), vector) / scale;
	}
	
	@Override
	/** Move to new position relative to zeroPos[] */
	public void set(double val) {
		double current = get();
		double shift = (val - current) * scale;
		element.shift(Util.mul(vector, shift));
	}
	
	
	
}
