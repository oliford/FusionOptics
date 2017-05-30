package fusionOptics.optimisationMulti;

import fusionOptics.Util;
import fusionOptics.surfaces.Dish;
import fusionOptics.types.Element;
import fusionOptics.types.Medium;

public class BendableDish extends Parameter {
	
	/** The dish to bend */
	private Dish dish;
	private double initDishNormal[];
	private Medium initFrontMedium, initBackMedium;
	
	public BendableDish(Dish dish) {
		this(dish, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
	}
	
	public BendableDish(Dish dish, double min, double max) {
		super(min, max); //parameter is relative to the initial bounding sphere centre, so starts at 0
		this.dish = dish;
		this.initDishNormal = dish.getDishNormal().clone();
		this.initFrontMedium = dish.getFrontMedium();
		this.initBackMedium = dish.getBackMedium();
		this.initVal = get();
	}
	
	public BendableDish(Dish dish, double min, double max, double scale) {
		super(min, max); //parameter is relative to the initial bounding sphere centre, so starts at 0
		this.dish = dish;
		this.initDishNormal = dish.getDishNormal().clone();
		this.initFrontMedium = dish.getFrontMedium();
		this.initBackMedium = dish.getBackMedium();
		this.initVal = get();
		this.setScale(scale);
	}
	
	public BendableDish(Dish dish, double fractionChange) {
		this.dish = dish; 
		this.initDishNormal = dish.getDishNormal().clone();
		this.initVal = get();
		this.min = (1.0 - fractionChange) * initVal;
		this.max = (1.0 + fractionChange) * initVal;
	}
	
	@Override
	public String toString() { 
		return "BendableDish["+dish.getName()+"]: init = " + initVal() + ", current="+get()+", min = " + min() + ", max = " + max() + ", scale = " + scale;
	}
	
	@Override
	public double get() { //distance from init pos along the given vector
		double R = dish.getRadiusOfCurv();
		boolean isFlipped = Util.dot(dish.getDishNormal(), initDishNormal) < 0; 
		
		return ((isFlipped ? -1 : 1) / R) / scale;
	}
	
	@Override
	public void set(double val) {
		boolean isFlipped = Util.dot(dish.getDishNormal(), initDishNormal) < 0; 
	
		if(val < 0){ //needs to be flipped
			if(!isFlipped){ //but isn't
				dish.setFrontMedium(initBackMedium);
				dish.setBackMedium(initFrontMedium);
				dish.setDishNormal(new double[]{ -initDishNormal[0], -initDishNormal[1], -initDishNormal[2] });
			}
			dish.setRadiusOfCurv(-1.0 / (val*scale));
		}else{ //shouldn't be flipped
			if(isFlipped){ //but is
				dish.setFrontMedium(initFrontMedium);
				dish.setBackMedium(initBackMedium);
				dish.setDishNormal(initDishNormal.clone());
			}
			dish.setRadiusOfCurv(1.0 / (val*scale));
		}	
		
	}
	
	
}
