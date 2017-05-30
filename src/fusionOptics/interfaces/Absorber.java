package fusionOptics.interfaces;

import java.util.List;

import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.Surface;


public final class Absorber implements Interface {

	private final static Absorber ideal = new Absorber();
	
	private Absorber(){ } //protect - because there's no need to keep creating this

	@Override
	public final void calcIntersection(Intersection newHit, double minIntensity) {
		//Nothing to do
	}

	/** Returns a perfect absorber interface (global single instance) */
	public final static Absorber ideal(){ return ideal; }

	@Override
	public final void checkCompatibility(Surface surface) {
		//everything thing is OK
	}
	
	//all absorbers are alike
	@Override
	public int hashCode() { return getClass().getName().hashCode(); }
	@Override
	public boolean equals(Object obj) { return (obj != null) && obj instanceof Absorber; }
}
