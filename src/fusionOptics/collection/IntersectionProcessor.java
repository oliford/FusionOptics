package fusionOptics.collection;

import fusionOptics.types.Intersection;


/** Callback for anything which handles hit points with a specific element */ 
public interface IntersectionProcessor {

	public void nextIntersection(Intersection hit);
}
