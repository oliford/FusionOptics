package fusionOptics.types;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import algorithmrepository.Algorithms;

/** Describes an intersection of a ray with a surface */
public class Intersection {
	
	/** Surface which was intersected */
	public Surface surface;
	
	/** Intersection position */
	public double pos[];
	
	/** Surface normal at intersection point (unit vec) */
	public double normal[];
	
	/** Medium that the ray is leaving (or reflecting into) */
	public Medium incidentMedium;
	
	/** Medium that the transmitted ray (if any) is entering. */
	public Medium transmissionMedium;
	
	/** Incoming ray that hit the surface */
	public RaySegment incidentRay;
	
	/** Various possible reflected or transmitted rays */
	public RaySegment reflectedOrdinary;
	public RaySegment reflectedExtraordinary;
	public RaySegment transmittedOrdinary;
	public RaySegment transmittedExtraordinary;
	
	@Override
	/** Copys arrays but not references to segements */
	protected Intersection clone()  {
		Intersection newHit = new Intersection();
				
		newHit.pos = pos.clone();
		newHit.normal = normal.clone();
		
		newHit.surface = surface;
		newHit.incidentMedium = incidentMedium;
		newHit.transmissionMedium = transmissionMedium;
		newHit.incidentRay = incidentRay;
		
		newHit.reflectedOrdinary = reflectedOrdinary;
		newHit.reflectedExtraordinary = reflectedExtraordinary;
		newHit.transmittedOrdinary = reflectedOrdinary;
		newHit.transmittedExtraordinary = reflectedExtraordinary;
				
		return newHit;
	}
	
	public final boolean isStart(){ return incidentRay == null; }

	public final boolean isEnd(){
		return reflectedOrdinary == null && reflectedExtraordinary == null && 
				transmittedOrdinary == null && transmittedExtraordinary == null;
	}
	
	public final int getNOutgoingRays(){
		return (reflectedOrdinary != null ? 1 : 0) +
				(reflectedExtraordinary != null ? 1 : 0) +
				(transmittedOrdinary != null ? 1 : 0) +
				(transmittedExtraordinary != null ? 1 : 0);
	}
	
	
	public final RaySegment[] getSortedOutgoingRaysArray() {
		boolean TO = transmittedOrdinary != null;
		boolean TE = transmittedExtraordinary != null;
		boolean RO = reflectedOrdinary != null;
		boolean RE = reflectedExtraordinary != null;
		
		//common simple cases
		if(TO && !TE && !RO && !RE)
			return new RaySegment[]{ transmittedOrdinary };
		else if(!TO && !TE && RO && !RE)
			return new RaySegment[]{ reflectedOrdinary };
		else if(TO && !TE && RO && !RE){
			if(reflectedOrdinary.startIntensity() > transmittedOrdinary.startIntensity())
				return new RaySegment[]{ reflectedOrdinary, transmittedOrdinary};
			else
				return new RaySegment[]{ transmittedOrdinary, reflectedOrdinary};
		}
			
		//othwerwise, do the slow build and sort
		RaySegment[] rays = new RaySegment[]{ 
				reflectedOrdinary,
				reflectedExtraordinary,
				transmittedOrdinary,
				transmittedExtraordinary
		};
		
		double I[] = new double[]{
				(RO ? reflectedOrdinary.startIntensity() : 0),
				(RE ? reflectedExtraordinary.startIntensity() : 0),
				(TO ? transmittedOrdinary.startIntensity() : 0),
				(TE ? transmittedExtraordinary.startIntensity() : 0), 
		};
		int idx[] = new int[]{ 0, 1, 2, 3 };
		
		Algorithms.quicksort(I, idx);
		
		return new RaySegment[]{ 
				rays[idx[3]],
				rays[idx[2]],
				rays[idx[1]],
				rays[idx[0]],
		};
		
		
	}
	
	
	/** Returns a list of outgoing rays from this intersection,
	 * sorted by total intensity, with the strongest first  */
	public final List<RaySegment> getSortedOutgoingRays() {
		ArrayList<RaySegment> rays = new ArrayList<RaySegment>(getNOutgoingRays());
		
		//get all new rays
		if(reflectedOrdinary != null)
			rays.add(reflectedOrdinary);
		if(reflectedExtraordinary != null)
			rays.add(reflectedExtraordinary);
		if(transmittedOrdinary != null)
			rays.add(transmittedOrdinary);
		if(transmittedExtraordinary != null)
			rays.add(transmittedExtraordinary);
		
		//if there are 2 or more, sort by intensity (strongest first)
		if(rays.size() > 1){
			Collections.sort(rays, new Comparator<RaySegment>() {
				@Override
				public int compare(RaySegment o1, RaySegment o2) {
					return Double.compare(o2.startIntensity(), o1.startIntensity());						
				}
			});
		}
		
		return rays;
	}
	
	@Override
	public String toString() {
		return "S="+((surface==null)?"null":surface.getName()) + ",O=" +
						((transmittedOrdinary != null) ? ",T" : "") +
						((transmittedExtraordinary != null) ? ",TE" : "") +
						((reflectedOrdinary != null) ? ",R" : "") +
						((reflectedExtraordinary != null) ? ",RE" : "");
	}
}
