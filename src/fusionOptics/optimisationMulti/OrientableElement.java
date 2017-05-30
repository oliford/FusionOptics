package fusionOptics.optimisationMulti;

import java.util.ArrayList;
import java.util.List;

import algorithmrepository.Algorithms;
import jafama.FastMath;
import fusionOptics.Util;
import fusionOptics.surfaces.Plane;
import fusionOptics.types.Element;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;

/** Plane rotations.
 * 
 *  Pan sets the rotation of the plane around the original up, regardless of the current tilt.
 *  Tilt sets the angle around the CURRENT 'right', regardless of the current pan.
 *  
 */
public class OrientableElement extends Parameter {
	
	/** The element to move */
	private Element element;
	
	/** We need a plane to use for directions */
	private Plane plane;
		
	private double[] rotationCentre;
	private double[] initNormal, initUp, initRight;
	
	private int param;
	private static final int PARAM_TILT = 1;
	private static final int PARAM_PAN = 2;
	
	public static List<Parameter> panAndTilt(Element element, double min, double max) {
		Plane plane = findPlane(element);
		List<Parameter> params = new ArrayList<Parameter>(2);
		
		params.add(new OrientableElement(element, plane, PARAM_PAN, min, max, Math.PI / 180.0));
		params.add(new OrientableElement(element, plane, PARAM_TILT, min, max, Math.PI / 180.0));
		
		return params;
	}
	
	/** Search the whole tree of elements to find any Plane */
	private static Plane findPlane(Element element){
		if(element instanceof Plane){
			return (Plane)element;
		}
		
		if(element instanceof Optic){
			for(Surface surface : ((Optic)element).getSurfaces()){
				if(surface instanceof Plane)
					return (Plane)surface;
			}
			
			for(Optic optic : ((Optic)element).getSubOptics()){
				Plane plane = findPlane(optic);
				if(plane != null)
					return plane;
			}
		}
		
		return null;
	}
	
	public OrientableElement(Element element, int param, double min, double max) {
		this(element, findPlane(element), param, min, max, 1.0);
	}
	
	public OrientableElement(Element element, int param, double min, double max, double scale) {
		this(element, findPlane(element), param, min, max, scale);
	}
	
	public OrientableElement(Element element, Plane orientationPlane, int param, double min, double max, double scale) {
		super(min, max, 0); //parameter is relative to the initial bounding sphere centre, so starts at 0
		
		if(orientationPlane == null)
			throw new IllegalArgumentException("Object "+element.getName()+" has not planes in it to use to know it's orientation"); 
		
		this.element = element;
		this.plane = orientationPlane;
		this.rotationCentre = element.getBoundarySphereCentre();
		this.initNormal = plane.getNormal();
		this.initUp = plane.getUp();
		this.initRight = plane.getRight();
		this.param = param;
		this.scale = scale;
	}
	
	@Override
	public String toString() { 
		double n[] = plane.getNormal();
		return "OrientableElement["+element.getName()+", "+((param == PARAM_TILT) ? "Tilt" : "Pan")+"]: init = " + initVal() + ", current="+get()+", min = " + min() + ", max = " + max() + ", scale = " + scale + ", norm={ "+n[0]+","+n[1]+","+n[2]+"}";
	}
	
	@Override
	/** Get position relative to zeroPos[], along vector */
	public double get() { //distance from init pos along the given vector
		
		double n[] = plane.getNormal();
		double nR = Util.dot(n, initRight);
		double nN = Util.dot(n, initNormal);
		double nU = Util.dot(n, initUp);		
		double tilt = FastMath.atan2(nU, FastMath.sqrt(nN*nN + nR*nR));
		double pan = FastMath.atan2(nR, nN);
		return ((param == PARAM_TILT) ? tilt : pan) / scale;
	}
	
	@Override
	/** Move to new position relative to zeroPos[] */
	public void set(double val) {
		//get current pan/tilt
		double n[] = plane.getNormal();
		double currentRight[] = plane.getRight();
		double nR = Util.dot(n, initRight);
		double nN = Util.dot(n, initNormal);
		double nU = Util.dot(n, initUp);		
		double currentTilt = FastMath.atan2(nU, FastMath.sqrt(nN*nN + nR*nR));
		double currentPan = FastMath.atan2(nR, nN);
		
		
		if(param == PARAM_TILT){
			double newTilt = val * scale;
			element.rotate(rotationCentre, Algorithms.rotationMatrix(currentRight, newTilt - currentTilt));
			
		}else{ 
			double newPan = val * scale;
			element.rotate(rotationCentre, Algorithms.rotationMatrix(initUp, -(newPan - currentPan)) );
			
		}
	}
	
	
	
}
