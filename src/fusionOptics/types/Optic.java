package fusionOptics.types;

import net.jafama.FastMath;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;


/** Describes a collection of surfaces, media and other optics.
 *  
 * Something with materials and surfaces.
 * 
 * @author oliford
 */
public class Optic implements Element {
	/** Surfaces of this optic (not including sub optics). Inheritors should set this. */
	private List<Surface> surfaces;
		
	/** Sub optics list. Surfaces of sub optics should not be in this optic's surfaces list */
	private List<Optic> subOptics;
			
	public boolean enableBoundingCheck = true;
	
	/** Identity of this optic */
	protected String name;
	
	public Optic(String name) {
		this.name = name;
		this.surfaces =  new LinkedList<Surface>();
		this.subOptics = new LinkedList<Optic>();
	}
		
	/** Instantiate an optic from an array of surfaces and/or optics */ 
	public Optic(String name, Element[] elements) {
		this.name = name;
		this.surfaces =  new LinkedList<Surface>();
		this.subOptics = new LinkedList<Optic>();
		for(Element element : elements) {
			addElement(element);
		}
	}
		
	/** Instantiate an optic from an array of surfaces and/or optics */ 
	public Optic(String name, Iterable<Element> elements) {
		this.name = name;
		this.surfaces =  new LinkedList<Surface>();
		this.subOptics = new LinkedList<Optic>();
		for(Element element : elements) {
			addElement(element);
		}
	}
		
	private double boundSphereCentre[] = null;
	private double boundSphereRadius = Double.NaN;
	
	private double geometricCentre[] = null;
	
	public String getName(){ return name; }
	
	/** Returns a list of sub optics of this optic */
	public List<Optic> getSubOptics() { return subOptics; }

	/** Returns a list of surfaces that make up this object */ 
	public List<Surface> getSurfaces() { return surfaces; }
	
	public List<Surface> getSurfacesAll() {
		List<Surface> allSurfaces = new ArrayList<Surface>(this.surfaces);
		for(Optic optic : subOptics)
			allSurfaces.addAll(optic.getSurfacesAll());
		
		return allSurfaces;
	}
	
	
	/** Returns a list of media of which this optic is made */ 
	//public List<Medium> getMedia() { return media; }
		
	/** Returns the approximate centre of the object (used for drawing and bounding sphere calc 
	 * This base method just takes maximum of all surfaces.
	 * Implementors should override with something faster if possible.	 */
	public double[] getBoundarySphereCentre() {
		if(boundSphereCentre != null)
			return boundSphereCentre;
		boundSphereCentre = new double[3];
		
		int n=0;
		for(Surface surface : surfaces) {
			double c[] = surface.getBoundarySphereCentre();
			boundSphereCentre[0] += c[0];
			boundSphereCentre[1] += c[1];
			boundSphereCentre[2] += c[2];
			n++;
		}
		for(Optic subOptic : subOptics) {
			double c[] = subOptic.getBoundarySphereCentre();
			boundSphereCentre[0] += c[0];
			boundSphereCentre[1] += c[1];
			boundSphereCentre[2] += c[2];
			n++;
		}
		if(n == 0)
			throw new RuntimeException("Trying to calculate bounding sphere for Optic '"+getName()+"' with no subelements");
		boundSphereCentre[0] /= n;
		boundSphereCentre[1] /= n;
		boundSphereCentre[2] /= n;
		return boundSphereCentre;
	}
	
	/** Returns the approximate centre of the object (as an arbitrary definition)
	 * This base method just takes the average of the centres of all sub-optics and surfaces
	 * and therefore doesn't necessarily mean anything sensible.	 */
	public double[] getCentre() {
		if(geometricCentre != null)
			return geometricCentre;
		geometricCentre = new double[3];
		
		int n=0;
		for(Surface surface : surfaces) {
			double c[] = surface.getCentre();
			geometricCentre[0] += c[0];
			geometricCentre[1] += c[1];
			geometricCentre[2] += c[2];
			n++;
		}
		for(Optic subOptic : subOptics) {
			double c[] = subOptic.getCentre();
			geometricCentre[0] += c[0];
			geometricCentre[1] += c[1];
			geometricCentre[2] += c[2];
			n++;
		}
		geometricCentre[0] /= n;
		geometricCentre[1] /= n;
		geometricCentre[2] /= n;
		return geometricCentre;
	}
	
	/** Returns the radius of the bounding sphere of this optic 
	 * This base method just takes maximum of all surfaces. 
	 * Implementors should override with something faster if possible.	 */
	public double getBoundarySphereRadius(){
		if(!Double.isNaN(boundSphereRadius))
			return boundSphereRadius;
		if(boundSphereCentre == null)
			boundSphereCentre = getBoundarySphereCentre(); //set, because it might be overridden
		
		boundSphereRadius = 0;
		for(Surface surface : surfaces) {
			double c[] = surface.getBoundarySphereCentre();
			double r = surface.getBoundarySphereRadius();
			double dx = (c[0]-boundSphereCentre[0]);
			double dy = (c[1]-boundSphereCentre[1]);
			double dz = (c[2]-boundSphereCentre[2]);
			double d = FastMath.sqrt(dx*dx+dy*dy+dz*dz);
			if( (d+r) > boundSphereRadius )
				boundSphereRadius = d+r;
		}
		for(Optic subOptic : subOptics) {
			double c[] = subOptic.getBoundarySphereCentre();
			double r = subOptic.getBoundarySphereRadius();
			double dx = (c[0]-boundSphereCentre[0]);
			double dy = (c[1]-boundSphereCentre[1]);
			double dz = (c[2]-boundSphereCentre[2]);
			double d = FastMath.sqrt(dx*dx+dy*dy+dz*dz);
			if( (d+r) > boundSphereRadius )
				boundSphereRadius = d+r;
		}
		return boundSphereRadius;
	}
	
	/** Test if the given ray enters the optic's bounding sphere */
	public boolean testBoundingSphere(RaySegment ray) {
		
		double c[] = getBoundarySphereCentre();
		double r = 1*getBoundarySphereRadius(); //sphere radius
	
		double ac0 = c[0] - ray.startPos[0];
		double ac1 = c[1] - ray.startPos[1];
		double ac2 = c[2] - ray.startPos[2];
	    //double[] abxac = Util.cross(ray.dir, ac);
		
		/* Some mostly unhelpful optimisation, only gives a per % at best
		double r2 = r*r;
		double distToCentre2 = ac0*ac0 + ac1*ac1 + ac2*ac2;
		if(distToCentre2 > r2){
			double distToCentre = FastMath.sqrt(distToCentre2);
			if((distToCentre - r) > ray.length){ 
				return false; //ray will never reach the sphere
			}
		}*/
				
	    
	  //shortest distance perp to line, to centre of bounding sphere
	    double perpDist2 = (ray.dir[1] * ac2 - ray.dir[2] * ac1) * (ray.dir[1] * ac2 - ray.dir[2] * ac1) +
					    (-ray.dir[0] * ac2 + ray.dir[2] * ac0) * (-ray.dir[0] * ac2 + ray.dir[2] * ac0) +
					    (ray.dir[0] * ac1 - ray.dir[1] * ac0) * (ray.dir[0] * ac1 - ray.dir[1] * ac0);

	    
	    if(perpDist2 > r*r) //if it's over, line never hits
	    	return false;
	    
	    //dist along ray of the point nearest the sphere centre
	    //double l = Util.dot(ac, ray.dir);
	    double l = ac0*ray.dir[0] + ac1*ray.dir[1] + ac2*ray.dir[2];
	    
	    //dist forward and back of that point, of the contact with the sphere
	    double dl = FastMath.sqrt(r*r - perpDist2);
	    
	    double l0 = l - dl;
	    double l1 = l + dl;
	    
	    return (l0 > 0 && l0 < ray.length) || //ray hits entering edge
	    	(l1 > 0 && l1 < ray.length) || //ray hist entering edge
	    	(l0 < 0 && l1 > ray.length); // ray is entirely inside sphere
	    
	    
	}

	/** Construct arrays of the sub-optics, surfaces and media references
	 * with each medium in only once (so rotations don't mess them up) */
	public List<Medium> getUniqueMediaList(){
		LinkedList<Medium> media = new LinkedList<Medium>();
		for(Surface s : getSurfaces()){
			List<Medium> mList = s.getMedia();
			if(mList != null){
				for(Medium m : mList) {
					if(m != null && !media.contains(m))
						media.add(m);
				}
			}
		}
		return media;
	}


	/** Tests if the given ray intersects any surface in this optic, and if so, fills some information 
	 * in to the given intersection object:  surface, pos[], normal[]
	 * This includes testing if it is within the current ray segment length, 
	 *   so calling this for multiple surfaces will leave it with the first intersection.
	 *    
	 * Implementors are not obliged to set the normal in any particular direction relative to the ray.
	 */
	public boolean findEarlierIntersection(RaySegment ray, Intersection hit) {
		if(enableBoundingCheck && !testBoundingSphere(ray))
			return false;
		
		boolean earlierFound = false;
		for(Optic optic : subOptics) {
			earlierFound |= optic.findEarlierIntersection(ray, hit);
		}
		for(Surface surface : surfaces) {
			earlierFound |= surface.findEarlierIntersection(ray, hit);
		}
		return earlierFound;
	}
	
	/** Applies a translation to this optic */
	public void shift(double dX[]) {
		for(Optic optic : subOptics)
			optic.shift(dX);
		for(Surface surface : surfaces)
			surface.shift(dX);
		//media don't need translations

		boundSphereCentre = null;
		boundSphereRadius = Double.NaN;
	}
	
	/** Applies a rotation to this optic */
	public void rotate(double point[], double matrix[][]) {
		for(Optic o : subOptics)
			o.rotate(point, matrix);

		for(Surface s : surfaces)
			s.rotate(point, matrix);
		//It is up to us to rotate the media, not the surfaces
		//so that we know that we do each only once
		for(Medium m : getUniqueMediaList())
			m.rotate(point, matrix);
		
		boundSphereCentre = null;
		boundSphereRadius = Double.NaN;
	}
	
	/** Returns lines (lists of, lists of double[x/y/z]s), describing the surfaces */ 
	public List<List<double[][]>> draw(){
		ArrayList<List<double[][]>> lineGroups = new ArrayList<List<double[][]>>();
		for(Optic optic : subOptics) {
			lineGroups.addAll(optic.draw());
		}
		for(Surface surface : surfaces) {
			lineGroups.add(surface.draw());
		}
		return lineGroups;
	}

	/** Adds a surface or element to the optic */
	public void addElement(Element element) {
		if( element instanceof Optic ){
			subOptics.add((Optic) element);
		}else if( element instanceof Surface ){
			surfaces.add((Surface) element);
			((Surface) element).setOptic(this);
		}else
			throw new IllegalArgumentException("Optics can only be made up of surfaces or other optics. Element '"+element+"' is neither.");

		boundSphereCentre = null;
		boundSphereRadius = Double.NaN;
	}
	
	public void removeElement(Element element){
		if( element instanceof Optic ){
			subOptics.remove(element);
		}else if( element instanceof Surface ){
			surfaces.remove(element);
		}else
			throw new IllegalArgumentException("Element '"+element+"' wasn't a surface or an optic.");
		
		boundSphereCentre = null;
		boundSphereRadius = Double.NaN;
	}

	/**  Returns a tree of List of List of ... of double[][] that describe the optic element
	 * drawings.  
	 * This is opposed to OpticCollection.draw() which will return a list of lines in each surface.
	 * (e.g. the tree to only one depth) 
	 * @return A list of objects, that can be also lists, or somewhere down the tree, a line: double[vertexIdx][x/y/z] */
	public List<Object> drawTreeGroups(){
		ArrayList<Object> drawTree = new ArrayList<Object>();  
		for(Optic optic : subOptics) {
			drawTree.add(optic.drawTreeGroups());
		}
		for(Surface surface : surfaces) {
			drawTree.add(surface.draw());
		}
		return drawTree;
	}

	public void setName(String name) { this.name = name; }

	@Override
	public int hashCode() {
		int r = 1;
		r = 31 * r + name.hashCode();
		r = 31 * r + subOptics.hashCode();
		r = 31 * r + surfaces.hashCode();
		return r;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		
		Optic other = (Optic) obj;
		return name.equals(other.name) 
				&& subOptics.equals(other.subOptics) 
				&& surfaces.equals(other.surfaces);
	}

	@Override
	public void setApproxDrawQuality(int approxDrawQuality) {
		for(Optic subOptic : subOptics)
			subOptic.setApproxDrawQuality(approxDrawQuality);
		
		for(Surface surface : surfaces)
			surface.setApproxDrawQuality(approxDrawQuality);							
	}
	
	
}

