package fusionOptics.surfaces;

import net.jafama.FastMath;

import java.util.ArrayList;


import java.util.List;

import fusionOptics.Util;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.RaySegment;



/** A rectangular planar surface */
public class Square extends Plane {
	
	private double height;
	private double width;
	private double boundRadius;
	
	/**@param name		Name of surface
	 * @param A	,B,C	Corner positions []{x,y,z}
	 * @param iface		Surface interface.	 */
	public Square(String name, double A[], double B[], double C[], Interface iface) {
		this(name, A, B, C, null, null, iface);
	}
		
	/**@param name		Name of surface
	 * @param A	,B,C	Corner positions []{x,y,z}
	 * @param iface		Surface interface.
	 * @param frontMedium	Medium on the +ve normal side. Normal = (B-A) x (C-B)  
	 * @param backMedium	Medium on the -ve normal side. Normal = (B-A) x (C-B) 
	 */
	public Square(String name, double A[], double B[], double C[], Medium frontMedium, Medium backMedium, Interface iface) {
		super(name, null, null, frontMedium, backMedium, iface);
		
		this.centre = new double[]{ (A[0] + C[0])/2, (A[1] + C[1])/2, (A[2] + C[2])/2 };
		
		this.up = new double[]{ B[0] - A[0], B[1] - A[1], B[2] - A[2] };
		this.height = Util.length(up);
		this.up = Util.reNorm(up);
		
		this.right = new double[]{ C[0] - B[0], C[1] - B[1], C[2] - B[2] };
		this.width = Util.length(right);
		this.right = Util.reNorm(right);
		
		this.normal = Util.cross(up, right);
		this.boundRadius = FastMath.sqrt(width*width + height*height)/2;		
	}
	
	/**@param centre	Position of rectangle centre
	 * @param normal	Unit normal to surface
	 * @param upDir		Direction of up (height) - should be perp to normal.
	 * @param height	Full size in 'up' direction.
	 * @param width		Full size in 'righ' direction.
	 * @param iface		Surface interface.	
	 */
	public Square(String name, double centre[], double normal[], double upDir[], double height, double width, Interface iface) {
		this(name, centre, normal, upDir, height, width, null, null, iface);
	}
		
	/**@param centre	Position of rectangle centre
	 * @param normal	Unit normal to surface
	 * @param upDir		Direction of up (height) - should be perp to normal.
	 * @param height	Full size in 'up' direction.
	 * @param width		Full size in 'righ' direction.
	 * @param frontMedium	Medium on the +ve normal side. Normal = (B-A) x (C-B)  
	 * @param backMedium	Medium on the -ve normal side. Normal = (B-A) x (C-B) 
	 * @param iface		Surface interface.	
	 */
	public Square(String name, double centre[], double normal[], double upDir[], double height, double width, 
			Medium frontMedium, Medium backMedium, Interface iface) {
		super(name, centre, normal, frontMedium, backMedium, iface);
		this.up = upDir;
		this.height = height;
		this.width = width;
		this.boundRadius = FastMath.sqrt(width*width + height*height)/2;
		right = Util.cross(normal,upDir);		
	}


	@Override
	public double[] getBoundarySphereCentre() { return centre.clone(); }

	@Override
	public double getBoundarySphereRadius() { return boundRadius; }

	@Override
	public boolean findEarlierIntersection(RaySegment ray, Intersection hit) {		
		double pos[] = new double[3];
		double dist = super.calcPlaneIntersection(ray, pos); //find the place where the ray hits us (as if we were an infinite plane)
		
		//Check that it actually hit the plane, (before us, not over length, and not parallel)
		if(dist < -Tracer.reHitTolerance || dist > ray.length || Double.isNaN(dist)) 
			return false;
		if(ray.startHit != null && ray.startHit.surface == this && dist < Tracer.reHitTolerance)
			return false; //rehit us within tolerance, ignore
		
		double paraInSquare = 0, perpInSquare = 0;
		
		for(int i=0;i<3;i++){
			paraInSquare += (pos[i] - centre[i]) * up[i];
			perpInSquare += (pos[i] - centre[i]) * right[i];			
		}
		
		
		if( Math.abs(paraInSquare) <= (height/2.0) && Math.abs(perpInSquare) <= (width/2.0)){
			hit.surface = this;
			hit.pos = pos;
			hit.normal = normal.clone();
			ray.length = dist;
			return true;
		}else
			return false;
	}
	
	@Override
	public List<double[][]> draw() {
		ArrayList<double[][]> lines = new ArrayList<double[][]>();
		
		//edge
		double eLine[][] = new double[3][5];
		double lLine[][] = new double[3][2];
		double wLine[][] = new double[3][2];
		for(int k=0; k < 3; k++){
			eLine[k][0] = centre[k] - height/2 * up[k] - width/2 * right[k];
			eLine[k][1] = centre[k] - height/2 * up[k] + width/2 * right[k];
			eLine[k][2] = centre[k] + height/2 * up[k] + width/2 * right[k];
			eLine[k][3] = centre[k] + height/2 * up[k] - width/2 * right[k];
			eLine[k][4] = centre[k] - height/2 * up[k] - width/2 * right[k];
			
			lLine[k][0] = centre[k] - height/2 * up[k];
			lLine[k][1] = centre[k] + height/2 * up[k];
			
			wLine[k][0] = centre[k] - width/2 * right[k];
			wLine[k][1] = centre[k] + width/2 * right[k];
		}
		lines.add(eLine);
		lines.add(lLine);
		lines.add(wLine);
		
		return lines;
	}
	
	public void setHeight(double height){ this.height = height; boundRadius = FastMath.sqrt(width*width + height*height)/2; }
	public void setWidth(double width){ this.width = width; boundRadius = FastMath.sqrt(width*width + height*height)/2; }

	@Override
	protected int planeBoundaryHashCode() {
		long t; int r = 1;
		t = Double.doubleToLongBits(width);		r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(height);	r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}

	@Override
	protected boolean planeBoundaryEquals(Plane obj) {
		return (this == obj) || ( super.equals(obj)
				&& Double.doubleToLongBits(width) == Double.doubleToLongBits(((Square)obj).width)
				&& Double.doubleToLongBits(height) == Double.doubleToLongBits(((Square)obj).height) );				
	}

	public double getWidth() { return width; }
	public double getHeight() { return height; }

}
