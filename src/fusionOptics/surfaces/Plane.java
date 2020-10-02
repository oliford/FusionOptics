package fusionOptics.surfaces;

import net.jafama.FastMath;

import java.util.Arrays;
import java.util.List;

import fusionOptics.Util;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;
import algorithmrepository.exceptions.NotImplementedException;




/** An infinite plane, this is also the superclass of elements which are finite planes, like Disc and Square */
public abstract class Plane extends Surface {
	
	protected double centre[];
	protected double normal[]; //normal vector
	protected double up[]; //two vectors in the plane
	protected double right[];
	
	public Plane(String name, double centre[], double normal[], Medium frontMedium, Medium backMedium, Interface iface){
		super(name, frontMedium, backMedium, iface);
		this.centre = (centre == null) ? null : centre.clone();
		this.normal = (normal == null) ? null : normal.clone();
	}
	
	protected void initArbitraryPerps(){
		//we need 2 perp vecrors to normal
		//Ideally, we'd like 'up' to be as 'up' as possible, so start by crossing Z with normal
		this.right = Util.cross(normal, new double[]{0,0,1} );
		if(right[0]*right[0] + right[1]*right[1] + right[2]*right[2] == 0) //oh well, normal is up
			right = Util.cross(normal, new double[]{0,1,0} );
		
		Util.reNorm(right);
		
		this.up = Util.cross(right, normal);
	}
	
	public void shift(double[] dX) {
		for(int i=0;i<3;i++)
			centre[i] += dX[i];
	}
	
	public void rotate(double point[], double matrix[][]){
		for(int i=0;i<3;i++)
			centre[i] -= point[i];
			
		double newCentre[] = new double[3];
		double newNormal[] = new double[3];
		double newUp[] = new double[3];
		double newRight[] = new double[3];
		
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++){
				newCentre[i] += matrix[i][j] * centre[j];
				newNormal[i] += matrix[i][j] * normal[j];
				newUp[i] += matrix[i][j] * up[j];
				newRight[i] += matrix[i][j] * right[j];
			}
				
		for(int i=0;i<3;i++)
			centre[i] = point[i] + newCentre[i];

		normal = Util.reNorm(newNormal);
		up = Util.reNorm(newUp);
		right = Util.reNorm(newRight);
	}
	
	/** Calculate the intersection of the given ray with the plane
	 * @param position vector to fill in
	 * 
	 * @return distance along ray
	 */
	public final double calcPlaneIntersection(RaySegment ray, double pos[]){
		
		double uDotN = Util.dot(ray.dir, normal);
		
		if(uDotN == 0)   // if // to plane, it never hits
			return Double.NaN;	
				
		double dist = (centre[0] - ray.startPos[0]) * normal[0] + 
						(centre[1] - ray.startPos[1]) * normal[1] + 
						(centre[2] - ray.startPos[2]) * normal[2];
							
		dist /= uDotN;
						
		pos[0] = ray.startPos[0] + dist * ray.dir[0];
		pos[1] = ray.startPos[1] + dist * ray.dir[1];
		pos[2] = ray.startPos[2] + dist * ray.dir[2];
		
		return dist;
	}

	
	/** Returns the position in ('up','right')  space of this plane, of a certesian XYZ position 
	 * @deprecated This was back-to-front, 2D coords should be (R,U) not (U,R) */
	public final double[] posXYZToPlaneUR(double posXYZ[]){
		return new double[]{
				(posXYZ[0] - centre[0]) * up[0] +
				(posXYZ[1] - centre[1]) * up[1] +
				(posXYZ[2] - centre[2]) * up[2],
				
				(posXYZ[0] - centre[0]) * right[0] +
				(posXYZ[1] - centre[1]) * right[1] +
				(posXYZ[2] - centre[2]) * right[2],
			};
	}

	/** Returns the position in ('right', 'up') space of this plane, of a certesian XYZ position */
	public final double[] posXYZToPlaneRU(double posXYZ[]){
		return new double[]{
				(posXYZ[0] - centre[0]) * right[0] +
				(posXYZ[1] - centre[1]) * right[1] +
				(posXYZ[2] - centre[2]) * right[2],
				
				(posXYZ[0] - centre[0]) * up[0] +
				(posXYZ[1] - centre[1]) * up[1] +
				(posXYZ[2] - centre[2]) * up[2]				
			};
	}
	
	/** Returns the position in 'right' 'up' space of this plane, of a certesian XYZ position */
	public final double[] planeRUToPosXYZ(double posRU[]){
		return new double[]{
				centre[0] + posRU[0]*right[0] + posRU[1]*up[0],
				centre[1] + posRU[0]*right[1] + posRU[1]*up[1],
				centre[2] + posRU[0]*right[2] + posRU[1]*up[2],
			};
	}
	
	public double[] getCentre() { return centre.clone(); }
	
	public double[] getNormal() { return normal.clone(); }

	public double[] getUp() { return up; }
	
	public double[] getRight() { return right; }
	
	public void setUp(double up[]){
		this.right = Util.reNorm(Util.cross(normal, up));
		this.up = Util.cross(right, normal);
	};

	/** Modifies the plane's normal. Up vector is set to the closest to what is was */
	public void setNormal(double[] n) {
		this.normal = Util.reNorm(n);
		
		//for now, set up and right to the nearest we can get to what they were
		this.right = Util.reNorm(Util.cross(normal, up));
		this.up = Util.cross(right, normal);
	}
	
	public void setCentre(double[] centre) { this.centre = centre; }

	@Override
	public int surfaceGeometryHashCode() {
			return ( ( ( 31 + planeBoundaryHashCode()
					     ) * 31 + Arrays.hashCode(centre)
					   ) * 31 + Arrays.hashCode(normal)
					 ) * 31 + Arrays.hashCode(up);
	}

	@Override
	public boolean surfaceGeometryEquals(Surface obj) {
		Plane other = (Plane) obj;
		return planeBoundaryEquals(other)
				&& Arrays.equals(centre, other.centre)
				&& Arrays.equals(normal, other.normal)
				&& Arrays.equals(up, other.up);
	}

	/** HashCode of surface plane boundary implementation */
	protected abstract int planeBoundaryHashCode();
	
	/** .equals() methods for plane boundary implementation */
	protected abstract boolean planeBoundaryEquals(Plane other);
	
	/** @return The angle of the 'up' of this plane from the given vector projected onto the plane */
	public double getOrientationAngle(double relativeTo[]){
		//project the  'relativeTo' into the plane
		double u = Util.dot(relativeTo, up);
		double r = Util.dot(relativeTo, right);
		
		//so the angle of our up, relative to that, is -ve that
		return -FastMath.atan2(r, u);
	}
}
