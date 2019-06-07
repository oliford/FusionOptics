package fusionOptics.surfaces;

import net.jafama.FastMath;

import java.util.ArrayList;
import java.util.Arrays;


import java.util.List;

import fusionOptics.Util;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.RaySegment;



/** Intersection code (and constructor pre-calcs for it) taken from S.Boz's code (IPP).
 * 
 * Normal = (B-A) x (C-A) so points A,B,C
 *  are anticlockwise looking at the surface from it's 'front' side.
 * 
 * @author oliford
 */
public class Triangle extends Plane {
	
	private double A[], B[], C[], AB[], AC[];
	
	private double boundRadius;
	
	private int dirIdx2, dirIdx3;
	private double areaNormK;
	
	public Triangle(String name, double A[], double B[], double C[], Interface iface) {
		this(name, A, B, C, null, null, iface);
	}
		
	public Triangle(String name, double A[], double B[], double C[], 
			Medium frontMedium, Medium backMedium, Interface iface) {
		super(name, null, null, frontMedium, backMedium, iface);
		this.A = A; this.B = B; this.C = C;
		
		this.centre = new double[]{
				(A[0] + B[0] + C[0]) / 3,
				(A[1] + B[1] + C[1]) / 3,
				(A[2] + B[2] + C[2]) / 3,
		};
		
		AB = new double[]{ B[0] - A[0], B[1] - A[1], B[2] - A[2] };
		AC = new double[]{ C[0] - A[0], C[1] - A[1], C[2] - A[2] };
		
	    //areaNorm = find normal, not normalized - this is used to get
	    //the area later without the need to recalculate the
	    //cross product
		double areaNorm[] = Util.cross(AB, AC);
		this.normal = Util.reNorm(areaNorm.clone());
		
	    //find direction in which normal has the maximum absolute value
	    double max = FastMath.abs(normal[0]);
	    int dirIdx1 = 0;
	    if (FastMath.abs(normal[1]) > max) {
	        max = FastMath.abs(normal[1]);
	        dirIdx1 = 1;
	    }
	    if (FastMath.abs(normal[2]) > max) {
	        max = FastMath.abs(normal[2]);
	        dirIdx1 = 2;
	    }
	    
	    if (max == 0)
	        throw new IllegalArgumentException("Degenerated triangle. A,B and C are on a line.");
	    
	    areaNormK = areaNorm[dirIdx1];
	    
	    dirIdx2 = (dirIdx1 + 1) % 3;
	    dirIdx3 = (dirIdx1 + 2) % 3;
	  
		this.up = Util.reNorm(Util.minus(A, centre));
		
		double lA = FastMath.sqrt(FastMath.pow2(A[0]-centre[0]) + FastMath.pow2(A[1]-centre[2]) + FastMath.pow2(A[1]-centre[2]));
		double lB = FastMath.sqrt(FastMath.pow2(B[0]-centre[0]) + FastMath.pow2(B[1]-centre[2]) + FastMath.pow2(B[1]-centre[2]));
		double lC = FastMath.sqrt(FastMath.pow2(C[0]-centre[0]) + FastMath.pow2(C[1]-centre[2]) + FastMath.pow2(C[1]-centre[2]));
		
		this.boundRadius = FastMath.max(FastMath.max(lA, lB), lC);
		right = Util.reNorm(Util.cross(normal, up));		
	}


	@Override
	public double[] getBoundarySphereCentre() { return centre.clone(); }

	@Override
	public double getBoundarySphereRadius() { return boundRadius; }

	@Override
	public boolean findEarlierIntersection(RaySegment ray, Intersection hit) {		
		double pos[] = new double[3];
		double dist = super.calcPlaneIntersection(ray, pos); //find the place where the ray hits us (as if we were an infinite plane)
		
		//Check that it actually hit the plane, (not before, not under tolerance, not over length, and not parallel)
		if(dist < -Tracer.reHitTolerance || dist > ray.length || Double.isNaN(dist)) 
			return false;
		if(ray.startHit != null && ray.startHit.surface == this && dist < Tracer.reHitTolerance)
			return false; //rehit us within tolerance, ignore
		
		double x1 = pos[dirIdx2] - A[dirIdx2];
		double x2 = pos[dirIdx3] - A[dirIdx3];
		
		//first coord along AB
		double u = (x1 * AC[dirIdx3] - x2 * AC[dirIdx2]) / areaNormK;
		if (u < 0 || u > 1)
		    return false;
		
		//second coord along AC
		double v = (AB[dirIdx2]*x2 - AB[dirIdx3]*x1) / areaNormK;
		if (v < 0 || v > 1)
		    return false;
		
		//finaly check sum u+v
		if ( (u+v) > 1)
		    return false;
		
		hit.surface = this;
		hit.pos = pos;
		hit.normal = normal.clone();
		ray.length = dist;
		return true;
		
	}
	
	@Override
	public List<double[][]> draw() {
		ArrayList<double[][]> lines = new ArrayList<double[][]>();
		
		lines.add(new double[][]{ 
				{ A[0], B[0], C[0], A[0] },
				{ A[1], B[1], C[1], A[1] },
				{ A[2], B[2], C[2], A[2] },
			});
		
		return lines;
	}
	
	@Override
	public void shift(double[] dX) {
		super.shift(dX);
		for(int i=0;i<3;i++){
			A[i] += dX[i];
			B[i] += dX[i];
			C[i] += dX[i];
		}
			
	}
	
	@Override
	public void rotate(double[] point, double[][] matrix) {
		super.rotate(point, matrix);
		for(int i=0;i<3;i++){
			A[i] -= point[i];
			B[i] -= point[i];
			C[i] -= point[i];
		}
		
		double newAB[] = new double[3];
		double newAC[] = new double[3];
		double newA[] = new double[3];
		double newB[] = new double[3];
		double newC[] = new double[3];
		
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++){
				newAB[i] += matrix[i][j] * AB[j];
				newAC[i] += matrix[i][j] * AC[j];
				newA[i] += matrix[i][j] * A[j];
				newB[i] += matrix[i][j] * B[j];
				newC[i] += matrix[i][j] * C[j];
			}
				
		for(int i=0;i<3;i++){
			A[i] = point[i] + newA[i];
			B[i] = point[i] + newB[i];
			C[i] = point[i] + newC[i];
		}

		AB = newAB;
		AC = newAC;
	}

	@Override
	protected int planeBoundaryHashCode() {
		int r = 1;
		r = 31 * r + Arrays.hashCode(A);
		r = 31 * r + Arrays.hashCode(B);
		r = 31 * r + Arrays.hashCode(C);
		return r;
	}

	@Override
	protected boolean planeBoundaryEquals(Plane obj) {
		Triangle other = (Triangle) obj;
		return Arrays.equals(A, other.A)
				&& Arrays.equals(B, other.B)
				&& Arrays.equals(C, other.C);
	}
	

}
