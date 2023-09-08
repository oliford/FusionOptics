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
import fusionOptics.types.Surface;

import algorithmrepository.Algorithms;



public class TruncatedCylinder extends Surface {
	
	private double centre[];
	private double radius;
	private double boundRadius;
	private double length;
	private double axis[];
	private double up[];
	private double right[];
	/** limits relative to cylinder axis in up and right directions */
	private double rMin, rMax, uMin, uMax; 
	
	/** Drawing properties */
	private int nPointsPerSection = -1;
	private int nSections = -1;
	
	public TruncatedCylinder(String name, double centre[], double axis[], double up[], double radius, double length, 
			double height, double width, Interface iface) {
		this(name, centre, axis, up, radius, length, height, width, null, null, iface);
	}

	/** right half of cylinder of given height symmetric around horizonal midplane */
	public TruncatedCylinder(String name, double centre[], double axis[], double up[], double radius, double length, double height, 
			Medium innerMedium, Medium outerMedium, Interface iface) {
		this(name, centre, axis, up, radius, length,
				0, Double.POSITIVE_INFINITY,
				-height/2,
				height/2,
				innerMedium, outerMedium, iface);
	}
	
	public TruncatedCylinder(String name, double centre[], double axis[], double up[], double radius, double length, 
			double height, double width,
			Medium innerMedium, Medium outerMedium, Interface iface) {
		this(name, centre, axis, up, radius, length,
				-width/2, width/2,
				-height/2, height/2,
				innerMedium, outerMedium, iface);
	}

	public TruncatedCylinder(String name, double centre[], double axis[], double up[], double radius, double length, 
			double rMin, double rMax, double uMin, double uMax, 
			Medium innerMedium, Medium outerMedium, Interface iface) {
		super(name, outerMedium, innerMedium, iface);
		this.centre =  centre;
		this.axis = axis;
		this.length = length;
		this.radius = radius;
		this.up = up;
		this.right = Util.reNorm(Util.cross(axis, up));
		this.rMin = rMin;
		this.rMax = rMax;
		this.uMin = uMin;
		this.uMax = uMax;
				
		this.boundRadius = FastMath.sqrt(radius*radius + length*length/4);
	}

	@Override
	public double[] getBoundarySphereCentre() { return centre.clone(); }

	@Override
	public double getBoundarySphereRadius() { return boundRadius; }	
	
	@Override
	public boolean findEarlierIntersection(RaySegment ray, Intersection hit) {
		
		double t[] = Algorithms.cylinderLineIntersection(ray.startPos, ray.dir, centre, axis, radius*radius);
		
		if(t.length <= 0)
			return false;
		
		
		//calc hit positions
		double P0[] = new double[]{
				ray.startPos[0] + t[0] * ray.dir[0],
				ray.startPos[1] + t[0] * ray.dir[1],
				ray.startPos[2] + t[0] * ray.dir[2],
		};
		double P1[] = new double[]{
				ray.startPos[0] + t[1] * ray.dir[0],
				ray.startPos[1] + t[1] * ray.dir[1],
				ray.startPos[2] + t[1] * ray.dir[2],
		};
		
		//calc dist along cylinder axis of each point: (P-C).N
		double d0 = (P0[0] - centre[0]) * axis[0] + (P0[1] - centre[1]) * axis[1] + (P0[2] - centre[2]) * axis[2]; 
		double d1 = (P1[0] - centre[0]) * axis[0] + (P1[1] - centre[1]) * axis[1] + (P1[2] - centre[2]) * axis[2];
		
		double u0 = (P0[0] - centre[0]) * up[0] + (P0[1] - centre[1]) * up[1] + (P0[2] - centre[2]) * up[2]; 
		double u1 = (P1[0] - centre[0]) * up[0] + (P1[1] - centre[1]) * up[1] + (P1[2] - centre[2]) * up[2];
		
		double r0 = (P0[0] - centre[0]) * right[0] + (P0[1] - centre[1]) * right[1] + (P0[2] - centre[2]) * right[2]; 
		double r1 = (P1[0] - centre[0]) * right[0] + (P1[1] - centre[1]) * right[1] + (P1[2] - centre[2]) * right[2];
		
		boolean lastSurfWasUs = (ray.startHit != null && ray.startHit.surface == this);
		double reHitTolerance =  lastSurfWasUs ? Tracer.reHitTolerance : 0;
			
		//work out if each contact point is on the actual ray and hits within the actual cylinder length
		boolean p0OnCylinder = t[0] > reHitTolerance 
								&& d0 >= -length/2 && d0 <= length/2
								&& u0 >= uMin && u0 < uMax
								&& r0 >= rMin && r0 < rMax;
		boolean p1OnCylinder = t[1] > reHitTolerance 
								&& d1 >= -length/2 && d1 <= length/2
								&& u1 >= uMin && u1 < uMax
								&& r1 >= rMin && r1 < rMax;
		
		//Intersection hit = new Intersection();
		double d;
		boolean hitOnExit;
		if(p0OnCylinder && (t[0] < t[1] || !p1OnCylinder)){
			//t[0] is the first or only contact with cylinder
			if(t[0] > ray.length)
				return false;
			
			hit.pos = P0;
			ray.length = t[0];
			d = d0;
			
			hitOnExit = t[0] > t[1];
			
		}else if(p1OnCylinder && (t[1] < t[0] || !p0OnCylinder)){
			//t[1] is the first or only contact with cylinder
			if(t[1] > ray.length)
				return false;
			
			hit.pos = P1;
			ray.length = t[1];
			d = d1;
			
			hitOnExit = t[1] > t[0];
		}else{
			return false;
		}
		
		hit.normal = Util.reNorm(new double[]{ // P - C - d.N
				hit.pos[0] - centre[0] - d * axis[0],
				hit.pos[1] - centre[1] - d * axis[1],
				hit.pos[2] - centre[2] - d * axis[2],
		});
		
		hit.surface = this;
		return true;

	}
	
	public void setDrawingDetails(int nSections, int nPointsPerSection){
		this.nPointsPerSection = nPointsPerSection;
		this.nSections = nSections;
	}
	
	@Override
	public List<double[][]> draw() {
		int nPointsPerSection = (this.nPointsPerSection > 0) ? this.nPointsPerSection : 5 + (approxDrawQuality * 10 / 100); //min 5, max 55  at def(10)
		int nSections = (this.nSections > 0) ? this.nSections : 1 + (approxDrawQuality * 10);
		
		ArrayList<double[][]> lines = new ArrayList<double[][]>();
		

		for(int i=0; i < nSections; i++){
			double phi0 = i * 2 * Math.PI / nSections;
			double phi1 = (i+1) * 2 * Math.PI / nSections;
			
			double dPhi = (phi1 - phi0) / (nPointsPerSection - 1);
			double line[][] = new double[3][nPointsPerSection*2 + 1];
			int nPointsOutOfRange = 0;
			for(int j=0; j < nPointsPerSection; j++){
				double subPhi = j * dPhi;

				double u = radius * Math.cos(phi0+subPhi);
				double r = radius * Math.sin(phi0+subPhi);
				if(u < uMin || u > uMax || r < rMin || r > rMax)
					nPointsOutOfRange++;
				
				u = FastMath.max(u, uMin);
				u = FastMath.min(u, uMax);
				r = FastMath.max(r, rMin);
				r = FastMath.min(r, rMax);
				
				for(int k=0; k < 3; k++){
					line[k][j] = centre[k] - length/2 * axis[k]
				                           + u * up[k]
				                           + r * right[k];
				}

				u = radius * Math.cos(phi1-subPhi);
				r = radius * Math.sin(phi1-subPhi);
				if(u < uMin || u > uMax || r < rMin || r > rMax)
					nPointsOutOfRange++;
				
				u = FastMath.max(u, uMin);
				u = FastMath.min(u, uMax);
				r = FastMath.max(r, rMin);
				r = FastMath.min(r, rMax);
				
				for(int k=0; k < 3; k++){
					line[k][nPointsPerSection+j] = centre[k] + length/2 * axis[k]
				                           + u * up[k]
				                           + r * right[k];
				}
			}

			line[0][nPointsPerSection*2] = line[0][0];
			line[1][nPointsPerSection*2] = line[1][0];
			line[2][nPointsPerSection*2] = line[2][0];
			
			if(nPointsOutOfRange == 0)
				lines.add(line);
		}
		
		return lines;
			
		
	}

	@Override
	public void shift(double[] dX) {
		centre[0] += dX[0];
		centre[1] += dX[1];
		centre[2] += dX[2];
	}

	@Override
	public void rotate(double[] point, double[][] matrix) {

		for(int i=0;i<3;i++)
			centre[i] -= point[i];
			
		double newCentre[] = new double[3];
		double newAxis[] = new double[3];
		
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++){
				newCentre[i] += matrix[i][j] * centre[j];
				newAxis[i] += matrix[i][j] * axis[j];
			}
				
		for(int i=0;i<3;i++)
			centre[i] = point[i] + newCentre[i];
		
		axis = newAxis;
	}

	public double[] getCentre() { return centre.clone(); }
	
	public void setCentre(double[] centre) { this.centre = centre;	}
	
	public double[] getAxis() { return axis.clone(); }

	public double getRadius() { return radius; }
	public double getLength() { return length; }
	public void setLength(double length){ this.length = length; }
	public void setRadius(double radius){ this.radius = radius; }

	@Override
	public int surfaceGeometryHashCode() {
		long t; int r = 1;
		r = 31 * r + Arrays.hashCode(axis);
		r = 31 * r + Arrays.hashCode(centre);
		t = Double.doubleToLongBits(length);	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(radius); 	r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}

	@Override
	public boolean surfaceGeometryEquals(Surface obj) {
		TruncatedCylinder other = (TruncatedCylinder) obj;
		return Arrays.equals(axis, other.axis)
				&& Arrays.equals(centre, other.centre)
				&& Double.doubleToLongBits(length) == Double.doubleToLongBits(other.length)
				&& Double.doubleToLongBits(radius) == Double.doubleToLongBits(other.radius);
	}
	
	
}
