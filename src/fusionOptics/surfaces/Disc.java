package fusionOptics.surfaces;

import java.util.ArrayList;
import java.util.List;

import fusionOptics.tracer.Tracer;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.RaySegment;


public class Disc extends Plane {
	
	double radius;
	public Disc(String name, double centre[], double normal[], double radius, Interface iface) {
		this(name, centre, normal, radius, null, null, iface);
		
	}
		
	public Disc(String name, double centre[], double normal[], double radius, Medium frontMedium, Medium backMedium, Interface iface) {
		super(name, centre, normal, frontMedium, backMedium, iface);
		this.radius = radius;
		
		initArbitraryPerps();
	}
	

	@Override
	public double[] getBoundarySphereCentre() { return centre.clone(); }

	@Override
	public double getBoundarySphereRadius() { return radius; }
	
	public boolean findEarlierIntersection(RaySegment ray, Intersection hit){
		double pos[] = new double[3];
		double dist = super.calcPlaneIntersection(ray, pos); //find the place where the ray hits us (as if we were an infinite plane)
		
		//Check that it actually hit the plane, (not before, not under tolerance, not over length, and not parallel)
		if(dist < -Tracer.reHitTolerance || dist > ray.length || Double.isNaN(dist)) 
			return false;
		if(ray.startHit != null && ray.startHit.surface == this && dist < Tracer.reHitTolerance)
			return false; //rehit us within tolerance, ignore
				
		double sum=0;
		
		for(int i=0;i<3;i++)
			sum += (centre[i] - pos[i]) * (centre[i] - pos[i]);
		
		
		if( sum <= (radius*radius)){
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
		int nPointsOnCirum = 45 + (approxDrawQuality * 500 / 100); //min 45, max 545, 95 at def(10)
		int nDiameterLines = 3 + (approxDrawQuality * 30 / 100); // min 3, max 33, 6 at def(10)
			
		ArrayList<double[][]> lines = new ArrayList<double[][]>();
		
		//line around disc circumference
		double line[][] = new double[3][nPointsOnCirum];
		for(int i=0; i < nPointsOnCirum; i++){
			double phi = i * 2 * Math.PI / (nPointsOnCirum - 1);
			
			for(int k=0; k < 3; k++)
				line[k][i] = centre[k] + radius * (   0 * normal[k] 
			                                       	+ Math.cos(phi) * up[k]
			                                       	+ Math.sin(phi) * right[k]);                                  			
		}
		lines.add(line);
		
		for(int i=0; i < nDiameterLines; i++){
			double phi = i * Math.PI / nDiameterLines;
			line = new double[3][2];
			
			for(int k=0; k < 3; k++){
				line[k][0] = centre[k] - radius * (   0 * normal[k] 
									                          + Math.cos(phi) * up[k]
									                          + Math.sin(phi) * right[k]);
				line[k][1] = centre[k] + radius * (   0 * normal[k] 
									                          + Math.cos(phi) * up[k]
									                          + Math.sin(phi) * right[k]);
				}
			
			lines.add(line);
		}
		
		
		return lines;
	}

	public final double getRadius() { return radius; }

	@Override
	protected int planeBoundaryHashCode() {
		long temp = Double.doubleToLongBits(radius);
		return 31 + (int) (temp ^ (temp >>> 32));		
	}

	@Override
	protected boolean planeBoundaryEquals(Plane obj) {
		return (this == obj) || ( super.equals(obj) 
					&& Double.doubleToLongBits(radius) == Double.doubleToLongBits(((Disc)obj).radius) );
	}
	
}
