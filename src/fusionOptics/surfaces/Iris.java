package fusionOptics.surfaces;

import java.util.ArrayList;
import java.util.List;

import fusionOptics.tracer.Tracer;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.RaySegment;


public class Iris extends Plane {
	
	double discRadius;
	double apatureRadius;
	
	public Iris(String name, double centre[], double normal[], double discRadius, double apatureRadius, Interface iface) {
		this(name, centre, normal, discRadius, apatureRadius, null, null, iface);
	}
		
	public Iris(String name, double centre[], double normal[], double discRadius, double apatureRadius, Medium frontMedium, Medium backMedium, Interface iface) {
		super(name, centre, normal, frontMedium, backMedium, iface);
		this.discRadius = discRadius;
		this.apatureRadius = apatureRadius;
		
		initArbitraryPerps();
	}
	
	@Override
	public double[] getBoundarySphereCentre() { return centre.clone(); }

	@Override
	public double getBoundarySphereRadius() { return discRadius; }

	@Override
	public boolean findEarlierIntersection(RaySegment ray, Intersection hit) {
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
				
		if( sum <= (discRadius*discRadius) && sum >= (apatureRadius*apatureRadius)){ //if inside the disc radius but outside the apature radius
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
		int nPointsPerSection = 11;
		int nSections = 8;
		
		ArrayList<double[][]> lines = new ArrayList<double[][]>();
		
		for(int i=0; i < nSections; i++){
			double phi0 = i * 2 * Math.PI / nSections;
			double phi1 = (i+1) * 2 * Math.PI / nSections;
			
			double dPhi = (phi1 - phi0) / (nPointsPerSection - 1);
			
			double line[][] = new double[3][nPointsPerSection*2 + 1];
			for(int j=0; j < nPointsPerSection; j++){
				double subPhi = j * dPhi;
				
				for(int k=0; k < 3; k++)
					line[k][j] = centre[k] + discRadius * (   0 * normal[k] 
				                                       	+ Math.cos(phi0+subPhi) * up[k]
				                                       	+ Math.sin(phi0+subPhi) * right[k]); 
				              
				for(int k=0; k < 3; k++)
					line[k][nPointsPerSection+j] = centre[k] + apatureRadius * (   0 * normal[k] 
				                                       	+ Math.cos(phi1-subPhi) * up[k]
				                                       	+ Math.sin(phi1-subPhi) * right[k]); 
				                                 
			}

			line[0][nPointsPerSection*2] = line[0][0];
			line[1][nPointsPerSection*2] = line[1][0];
			line[2][nPointsPerSection*2] = line[2][0];
			
			lines.add(line);
		}
		
		return lines;
	}

	public void setApatureRadius(double apatureRadius) { this.apatureRadius = apatureRadius; }
	public double getApatureRadius() { return this.apatureRadius; }
	public void setDiscRadius(double discRadius) { this.discRadius = discRadius; }
	public double getDiscRadius() { return this.discRadius; }

	public void setCentre(double[] centre) { this.centre = centre;	}

	@Override
	protected int planeBoundaryHashCode() {
		long t; int r = 1; 
		t = Double.doubleToLongBits(discRadius);	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(apatureRadius); r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}

	@Override
	protected boolean planeBoundaryEquals(Plane obj) {
		return Double.doubleToLongBits(discRadius) == Double.doubleToLongBits(((Iris)obj).discRadius)
				&& Double.doubleToLongBits(apatureRadius) == Double.doubleToLongBits(((Iris)obj).apatureRadius);
	}

}
