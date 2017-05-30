package fusionOptics.collection;

import fusionOptics.surfaces.Plane;
import fusionOptics.types.Intersection;
import binaryMatrixFile.BinaryMatrixFile;
import oneLiners.OneLiners;

/** Intersection processor which simply builds and image
 * which can be written to a binary file. 
 *  
 * @author oliford */
public class ImageCollector implements IntersectionProcessor {

	public ImageCollector(double x0, double x1, int nX, double y0, double y1, int nY) {
		this.x0 = x0; this.x1 = x1; this.nX = nX; this.dx = (x1 - x0) / (nX - 1);
		this.y0 = y0; this.y1 = y1; this.nY = nY; this.dy = (y1 - y0) / (nY - 1);
		image = new double[nY][nX];
	}
	
	private int nX,nY;
	private double x0,x1,dx;
	private double y0,y1,dy;
	private double image[][];
	
	@Override
	public void nextIntersection(Intersection hit) {
		Plane imagePlane = (Plane)hit.surface;
		double imgPos[] = imagePlane.posXYZToPlaneUR(hit.pos);
		
		int iX = (int)((imgPos[0] - x0)/dx);
		int iY = (int)((imgPos[1] - y0)/dy);
		
		if(iX >= 0 && iX < nX && iY >=0 && iY < nY)
			image[iY][iX] += hit.incidentRay.endIntensity();
		
	}
	
	public double[][] getImage() { return image; }

	public void writeImage(String fileName) {
		BinaryMatrixFile.mustWrite(fileName, 
							OneLiners.linSpace(x0, x1, dx),
							OneLiners.linSpace(y0, y1, dy),
							image, false);
	}
}
