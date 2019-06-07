package fusionOptics;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import fusionOptics.interfaces.NullInterface;
import fusionOptics.surfaces.Iris;
import fusionOptics.types.Element;
import fusionOptics.types.Material;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;

import net.jafama.FastMath;

/** Small misc static helper methods (mostly vector maths) */
public abstract class Util {
	
	public final static void rotateOnX(Element element, double point[], double ang){
		double mat[][] = new double[][]{ 	
				{1, 0, 0},
				{0, Math.cos(ang), -Math.sin(ang)},
				{0, Math.sin(ang), Math.cos(ang)} };
				
		element.rotate(point,mat);
	}
	public final static void rotateOnY(Element element, double point[], double ang){
		double mat[][] = new double[][]{ 	
				{Math.cos(ang), 0, -Math.sin(ang)},
				{0, 1, 0},
				{Math.sin(ang), 0, Math.cos(ang)} };
				
		element.rotate(point,mat);
	}
	public final static void rotateOnZ(Element element, double point[], double ang){
		double mat[][] = new double[][]{ 	
				{Math.cos(ang), -Math.sin(ang), 0},
				{ Math.sin(ang), Math.cos(ang), 0},
				{0, 0, 1} };				
		element.rotate(point,mat);
	}	
		
	public final static double length(double x[]){ 
		return FastMath.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); 
	}
		
	/** Renormalise the given vector */
	public final static double[] reNorm(double x[]){
		double sum=0;
		
		for(int i=0;i<3;i++)
			sum += x[i] * x[i];
		
		sum = FastMath.sqrt(sum);
		
		for(int i=0;i<3;i++)
			x[i] /= sum;
		
		return x;
	}
	
	/** Calculates the cross production of A and B*/
	public final static double[] cross(double A[], double B[]){
		double AxB[] = new double[3];
		
		AxB[0] =  A[1] * B[2] - A[2] * B[1];
		AxB[1] = -A[0] * B[2] + A[2] * B[0];
		AxB[2] =  A[0] * B[1] - A[1] * B[0];
		
		return AxB;
	}

	/** Dot product */
	public final static double dot(double A[], double B[]){
		return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
	}
	
	/** @return a - b */
	public final static double[] minus(double[] a, double[] b) { 
		return new double[]{ a[0] - b[0], a[1] - b[1], a[2] - b[2] };
	}
	
	/** @return a + b */
	public final static double[] plus(double[] a, double[] b) { 
		return new double[]{ a[0] + b[0], a[1] + b[1], a[2] + b[2] };
	}
	
	/** @return a * b */
	public final static double[] mul(double[] a, double b) { 
		return new double[]{ a[0] * b, a[1] * b, a[2] * b };
	}
	
	/** Create vector perpendicular to given vector. */ 
	public final static double[] createPerp(double a[]){
		if(a[1] == 0 && a[2] == 0)
			return new double[]{ 0, 0, 1 };
		else
			return reNorm(cross(a, new double[]{1,0,0} ));
	}
	
	/** Find an element by name in an arrays of elements */
	public final static Element findElement(Element elems[], String name){
		for(int i=0;i<elems.length;i++)
			if(name.equals(elems[i].getName()))
				return elems[i];
		
		return null;
	}
	
	/** Find an element by name within an optics subOptics/surfaces tree */
	public final static Element findElement(Optic optic, String name){
		for(Surface surface : optic.getSurfaces()) {
			if(name.equals(surface.getName()))
				return surface;
		}
		for(Optic subOptic : optic.getSubOptics()) {
			Element elem = findElement(subOptic, name);
			if(elem != null)
				return elem;
		}
		
		return null;
	}
	
	/** Convert the optic drawn line/list lists to points 
	 * 
	 * @param drawing
	 * @param lineDL	Maximum line segement length. Bigger ones will be broken.
	 * @return
	 */
	public static double[][] drawingToPoints(List drawing, double lineDL) {
		LinkedList<double[]> pts = new LinkedList<double[]>();
		drawToPtsRecursion(pts, drawing, lineDL);
		double ret[][] = new double[pts.size()][];
		int k = 0;
		for(double p[] : pts)
			ret[k++] = p;
		return ret;
	}
	
	
	private static void drawToPtsRecursion(List<double[]> pts, List drawing, double lineDL){
		
		for(Object o : drawing){
			
			if(o instanceof double[][]){
				double p[][] = (double[][])o;
				
				for(int i=0; i < (p[0].length-1); i++){
					double pl = FastMath.sqrt(
						FastMath.pow2(p[0][i+1] - p[0][i]) +
						FastMath.pow2(p[1][i+1] - p[1][i]) +
						FastMath.pow2(p[2][i+1] - p[2][i]));
				
					if(pl < lineDL){
						pts.add(new double[]{ p[0][i], p[1][i], p[2][i] });
					}else{ //need to break up  the line
						int nPoints = ((int)(pl / lineDL + 0.5)) + 1;
						double d = 1.0 / (nPoints - 1);
						for(int j=0; j < nPoints; j++){
							pts.add(new double[]{ 
									p[0][i] + j * (p[0][i+1] - p[0][i]) * d,
									p[1][i] + j * (p[1][i+1] - p[1][i]) * d,
									p[2][i] + j * (p[2][i+1] - p[2][i]) * d });
						}
					}
				}
					
				
			}else if(o instanceof List){
				drawToPtsRecursion(pts, (List)o, lineDL);
			}else{
				throw new IllegalArgumentException("Drawing should be tree of Lists and double[]s");
			}
		}
		
	}
	
	/** Calculates the required thickness of a waveplate to give a single 
	 * complete wavelength difference at the given frequency.
	 * (This assumes a ray travelling perpendicular to the optic axis)  
	 * 
	 * @param mat
	 * @param wavelength
	 * @return
	 */
	public static double calcWaveplateFullWaveThickness(Material mat, double wavelength) {
		
		double wavePlateNO = mat.getRefractiveIndex(0, wavelength, 300);
		double wavePlateNE = mat.getRefractiveIndex(1, wavelength, 300);
		return wavelength / (FastMath.abs(wavePlateNE - wavePlateNO));
	}
	
	/** Returns the bounding box { x0,y0,z1, x1,y1,z1 } of the min/max coordinates of the given element */
	public static double[] getBoundingBox(Element elem) {
		double c[] = elem.getBoundarySphereCentre();
		double r = elem.getBoundarySphereRadius();
		
		return new double[]{ 
				c[0]-r, c[1]-r, c[2]-r, 
				c[0]+r, c[1]+r, c[2]+r
			};
	}
	
	/** Mean and covariance mat of vector 
	 * @return double[]{ <x>, <y>, <z>, <θ>, <x.x>, <x.y>, <x.z>, <y.y>, <y.z>, <z.z> */
	public static double[] vecStats(double[][] vec, int n) {
		double stats[] = new double[10];
		
		for(int i=0; i < n; i++){
			stats[0] += vec[i][0];
			stats[1] += vec[i][1];
			stats[2] += vec[i][2];
		}
		stats[0] /= n;
		stats[1] /= n;
		stats[2] /= n;
		
		for(int i=0; i < n; i++){
			int e = 4;
			stats[3] += FastMath.acos(Util.dot(vec[i], stats));
			for(int j=0; j < 3; j++){
				for(int k=j; k < 3; k++){
					stats[e++] += (vec[i][j] - stats[j]) * (vec[i][k] - stats[k]);
				}
			}
		}
		for(int i=3; i < 10; i++)
			stats[i] /= (n - 1);
		
		return stats;
	}
	
	/** Strip elements with NullInterface
	 * 
	 * Such surfaces are usually only used for debugging and can be thrown
	 * away
	 *  
	 * @param surfList
	 */
	public static void stripNullInterfaceSurfaces(LinkedList<Surface> surfList) {
		Iterator<Surface> iter = surfList.iterator();
		while(iter.hasNext()){							
			if(iter.next().getInterface() instanceof NullInterface)
				iter.remove();
		}
	}

	/** Moves all Irises in a given optic tree a long way away.
	 * Useful for pretty drawings of lenses without their ray catching irises. 
	 * 
	 *  @param shiftOut How far away to throw the iris,
	 */
	public static void throwAwayIrises(Optic optic, double shiftOut){
		for(Optic o : optic.getSubOptics())
			throwAwayIrises(o, shiftOut);
		
		for(Surface s : optic.getSurfaces())
			if(s instanceof Iris)
				s.shift(new double[]{ shiftOut, 0, 0 });
	}
	
	
}
