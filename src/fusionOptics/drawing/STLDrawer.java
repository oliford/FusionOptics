package fusionOptics.drawing;

import java.lang.annotation.ElementType;
import java.util.ArrayList;
import java.util.List;

import fusionOptics.Util;
import fusionOptics.types.Element;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;
import uk.co.oliford.jolu.BinarySTLFile;
import uk.co.oliford.jolu.OneLiners;
import uk.co.oliford.jolu.BinarySTLFile.Triangles;

/** Make STL MESH from drawings (surfaces only)
 */
public class STLDrawer {
	
	private String fileName;
	private ArrayList<double[][]> triangles;
	private ArrayList<Class<Element>> ignoreElementTypes = new ArrayList<Class<Element>>();
	
	private double globalRotMat[][] = null;
	
	
	public STLDrawer(String fileName) {
		this.fileName = fileName;
		this.triangles = new ArrayList<double[][]>();
	}

	public void drawOptic(Optic optic) {		
		passTree(optic);
	}
	
	private void passTree(Optic optic) {		
		for(Optic subOptic : optic.getSubOptics()){
			passTree(subOptic);
		}
		
nextSurf:
		for(Surface surface : optic.getSurfaces()){
			for(Class<Element> elementType : ignoreElementTypes){
				if(elementType.isInstance(surface))
					continue nextSurf;
			}
			drawSurface(surface);
		}		
	}
	
	public void drawSurface(Surface surface){
		for(double polygon[][] : surface.draw())
			drawPolygon(polygon);
	}
	
	public void ignoreElementType(Class type){
		ignoreElementTypes.add(type);
	}
	
	private void drawPolygon(double[][] vertices) {
		//assume convex for now
		double v[][] = OneLiners.transpose(vertices);
		
		int i0 = 0;		
		for(int i2=2; i2 < v.length; i2++){
			int i1 = i2-1;
			
			triangles.add(new double[][]{ v[i0], v[i1], v[i2] });			
		}
	}

	public void destroy(){
		Triangles stlTriangles = new Triangles(triangles.size());
		
		for(int i=0; i < triangles.size(); i++){
			double triangle[][] = triangles.get(i);			
			stlTriangles.vertex1[i] = rotVec(triangle[0]);			
			stlTriangles.vertex2[i] = rotVec(triangle[1]);			
			stlTriangles.vertex3[i] = rotVec(triangle[2]);
			stlTriangles.normal[i] = Util.reNorm(rotVec(Util.cross(Util.minus(triangle[1], triangle[0]), Util.minus(triangle[2], triangle[0]))));
		}
		
		BinarySTLFile.mustWrite(fileName, stlTriangles);
	}
	

	private final double[] rotVec(double v[]){
		if(globalRotMat == null)
			return v;
		
		return new double[]{
				globalRotMat[0][0] * v[0] + globalRotMat[0][1] * v[1] + globalRotMat[0][2] * v[2],
				globalRotMat[1][0] * v[0] + globalRotMat[1][1] * v[1] + globalRotMat[1][2] * v[2],
				globalRotMat[2][0] * v[0] + globalRotMat[2][1] * v[1] + globalRotMat[2][2] * v[2], 
		};
	}

	/** Sets a global rotation matrix to be applied to all drawing */ 
	public void setTransformationMatrix(double[][] rotationMatrix) {
		this.globalRotMat = rotationMatrix;
	}

	public void drawElement(Element elem) {
		if(elem instanceof Optic){
			drawOptic((Optic)elem);
		}else{
			drawSurface((Surface)elem);
		}
	}
	
	
}
