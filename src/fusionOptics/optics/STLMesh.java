package fusionOptics.optics;

import java.io.IOException;

import otherSupport.BinarySTLFile;
import otherSupport.BinarySTLFile.Triangles;
import fusionOptics.Util;
import fusionOptics.interfaces.Absorber;
import fusionOptics.surfaces.Triangle;
import fusionOptics.types.Optic;

public class STLMesh  extends Optic {
			
			public STLMesh(String fileName) {				
				super(fileName.replaceAll(".*/", ""));
				mustAddTrianglesInRadius(this, fileName, null, 0);
			}
			
			public STLMesh(String name, String fileName) {
				super(name);
				mustAddTrianglesInRadius(this, fileName, null, 0);
			}
			
			public final static void addTriangles(Optic optic, String fileName) throws IOException{
				addTrianglesInRadius(optic, fileName, null, 0);
			}
			
			public STLMesh(String name, String fileName, double[] centre, double radius) {
				super(name);
				mustAddTrianglesInRadius(this, fileName, centre, radius);
			}
			public final static void addTriangles(Optic optic, String fileName, double[] centre, double radius) throws IOException{
				addTrianglesInRadius(optic, fileName, centre, radius);
			}
			
			public final static void mustAddTrianglesInRadius(Optic optic, String fileName, double[] centre, double radius){
				try {
					addTrianglesInRadius(optic, fileName, centre, radius);
				} catch (IOException e) { 
					throw new RuntimeException(e);
				}
				
			}
			
			public final static void addTrianglesInRadius(Optic optic, String fileName, double[] centre, double radius) throws IOException{
				
				Triangles triangles = BinarySTLFile.mustRead(fileName);
				double scale= 0.001;
				for(int i=0; i < triangles.count; i++){
					
					for(int j=0; j < 3; j++){
						triangles.vertex1[i][j] *= scale;
						triangles.vertex2[i][j] *= scale;
						triangles.vertex3[i][j] *= scale;
					}
					
					if(centre == null || ((Util.length(Util.minus(triangles.vertex1[i], centre)) < radius) && 
										  (Util.length(Util.minus(triangles.vertex2[i], centre)) < radius) &&
										  (Util.length(Util.minus(triangles.vertex3[i], centre)) < radius))){
						
							Triangle triangle = new Triangle("triangle_" + i, triangles.vertex1[i], triangles.vertex2[i], triangles.vertex3[i], Absorber.ideal());
							optic.addElement(triangle);			
					}
				}
				
				
				
			}
				
}
