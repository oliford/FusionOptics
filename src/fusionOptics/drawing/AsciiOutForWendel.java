package fusionOptics.drawing;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import fusionOptics.types.Optic;
import fusionOptics.types.RaySegment;

import oneLiners.OneLiners;



/** Outputs ascii format loaded by the Sergey Bozhenkov's 3D viewing program 'Wendel' 
 * @deprecated Replaced by VRMLViewer, since wendel reads those.
 * */
public class AsciiOutForWendel extends LineDrawer {
	private FileOutputStream fos;
	private PrintWriter writer;
	private String fileName;
	private int nLines = 0;
	private double currentColour[] = new double[]{ 0, 0, 0 };
	
	public AsciiOutForWendel(String fileName, boolean drawWholePaths, double smallLineLength) {
		super(drawWholePaths, smallLineLength);
		this.fileName = fileName;
		OneLiners.makePath(fileName);
		try {
			fos = new FileOutputStream(fileName);
			writer = new PrintWriter(fos);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
	}
	
	public void setCurrentColour(double[] currentColour) {
		this.currentColour = currentColour;
	}
	

	@Override
	public void destroy() {
		try {
			writer.flush();
			fos.flush();
			fos.close();
		} catch (IOException e) {
			System.err.println("WARNING: Error closing file '"+fileName+"'.");
		}
	}
	
	public void drawRay(RaySegment ray, double colour[]) {
		this.currentColour = colour;
		super.drawRay(ray);
	}
		
	public void drawElement(Optic optic, double colour[]) {
		this.currentColour = colour;
		super.drawElement(optic);
	}

	@Override
	public void drawSinglePath(List<double[]> path) {
		for(double vertex[] : path)
			addVertex(vertex[0], vertex[1], vertex[2], (vertex.length > 3) ? vertex[3] : 1.0);
		nLines++;
				
	}

	@Override
	public void drawSinglePath(double[][] path) {
		for(int i=0; i < path[0].length; i++)
			addVertex(path[0][i], path[1][i], path[2][i], (path.length > 3) ? path[3][i] : 1.0);
		nLines++;
	}
	
	private void addVertex(double x, double y, double z, double intensity) {
		writer.println(nLines + "\t" + x + "\t" + y+ "\t" + z + "\t" + 
						currentColour[0]*intensity + "\t"+currentColour[1]*intensity + "\t" + currentColour[2]*intensity);
		
	}

	@Override
	public void startGroup(String name) { }
	@Override
	public void endGroup() { }

}
