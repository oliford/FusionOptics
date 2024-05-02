package fusionOptics.optics;

import java.io.FileReader;
import java.io.IOException;
import java.util.Scanner;

import fusionOptics.Util;
import fusionOptics.interfaces.Absorber;
import fusionOptics.surfaces.Triangle;
import fusionOptics.types.Optic;
import uk.co.oliford.jolu.OneLiners;

public class NodesAndElementsMesh extends Optic {
	
	public NodesAndElementsMesh(String name, String filePrefix) {
		super(name);
		mustAddTrianglesInRadius(this, filePrefix + ".elements", filePrefix + ".nodes", null, 0);
	}
	
	public NodesAndElementsMesh(String name, String elementsFileName, String nodesFileName) {
		super(name);
		mustAddTrianglesInRadius(this, elementsFileName, nodesFileName, null, 0);
	}
	
	public final static void addTriangles(Optic optic, String elementsFileName, String nodesFileName) throws IOException{
		addTrianglesInRadius(optic, elementsFileName, nodesFileName, null, 0);
	}
	
	public NodesAndElementsMesh(String name, String filePrefix, double[] centre, double radius) {
		super(name);
		mustAddTrianglesInRadius(this, filePrefix + ".elements", filePrefix + ".nodes", centre, radius);
	}
	public final static void addTriangles(Optic optic, String elementsFileName, String nodesFileName, double[] centre, double radius) throws IOException{
		addTrianglesInRadius(optic, elementsFileName, nodesFileName, centre, radius);
	}
	
	public final static void mustAddTrianglesInRadius(Optic optic, String elementsFileName, String nodesFileName, double[] centre, double radius){
		try {
			addTrianglesInRadius(optic, elementsFileName, nodesFileName, centre, radius);
		} catch (IOException e) { 
			throw new RuntimeException(e);
		}
		
	}
	
	public final static void addTrianglesInRadius(Optic optic, String elementsFileName, String nodesFileName, double[] centre, double radius) throws IOException{
		
		//******************************** Elements file initialization ********************************
		
		Scanner inEl = new Scanner(new FileReader(elementsFileName));
		int totElLines = 0;
		int iEl = 0;
		
		// count the number of lines in the elements file
		while(inEl.hasNextLine()){
			inEl.nextLine();
			totElLines ++;
		}
		inEl.close();

		int[][] Elements = new int[totElLines][4];
		
		// insert elements values in a 2D-Array
		inEl = new Scanner(new FileReader(elementsFileName));
		while(iEl<totElLines-1){
			Elements[iEl][0] = inEl.nextInt();
			Elements[iEl][1] = inEl.nextInt();
			Elements[iEl][2] = inEl.nextInt();
			Elements[iEl][3] = inEl.nextInt();
			iEl ++;
		}
		inEl.close();
		
		///********************************** Nodes file initialization **********************************
		
		Scanner inNo = new Scanner(new FileReader(nodesFileName));
		int totNoLines = 0;
		int iNo = 0;
		
		// count the number of lines in the nodes file
		while(inNo.hasNextLine()){
			inNo.nextLine();
			totNoLines ++;
		}
		inNo.close();

		float[][] Nodes = new float[totNoLines][4];
		
		// insert nodes values in a 2D-Array
		inNo = new Scanner(new FileReader(nodesFileName));
		while(iNo<totNoLines-1){
			Nodes[iNo][0] = inNo.nextFloat();
			Nodes[iNo][1] = inNo.nextFloat();
			Nodes[iNo][2] = inNo.nextFloat();
			Nodes[iNo][3] = inNo.nextFloat();
			iNo ++;
		}
		inNo.close();
		
		//************************************* Triangles creation ****************************************
		
		// create triangles using fusionOptics.surfaces.Triangle and set them as ideal absorber
		/** (please note that here i used the laziest solution: a correct implementation of this method 
		 * should not rely on the fact that the node id is equal to the node line number + 1)*/
		iEl=0;
		while(iEl<totElLines-1){
			int idxA = findNode(Nodes, Elements[iEl][1]); double[] A = new double[]{Nodes[idxA][1],Nodes[idxA][2],Nodes[idxA][3]};
			int idxB = findNode(Nodes, Elements[iEl][2]); double[] B = new double[]{Nodes[idxB][1],Nodes[idxB][2],Nodes[idxB][3]};
			int idxC = findNode(Nodes, Elements[iEl][3]); double[] C = new double[]{Nodes[idxC][1],Nodes[idxC][2],Nodes[idxC][3]};
				if(centre == null || 
					((Util.length(Util.minus(A, centre)) < radius)&&(Util.length(Util.minus(B, centre)) < radius)&&(Util.length(Util.minus(C, centre)) < radius))){
				Triangle triangle = new Triangle(String.valueOf(Elements[iEl][0]), A, B, C, Absorber.ideal());
				optic.addElement(triangle);
			}
			iEl ++;
		}
	}
	
	private final static int findNode(float Nodes[][], int id){
		if(Nodes[id-1][0] == id)
			return id-1;
		for(int i=0; i < Nodes.length; i++){
			if(Nodes[i][0] == id)
				return i;
		}
		throw new RuntimeException("Node " + id + " not found.");
		
	}

}
