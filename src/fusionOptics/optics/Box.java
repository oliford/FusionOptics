package fusionOptics.optics;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import fusionOptics.surfaces.Square;
import fusionOptics.types.Element;
import fusionOptics.types.Interface;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;



public class Box extends Optic {
	
	//private double centre[];
	//private double maxRadius;
	
	public Box(String name, double centre[], double wx, double wy, double wz, Medium innerMedium, Interface iface) {
		this(name, 
				new double[]{ centre[0] - wx/2, centre[1] - wy/2, centre[2] - wz/2 },
				new double[]{ centre[0] + wx/2, centre[1] + wy/2, centre[2] + wz/2 },
				innerMedium, iface);
	}
	
	/** Create box aligned to axes from corner1 to corner2.
	 *  The first optic axis of the material, if it has one, will be aligned to the x axis, and the second to the y axis */
	public Box(String name, double corner1[], double corner2[], Medium innerMedium, Interface iface) {
		super(name);
		
		double centre[] = new double[]{
			(corner1[0] + corner2[0]) / 2.0,
			(corner1[1] + corner2[1]) / 2.0,
			(corner1[2] + corner2[2]) / 2.0
		};
		
		double w = corner2[0] - corner1[0]; //width (x) 
		double d = corner2[1] - corner1[1]; //depth (y)
		double h = corner2[2] - corner1[2]; //height (z)
		
		addElement(new Square("face+x", new double[]{ centre[0] + (w/2.0), centre[1], centre[2] }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, d, h, null, innerMedium, iface)); // right (+x)
		addElement(new Square("face-x", new double[]{ centre[0] - (w/2.0), centre[1], centre[2] }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, d, h, null, innerMedium, iface)); // left (-x)
		addElement(new Square("face+y", new double[]{ centre[0], centre[1] + (d/2.0), centre[2] }, new double[]{ 0, 1, 0 }, new double[]{ 1, 0, 0 }, w, h, null, innerMedium, iface)); // back (+y)
		addElement(new Square("face-y", new double[]{ centre[0], centre[1] - (d/2.0), centre[2] }, new double[]{ 0, -1, 0 }, new double[]{ 1, 0, 0 }, w, h, null, innerMedium, iface)); // front (-y)
		addElement(new Square("face+z", new double[]{ centre[0], centre[1], centre[2] + (h/2.0) }, new double[]{ 0, 0, 1 }, new double[]{ 1, 0, 0 }, w, d, null, innerMedium, iface)); // top (+z)
		addElement(new Square("face-z", new double[]{ centre[0], centre[1], centre[2] - (h/2.0) }, new double[]{ 0, 0, -1 }, new double[]{ 1, 0, 0 }, w, d, null, innerMedium, iface)); // bottom (-z)

	}
	
}
