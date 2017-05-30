package fusionOptics.optics;

import java.util.ArrayList;

import algorithmrepository.Algorithms;
import fusionOptics.Util;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.materials.IsotropicLinearDispersiveGlass;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
import fusionOptics.types.Element;
import fusionOptics.types.Interface;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;



/** f=100mm f/2 Focal length Double Gauss Lens from 
 * [ C.Kolb - "A Realistic Camera Model for Computer Graphics"
 *    http://www.cs.utexas.edu/~fussell/courses/cs395t/lens.pdf ]
 * 
 * Can be scaled to a different focal length, but is always f/2
 * 
 * @author oliford
 *
 */
public class DoubleGaussLens extends Optic {
	

	public Interface iface = new IsoIsoInterface(0.00);
	private Material[] mats = new Material[]{
			new IsotropicLinearDispersiveGlass(1.670, 47.1),
			new IsotropicLinearDispersiveGlass(1.670, 47.1),
			new IsotropicLinearDispersiveGlass(1.699, 30.1),
			
			new IsotropicLinearDispersiveGlass(1.603, 38.0),
			new IsotropicLinearDispersiveGlass(1.658, 57.3),
			new IsotropicLinearDispersiveGlass(1.717, 48.0)
		};
	
	public Medium media[] = new Medium[]{
			new Medium(mats[0]),
			new Medium(mats[1]),
			new Medium(mats[2]),
			
			new Medium(mats[3]),
			new Medium(mats[4]),
			new Medium(mats[5]),
	};

	private Dish dishes[] = new Dish[]{
			new Dish("s0", new double[]{ -0.033770, 0,0 }, new double[]{  1,0,0 }, 0.058950, 0.0504/2, media[0], null, iface),
			new Dish("s1", new double[]{ -0.026250, 0,0 }, new double[]{  1,0,0 }, 0.169660, 0.0504/2, null, media[0], iface),
			new Dish("s2", new double[]{ -0.026010, 0,0 }, new double[]{  1,0,0 }, 0.038550, 0.0460/2, media[1], null, iface),
			new Dish("s3", new double[]{ -0.017960, 0,0 }, new double[]{  1,0,0 }, 0.081540, 0.0460/2, media[2], media[1], iface),
			new Dish("s4", new double[]{ -0.011410, 0,0 }, new double[]{  1,0,0 }, 0.025500, 0.0360/2, null, media[2], iface),
			// iris goes here @ 33.770
			new Dish("s5", new double[]{  0.009000, 0,0 }, new double[]{ -1,0,0 }, 0.028990, 0.0340/2, null, media[3], iface),
			new Dish("s6", new double[]{  0.011360, 0,0 }, new double[]{  1,0,0 }, 0.081540, 0.0400/2, media[4], media[3], iface),
			new Dish("s7", new double[]{  0.023490, 0,0 }, new double[]{ -1,0,0 }, 0.040770, 0.0400/2, media[4], null, iface),
			new Dish("s8", new double[]{  0.023870, 0,0 }, new double[]{  1,0,0 }, 0.874130, 0.0400/2, media[5], null, iface), 
			new Dish("s9", new double[]{  0.030310, 0,0 }, new double[]{ -1,0,0 }, 0.079460, 0.0400/2, media[5], null, iface),
			// ?? @ 136.308
		};
	
	private Iris irises[] = new Iris[]{
			new Iris("apatureIris", new double[]{ 0.000000, 0, 0 }, new double[]{ 1,0,0 }, 0.0610/2, 0.0342/2, null, null, Absorber.ideal()),
			new Iris("iris1", new double[]{ -0.033770, 0, 0 }, new double[]{ 1,0,0 }, 0.0610/2, 0.0494/2, null, null, Absorber.ideal()),
			new Iris("iris2", new double[]{  0.030310, 0, 0 }, new double[]{ 1,0,0 }, 0.0610/2, 0.0390/2, null, null, Absorber.ideal()),
			new Iris("bigIris", new double[]{  0.0, 0, 0 }, new double[]{ 1,0,0 }, 0.2000/2, 0.0610/2, null, null, Absorber.ideal()),
		};
	
	private Cylinder cyld = new Cylinder("cyld", new double[]{-0.00173, 0, 0}, new double[]{1,0,0}, 0.0600/2, 0.06408, Absorber.ideal());
	
	
	public DoubleGaussLens() {
		this("DoubleGaussLens");
	}
	
	public DoubleGaussLens(String name, double centre[], double focalLength) {
		this(name, centre, null, focalLength);
	}
	
	public DoubleGaussLens(String name, double centre[], double axis[], double focalLength) {
		this(name);
		scale(focalLength / 0.100);
		shift(centre);
		
	
		if(axis != null){
			double y[] = Util.createPerp(axis);
			double z[] = Util.reNorm(Util.cross(axis, y));
			rotate(centre, Algorithms.rotationMatrix(new double[][]{ axis, y, z }));
		}
		
	}
	
	public DoubleGaussLens(String name) {
		super(name);
		
		for(Dish dish : dishes)
			addElement(dish);
		
		for(Iris iris : irises)
			addElement(iris);
		
		addElement(cyld);
	}
	
	public static void main(String[] args) {
		VRMLDrawer vrmlOut = new VRMLDrawer("/tmp/doubleGauss.vrml", 0.01);
		SVGRayDrawing svgOut = new SVGRayDrawing("/tmp/doubleGauss.svg", new double[]{ -0.1, -0.1, -0.1, 0.1, 0.1, 0.1 }, true);
		vrmlOut.drawElement(new DoubleGaussLens());
		svgOut.drawElement(new DoubleGaussLens());
		svgOut.destroy();
		vrmlOut.destroy();
	}
	
	public void scale(double scale){

		double lens0[] = irises[0].getCentre();

		for(Dish dish : dishes){
			double c[] = dish.getCentre();
			dish.setCentre(new double[]{
					lens0[0] + scale * (c[0] - lens0[0]),
					lens0[1] + scale * (c[1] - lens0[1]),
					lens0[2] + scale * (c[2] - lens0[2]),
				});

			dish.setRadiusOfCurv(dish.getRadiusOfCurv() * scale);
			dish.setDishDiameter(dish.getDishDiameter() * scale);
		}
		for(Iris iris : irises){
			double c[] = iris.getCentre();
			iris.setCentre(new double[]{
					lens0[0] + scale * (c[0] - lens0[0]),
					lens0[1] + scale * (c[1] - lens0[1]),
					lens0[2] + scale * (c[2] - lens0[2]),
				});
			
			iris.setApatureRadius(scale * iris.getApatureRadius());
			iris.setDiscRadius(scale * iris.getDiscRadius());
		}

		cyld.setLength(scale * cyld.getLength());
		cyld.setRadius(scale * cyld.getRadius());
		double c[] = cyld.getCentre();
		cyld.setCentre(new double[]{
				lens0[0] + scale * (c[0] - lens0[0]),
				lens0[1] + scale * (c[1] - lens0[1]),
				lens0[2] + scale * (c[2] - lens0[2]),
			});
	}

	/** Returns the radius of the lens case cylinder */
	public double getCaseRadius() { return cyld.getRadius(); }	
	/** Returns the length of the lens case cylinder */
	public double getCaseLength() { return cyld.getLength(); }
	
}
