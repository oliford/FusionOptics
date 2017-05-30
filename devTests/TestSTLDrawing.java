import fusionOptics.drawing.STLDrawer;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.SimplePlanarConvexLens;
import fusionOptics.surfaces.Disc;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.RaySegment;


public class TestSTLDrawing {

	public static void main(String[] args) {
		//Disc disc = new Disc("disc", new double[]{ 0, 0, 0 }, new double[]{ 1, 0, 0}, 1.00, NullInterface.ideal());
		SimplePlanarConvexLens lens = SimplePlanarConvexLens.fromFocalLengthAndCentreThickness("lens", 
				new double[]{ 0, 0, 0 },
				new double[]{ 1, 0, 0}, 
				0.050,
				0.300,
				0.010, 
				new Medium(new IsotropicFixedIndexGlass(1.5)), IsoIsoInterface.ideal(), 500e-9);
		Optic sys = new Optic("sys", new Element[]{ lens });
		STLDrawer stlDrawer = new STLDrawer("/tmp/disc.stl");
		stlDrawer.drawOptic(sys);
		stlDrawer.destroy();
		
	}
}
