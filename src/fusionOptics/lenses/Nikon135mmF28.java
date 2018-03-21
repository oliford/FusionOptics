package fusionOptics.lenses;

import algorithmrepository.Algorithms;
import fusionOptics.Util;
import fusionOptics.optics.SequentialLensSeries;

/** This is a 100mm/2.8 from Nikkor patent 4,057,330 (1977), ex#2, upscaled to be f=135mm
 * as shown at http://www.pierretoscani.com/echo_shortpres.html */
public class Nikon135mmF28 extends SequentialLensSeries {
	
	public Nikon135mmF28(){ this(1.0); }
	
	public Nikon135mmF28(double scale) {
		super("Nikon135mmF28", new double[][] { //This is 100mm/2.8 from Nikkor patent 4,057,330 (1977), ex#2			 
					{  0, 		6.44, 		0, 		0,		0 		}, //move lens plane so that focal plane is at +135mm
					{  41.630,  5.407,   18.4,	1.62041, 60.3,	},	
					{ 300.296,  0.444,   18.4,	0,       0		},
					{  28.778,  9.630,   17.4,	1.62041, 60.3	},
					{  61.593,  2.074,   13.8,	0,       0		},
					{ 157.876,  3.481,   14.2,	1.78470, 26.1	},
					{ -58.637,  1.111,   14.2,	1.74000, 28.2	},
					{  19.616,  25.704,  10.8,	0,       0		},
					{  62.222,  1.630,   10.4,	1.72825, 28.3	},
					{ 186.713,  0,       10.4,	0,       0		} 
				}, 
				0.001 * 1.35 * scale, //mm --> m, and make it 135mm instead of 100mm
				0.030);		
	}

	public Nikon135mmF28(double[] centre) { this(centre, 1.0); }
	public Nikon135mmF28(double[] centre, double scale) { this(centre, scale, null); }
		
	public Nikon135mmF28(double[] centre, double scale, double axis[]) {
		this(scale);
		shift(centre);
		if(axis != null){
			double y[] = Util.createPerp(axis);
			double z[] = Util.reNorm(Util.cross(axis, y));
			rotate(centre, Algorithms.rotationMatrix(new double[][]{ axis, y, z }));
		}
	}
	
	public double getCaseRadius() { return 0.030; }
}
