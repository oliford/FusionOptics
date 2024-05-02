package fusionOptics;

import java.io.*;
import java.util.Properties;

import uk.co.oliford.jolu.SettingsManager;

/** Provides settings for paths and things
 * 
 * This is mostly copied from minerva so that minerva-optics can be indepedant of Minerva's core.
 * 
 */
public abstract class MinervaOpticsSettings {
	
	static{ new SettingsManager("minerva", true); }
	
	/* Really common settings requests are below: */
		
	/** @return apth for tests outputs (typically want volatile storage)
	 * Default is $TEMP/minerva/tests/
	 */	
	public static String getTestsOutputPath(){
		return SettingsManager.defaultGlobal().getPathProperty(
				"Tests.OutputPath", 
				System.getProperty("java.io.tmpdir") + "/minerva/tests/");
	}
	
	/** @return path for minera applications results (data, svgs etc)
	 * Default is $TEMP/minerva/results/
	 */
	public static String getAppsOutputPath(){
		return SettingsManager.defaultGlobal().getPathProperty(
				"minerva.apps.resultpath", 
				System.getProperty("java.io.tmpdir") + "/minerva/results/");
	}

	/** @return path to the root of the Minerva java source where all the project directories are
	 * placed. This is useful for guessing the path to data files stored inside the SVN. */
	public static String getMinervaSourcePath(){
		return SettingsManager.defaultGlobal().getPathProperty(
				"minerva.source.path", 
				System.getProperty("java.io.tmpdir") + "/minerva/source/");		
	}

	/** For operations that are sometimes rapid and sometime stake ages (typically webservices with caching layers)
	 * If the overall operation is taking over this long, start displaying status info.
	 *  'minerva.user.attentionSpanMilisecs'
	 * @return The user's attention span in ms */
	public static int getUserAttentionSpan(){
		return Integer.parseInt(SettingsManager.defaultGlobal().getProperty("minerva.user.attentionSpanMilisecs", "5000"));
	}	
}
