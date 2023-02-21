package precompileMotifs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import motifEnumeration.MotifDegeneration;


public class precompileMain {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		//		System.out.println("Loading parameters file \n");
		//		Properties params = new Properties();
		//		params.load(new FileInputStream(args[0]));
		//
		//		/* Command line arguments */ 
		//		String motifsToTestFile = args[1];
		//		String outputFile = args[2];

		//boolean enumerateDegenMotifs = Boolean.parseBoolean(params.getProperty("enumerateDegenMotifs"));
		//boolean enumeratePossibleDegenMotifs = Boolean.parseBoolean(params.getProperty("enumeratePossibleDegenMotifs"));
		//int motifLength = Integer.parseInt(params.getProperty("motifLength"));
		//int maxDegenThreshold = Integer.parseInt(params.getProperty("maxDegenThreshold"));

		int motifLength = 8;
		int maxDegenThreshold = 7;

		//		String wd = params.getProperty("working_directory");
		String wd = args[0];

		File dir = new File(wd + "/motif_enumeration/");
		if(! dir.exists()) {
			dir.mkdir();
		}

		File dir2 = new File(wd + "/motif_enumeration/mapDegendMotifs");
		if (! dir2.exists()){
			dir2.mkdir();
		}

		File dir3 = new File(wd + "/motif_enumeration/degenMotifSet");
		if (! dir3.exists()){
			dir3.mkdir();
		}

		//		String degenMotifSetPrefix = wd + params.getProperty("degenMotifsToTestPrefix").replace("\\s+", "");
		String degenMotifSetPrefix = wd + "/motif_enumeration/degenMotifSet/degenMotifs_";
		//		String motifSetToTest = wd + motifsToTestFile; // if enumerating degen motifs = list of non degen motifs, if mapping degen motifs = list of possible degen motifs
		String mapOfDegenMotifs = wd + "/motif_enumeration/mapDegendMotifs/mapDegenMotifs_";

		// Enumerating possible degen motifs only needs to be run once; output all possibilities for motifs of length l, with an alphabet of x characters. 
		//if(enumeratePossibleDegenMotifs) {
		
		if(dir3.list().length < 1000) {
			System.out.println("**Generating all possible degen motifs**");
			MotifDegeneration d1 = new MotifDegeneration(motifLength, maxDegenThreshold);
			d1.generateAllPossibleMotifs(degenMotifSetPrefix);
		}
		
		//}

		// MOTIF DEGENERATION SHOULD BE RUN REMOTELY //
		//		if(enumerateDegenMotifs) { 
		System.out.println("**Enumerating motifs from degenerate motifs**");
		MotifDegeneration d = new MotifDegeneration(motifLength, maxDegenThreshold);
		d.enumerateNonDegenerateMotifs(degenMotifSetPrefix, Integer.parseInt(args[1]), mapOfDegenMotifs);
		//		}

	}

}
