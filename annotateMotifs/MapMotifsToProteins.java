import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Properties;

import org.apache.commons.cli.*;

public class MapMotifsToProteins {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		/* Command line options */
		Options options = new Options();

		options.addOption("p", "properties", true, "properties file");
		options.addOption("s", "step", true, "step to execute");
		options.addOption("n", "file_number", true, "number of file to map");
		options.addOption("h", "help", false, "show help");

		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();

		CommandLine cmd = null;

		try {
			cmd = parser.parse(options, args);

			if (cmd.hasOption("h")) {
				formatter.printHelp("MapMotifsToProteins", options);
			}

			System.out.println("Loading parameters file\n");
			Properties params = new Properties();
			params.load(new FileInputStream(cmd.getOptionValue("p")));	

			String wd = params.getProperty("annotation_directory");

			createDir(wd + "/motif_enumeration/");
			createDir(wd + "/motif_enumeration/annotations/");
			createDir(wd + "/motif_enumeration/motifs/");

			String motifPrefixFile = wd + "/motif_enumeration/motifs/motifs_";

			// Get the step number
			int step = Integer.parseInt(cmd.getOptionValue("step"));

			switch(step) {
			case 1: 

				int motifLength = 8;
				int maxDegenThreshold = 7;

				System.out.println("Running Step 1: Enumerate motifs");
				System.out.println("Enumerated motifs will be stored under: " + wd + "motif_enumeration/motifs/ \n");
				
				EnumerateMotifs d1 = new EnumerateMotifs(motifLength, maxDegenThreshold);
				d1.generateAllPossibleMotifs(motifPrefixFile);
				break;

			case 2: 
				if (!cmd.hasOption("n")) {
					throw new MissingOptionException("The '--file_number' option is required for the motif mapping step");
				}

				String fastaFile = wd + params.getProperty("fastaFile");
				String motifDir = wd + "/motif_enumeration/degenMotifSet/";
				String ids = wd + params.getProperty("geneIdsFile");

				String annotationFile = wd + "/motif_enumeration/annotations/annotation_";

//				String f = cmd.getOptionValue("file_number");
//				int fileNumber = Integer.parseInt(cmd.getOptionValue("f"));
				
				System.out.println("Running Step 2: Generating annotation file - " + cmd.getOptionValue("n"));
				System.out.println("Annotation files will be stored under: " + wd + "motif_enumeration/annotations/ \n");
				
				MotifMapper m = new MotifMapper(ids, annotationFile, fastaFile, motifDir, motifPrefixFile);
				m.mapMotifsToTheirAssociatedProteins(Integer.parseInt(cmd.getOptionValue("n")));

				break;
			}
		} catch (ParseException e) {
			System.out.println("Error parsing arguments: " + e.getMessage());
			formatter.printHelp("MapMotifsToProteins", options);
		}
	}

	private static void createDir(String path) {
		File dir = new File(path);
		if (!dir.exists()) {
			dir.mkdir();
		}
	}
}
