import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Properties;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class ClusteredMotifsMain {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		/* ------- Command line options --------*/
		System.out.println("Parsing commandline arguments");

		Options options = initializeOptions();

		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();

		try {
			CommandLine cmd = parser.parse(options, args);

			/* ----- help flag ----- */
			if (cmd.hasOption("h")) {
				formatter.printHelp("localEnrich.main", options);
			}

			/* required ARG */
			if (!cmd.hasOption("t")) {
				throw new MissingOptionException("The '-t' p-value threshold is required");
			}

			System.out.println("Loading parameters file");
			Properties params = new Properties();
			params.load(new FileInputStream(cmd.getOptionValue("p")));		

			/* max tolerated significance score */
			double pvalThreshold = Double.parseDouble(cmd.getOptionValue("t"));

			String wd = params.getProperty("working_directory");
			String annotationWD = params.getProperty("annotation_directory");
			String networkName = params.getProperty("project_name");

			int clusteringMeasure = Integer.parseInt(params.getProperty("clusteringMeasure", "0"));
			double percentThreshold = Double.parseDouble(params.getProperty("percentThreshold", "0.2"));

			String clusteringName = "";

			switch(clusteringMeasure) {
			case 0: clusteringName = "_TPD";
			break;
			case 1: clusteringName = "_TPPD_p" + percentThreshold;
			break;
			case 2: clusteringName = "_coreTPD_p" + percentThreshold;
			break;
			}

			String fastaFile = wd + params.getProperty("fastaFile");
			String motifClusteringPrefix = wd + "motifClustering/" + networkName + clusteringName + "_MotifClustering_";

			String significantMotifsFile = wd +  networkName + clusteringName +"_signficantMotifs_p" + pvalThreshold +".tsv";

			String annotationPrefixFile = annotationWD + "/motif_enumeration/annotations/annotation_";
			String protAnnotationFreqFile = wd + networkName + "_protFreqAnnotation.tsv";

			String extractedAnnotationsFile = wd + networkName + clusteringName + "_annotationSubset_p" + pvalThreshold +".tsv";
			String corePorteinsFile = wd + networkName + clusteringName +"_coreProteinsByMotif_p" + pvalThreshold + ".tsv";

			String proteinToRefSeqIdFile = wd + networkName + "_proteinsInNetwork_info.tsv";

			/* directories */
			File motifClusteringDir =  new File(wd + "motifClustering/");
			File annotationDir = new File(annotationWD + "/motif_enumeration/annotations/");

			/* ------- Get the step number  ---------------*/
			int step = 0;

			if(cmd.hasOption("step")) {
				step = Integer.parseInt(cmd.getOptionValue("step"));
			}

			String mode = "";
			String condition = "";
			if(cmd.hasOption("m")) {
				mode = cmd.getOptionValue("m");
				System.out.println("Running mode : " + mode);
				switch(mode) {
				case "all":
					condition = "all/";
					break;
				case "core":
					condition = "core/";
				}
			} else {
				/* required ARG */
				if(step !=1) {
					throw new MissingOptionException("The '-m' mode is required");
				}
			}
			System.out.println("Next step: " + step + "\n");
			switch(step) {
			case 1: 
				/* ----------------------------------------------------------------------------------------------------------- 
				 * Part 1 - Obtain significantly clustered motifs and calculate their similarity
				 * ----------------------------------------------------------------------------------------------------------- */
				System.out.println("Running Step 1: Cluster similar significantly clustered motifs");
				/* -- Identify motifs that pass significant threshold: (1) print details to separate file, (2) store motif and file # in map --*/ 
				System.out.println("+ Identifying significant motifs");
				IdentifyMotifs.getSignificantMotifs(motifClusteringPrefix, motifClusteringDir, pvalThreshold, significantMotifsFile);

				System.out.println("+ Loading significant motifs");
				HashMap<String, Integer> motifMapOfFileIdx = IdentifyMotifs.loadSignificantMotifs(significantMotifsFile);
				System.out.println("Number of significant motifs: " + motifMapOfFileIdx.size() + "\n");

				/* -- Get annotations (motif = protein1|protein2|...|proteinN) of significant motifs for easy access in future analysis -- */ 
				System.out.println("+ Identifying annotated proteins");
				IdentifyMotifs.getAnnotatedProteinInfo(motifMapOfFileIdx, protAnnotationFreqFile, annotationPrefixFile,	extractedAnnotationsFile, annotationDir);
				System.out.println("extracted annotations are stored: " + wd);

				/* -- For CoreTPD and TPPD; load annotation subset and determine core proteins; print to separate file -- */
				if(clusteringMeasure == 1 || clusteringMeasure == 2) {

					System.out.println("+ Identifying Core Proteins**");

					String distanceMatrixFile = wd + networkName + "_removedOverConnectedProteins_" + params.getProperty("maxInteractions") + "_distanceMatrix2.txt";

					IdentifyCoreProteins.getCoreProteins(extractedAnnotationsFile, Double.parseDouble(params.getProperty("percentThreshold")), 
							proteinToRefSeqIdFile, distanceMatrixFile, corePorteinsFile);
					System.out.println("extracted annotations are stored: " + wd);
				}

				/* -- Compute similarity between significant proteins for determination of motif family -- */
				createDir(wd + "/motifFamilies/");
				createDir(wd + "/motifFamilies/all/");

				String motifsInMatrixFile = wd +  "motifFamilies/all/" + networkName + clusteringName + "_p" + pvalThreshold + "_MotifsInMatrix.tsv";
				String similarityMatrix = wd + "motifFamilies/all/" + networkName + clusteringName + "_p" + pvalThreshold + "_DistanceMatrix.tsv" ;
				System.out.println("+ Loading annotation info");
				/* Search through annotation Files to get proteins annotated by significant motifs: (1) store in map for similarity measuring, (2) print to file for local testing */
				//HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.getAnnotatedProteinInfo(motifMapOfFileIdx, annotationPrefixFile);
				HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.loadAnnotatedProteinInfo(extractedAnnotationsFile);
				System.out.println("Found motif info: " + motifMapOfAnnotatedProteins.size() + "\n");


				System.out.println("+ Computing similarity");
				Similarity.computeMotifSimilary(motifMapOfAnnotatedProteins, motifsInMatrixFile, similarityMatrix);
				System.out.println("\nCreated similarity matrix for all proteins: " + similarityMatrix);

				/* For CoreTPD and TPPD; load annotation subset and determine core proteins; print */
				if(clusteringMeasure == 1 || clusteringMeasure == 2) {

					createDir(wd + "/motifFamilies/core/");

					motifsInMatrixFile = wd +  "motifFamilies/core/" + networkName + clusteringName + "_CoreProteins_p" + pvalThreshold + "_MotifsInMatrix.tsv";
					similarityMatrix = wd + "motifFamilies/core/" + networkName + clusteringName + "_CoreProteins_p" + pvalThreshold + "_DistanceMatrix.tsv" ;

					System.out.println("+ Similarity calculation for Core Proteins");
					motifMapOfAnnotatedProteins = IdentifyMotifs.loadAnnotatedProteinInfo(corePorteinsFile);
					System.out.println("Found motif info: " + motifMapOfAnnotatedProteins.size() + "\n");

					System.out.println("+ Computing similarity");
					Similarity.computeMotifSimilary(motifMapOfAnnotatedProteins, motifsInMatrixFile, similarityMatrix);
					System.out.println("\nCreated similarity matrix for core proteins: " + similarityMatrix);

				}
				System.out.println("--- Completed step 1 ---");
				System.out.println("Proceed to hierarchical clustering in R / Jupiter notebook");
				break;
			case 2:
				/* ----------------------------------------------------------------------------------------------------------- 
				 * Part 2 - Choose representative motif from similar clustered group
				 * ----------------------------------------------------------------------------------------------------------- */
				System.out.println("Running Step 2: Selecting representative motif for each family");

				/* required ARG */
				if (!cmd.hasOption("c")) {
					throw new MissingOptionException("The '-c' cutree height is required");
				}
				/* Assess motif families */ 
				String height = cmd.getOptionValue("c");

				createDir(wd + "/motifFamilies/" + condition + "pwm/");

				String motifFamilyFilePrefix = wd +  "motifFamilies/" + condition + "MotifFamily_h" + height +"_group";
				String motifInfoFile = wd +  "motifFamilies/" + condition + networkName + clusteringName + "_h" + height+ "_motifFamiliesInfo.tsv";

				File motifsDir = new File(wd + "motifFamilies/" + condition);
				System.out.println("+ Assessing motif families");
				MotifFamily.setMotifFamilyRepresentative(motifFamilyFilePrefix, significantMotifsFile, motifsDir, motifInfoFile, height);
				System.out.println("\nRepresentative motifs are found under : " + motifInfoFile);
				
				
				
				System.out.println("--- Completed step 2 ---");
				break;

			case 3:
				/* ----------------------------------------------------------------------------------------------------------- 
				 * Part 3 - Evaluate motif families print info for sequence logo generation
				 * ----------------------------------------------------------------------------------------------------------- */

				/* required ARG */
				if (!cmd.hasOption("c")) {
					throw new MissingOptionException("The '-c' cutree height is required");
				}

				/* required ARG */
				if (!cmd.hasOption("n")) {
					throw new MissingOptionException("The '-n' motif family number is required");
				}
				/* Assess motif families */ 
				height = cmd.getOptionValue("c");

				String motifInstancesPrefix = wd + "motifFamilies/" + condition + "pwm/" +  networkName + clusteringName + "_h" + height + "_motifInstances_motifFamily";
				String motifPPMPrefix = wd +  "motifFamilies/" +  condition + "pwm/" + networkName + clusteringName+ "_h" + height + "_ppm_motifFamilyGroup";
				motifInfoFile = wd +  "motifFamilies/" + condition + networkName + clusteringName + "_h" + height+ "_motifFamiliesInfo.tsv";
				MotifFamily.assessMotifFamilies(motifInfoFile, Integer.parseInt(cmd.getOptionValue("n")), proteinToRefSeqIdFile, motifInstancesPrefix, motifPPMPrefix, extractedAnnotationsFile, fastaFile);
				System.out.println("\nPWM are stored : " + wd + "motifFamilies/" + condition + "pwm/");

				
				String wd4r = wd + "motifFamilies/" + condition;
				String projectName = networkName + clusteringName;
				motifsDir = new File(wd + "motifFamilies/" + condition);

				String[] argsR = new String[6];
				argsR[0] = "/usr/local/bin/Rscript";
				argsR[1] = "seqLogo.R";
				argsR[2] = wd4r;
				argsR[3] = projectName;
				argsR[4] = String.valueOf(motifsDir.list().length);
				argsR[5] = String.valueOf(height);

				Runtime.getRuntime().exec(argsR);

				System.out.println("--- Completed step 3 ---");
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
			System.out.println("+ creating directory: " + path);
			dir.mkdir();
		}
	}

	private static Options initializeOptions() {

		Options options = new Options();

		options.addOption("p", "properties", true, "properties file");
		options.addOption("t", "threshold", true, "p-value threshold, range [0.0-1.0] "); 
		options.addOption("s", "step", true, "step to execute, options [1-3]");
		options.addOption("m", "mode", true, "values are 'all' or 'core', required for step 2 and 3");
		options.addOption("c", "cutree", true, "cutree height, required for step 3");
		options.addOption("n", "motifNumber", true, "number corresponding to motif family to evaluate, required for step 3");
		options.addOption("h", "help", false, "show help");

		return options;
	}



}
