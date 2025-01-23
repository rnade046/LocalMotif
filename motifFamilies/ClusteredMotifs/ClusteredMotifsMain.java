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

			String proteinToRefSeqIdFile = wd +  networkName + "_proteinsInNetwork_info.tsv";

			/* directories */
			File motifClusteringDir =  new File(wd + "motifClustering/");
			File annotationDir = new File(annotationWD + "/motif_enumeration/annotations/");

			/* ------- Get the step number  ---------------*/
			int step = 0;

			if(cmd.hasOption("step")) {
				step = Integer.parseInt(cmd.getOptionValue("step"));
			}
			System.out.println("Next step: " + step + "\n");
			switch(step) {
			case 1: 
				/* ----------------------------------------------------------------------------------------------------------- 
				 * Part 1 - Obtain significantly clustered motifs and calculate their similarity
				 * ----------------------------------------------------------------------------------------------------------- */

				/* -- Identify motifs that pass significant threshold: (1) print details to separate file, (2) store motif and file # in map --*/ 
				System.out.println("**Identifying significant motifs**");
				IdentifyMotifs.getSignificantMotifs(motifClusteringPrefix, motifClusteringDir, pvalThreshold, significantMotifsFile);

				System.out.println("**Loading significant motifs**");
				HashMap<String, Integer> motifMapOfFileIdx = IdentifyMotifs.loadSignificantMotifs(significantMotifsFile);
				System.out.println("Number of significant motifs: " + motifMapOfFileIdx.size() + "\n");

				/* -- Get annotations (motif = protein1|protein2|...|proteinN) of significant motifs for easy access in future analysis -- */ 
				System.out.println("**Identifying annotated proteins**");
				IdentifyMotifs.getAnnotatedProteinInfo(motifMapOfFileIdx, protAnnotationFreqFile, annotationPrefixFile,	extractedAnnotationsFile, annotationDir);

				/* -- For CoreTPD and TPPD; load annotation subset and determine core proteins; print to separate file -- */
				if(clusteringMeasure == 1 || clusteringMeasure == 2) {

					System.out.println("**Identifying Core Proteins**");

					String distanceMatrixFile = wd + networkName + "_removedOverConnectedProteins_" + params.getProperty("maxInteractions") + "_distanceMatrix2.txt";

					IdentifyCoreProteins.getCoreProteins(extractedAnnotationsFile, Double.parseDouble(params.getProperty("percentThreshold")), 
							proteinToRefSeqIdFile, distanceMatrixFile, corePorteinsFile);
				}

				/* -- Compute similarity between significant proteins for determination of motif family -- */
				File directory2 = new File(wd + "/motifFamilies/"); 
				if (! directory2.exists()){
					System.out.println("creating directory: motifFamilies/");
					directory2.mkdir();
				}

				String motifsInMatrixFile = wd +  "motifFamilies/" + networkName + clusteringName + "_p" + pvalThreshold + "_MotifsInMatrix.tsv";
				String similarityMatrix = wd + "motifFamilies/"+ networkName + clusteringName + "_p" + pvalThreshold + "_DistanceMatrix.tsv" ;
				System.out.println("**Loading annotation info**");
				/* Search through annotation Files to get proteins annotated by significant motifs: (1) store in map for similarity measuring, (2) print to file for local testing */
				//HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.getAnnotatedProteinInfo(motifMapOfFileIdx, annotationPrefixFile);
				HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.loadAnnotatedProteinInfo(extractedAnnotationsFile);
				System.out.println("Found motif info: " + motifMapOfAnnotatedProteins.size() + "\n");


				System.out.println("**Computing similarity**");
				Similarity.computeMotifSimilary(motifMapOfAnnotatedProteins, motifsInMatrixFile, similarityMatrix);
				
				/* For CoreTPD and TPPD; load annotation subset and determine core proteins; print */
				if(clusteringMeasure == 1 || clusteringMeasure == 2) {

					motifsInMatrixFile = wd +  "motifFamilies/" + networkName + clusteringName + "_CoreProteins_p" + pvalThreshold + "_MotifsInMatrix.tsv";
					similarityMatrix = wd + "motifFamilies/" + networkName + clusteringName + "_CoreProteins_p" + pvalThreshold + "_DistanceMatrix.tsv" ;

					System.out.println("**Similarity calculation for Core Proteins**");
					motifMapOfAnnotatedProteins = IdentifyMotifs.loadAnnotatedProteinInfo(corePorteinsFile);
					System.out.println("Found motif info: " + motifMapOfAnnotatedProteins.size() + "\n");

					System.out.println("**Computing similarity**");
					Similarity.computeMotifSimilary(motifMapOfAnnotatedProteins, motifsInMatrixFile, similarityMatrix);
				}
				break;
			case 3:
				/* ----------------------------------------------------------------------------------------------------------- 
				 * Part 3 - Evaluate motif families : choose representative motif and print info for sequence logo generation
				 * ----------------------------------------------------------------------------------------------------------- */
				
				/* required ARG */
				if (!cmd.hasOption("c")) {
					throw new MissingOptionException("The '-c' cutree height is required");
				}
				/* Assess motif families */ 
				String height = cmd.getOptionValue("c");
				String condition = "Groups_h" + height + "/";

				String motifInstancesPrefix = wd + "motifFamilies/" + condition + networkName + clusteringName + "_h" + height + "_motifInstances_motifFamily";
				String motifPPMPrefix = wd +  "motifFamilies/" +  condition + networkName + clusteringName+ "_h" + height + "_ppm_motifFamilyGroup";
				String motifInfoFile = wd +  "motifFamilies/" + condition + networkName + clusteringName + "_h" + height+ "_motifFamiliesInfo.tsv";

				String motifFamilyFilePrefix = wd +  "motifFamilies/" + condition + "MotifFamily_h" + height +"_group";

				File motifsDir = new File(wd + "motifFamilies/" + condition);
				System.out.println("**Assessing motif families**");
				MotifFamily.assessMotifFamilies(motifFamilyFilePrefix, motifsDir, significantMotifsFile, proteinToRefSeqIdFile, motifInstancesPrefix, motifPPMPrefix, motifInfoFile, extractedAnnotationsFile, fastaFile);
				
				String wd4r = wd + "motifFamilies/" + condition;
				String projectName = networkName + clusteringName;

				String[] argsR = new String[6];
				argsR[0] = "/usr/local/bin/Rscript";
				argsR[1] = "seqLogo.R";
				argsR[2] = wd4r;
				argsR[3] = projectName;
				argsR[4] = String.valueOf(motifsDir.list().length);
				argsR[5] = String.valueOf(height);

				Runtime.getRuntime().exec(argsR);

				if(clusteringMeasure == 1 || clusteringMeasure == 2) {
					height = params.getProperty("coreHeight", "0");
					double height2 = Double.parseDouble(params.getProperty("coreHeight", "0"));
					condition = "Groups_CoreProteins_h" + height + "/";

					motifInstancesPrefix = wd + "motifFamilies/" + condition + networkName + clusteringName + "_coreProteins_h" + height2 + "_motifInstances_motifFamily";
					motifPPMPrefix = wd +  "motifFamilies/" + condition + networkName + clusteringName + "_coreProteins_h" + height2 + "_ppm_motifFamilyGroup";
					motifInfoFile = wd +  "motifFamilies/" + condition + networkName + clusteringName + "_coreProteins_h" + height2 + "_motifFamiliesInfo.tsv";

					motifFamilyFilePrefix = wd +  "motifFamilies/" + condition + networkName + clusteringName + "_CoreProteins_p" + pvalThreshold + "_ward.D2_h"+ height2 +"_group";
					
					motifsDir = new File(wd + "motifFamilies/" + condition);
					System.out.println("**Assessing motif instances for Core Proteins**");
					MotifFamily.assessMotifFamilies(motifFamilyFilePrefix, motifsDir, significantMotifsFile, proteinToRefSeqIdFile, motifInstancesPrefix, motifPPMPrefix, motifInfoFile, extractedAnnotationsFile, fastaFile);

					/* Generate */ 
					wd4r = wd + "motifFamilies/" + condition;
					projectName = networkName + clusteringName + "_CoreProteins";

					argsR[0] = "/usr/local/bin/Rscript";
					argsR[1] = "seqLogo.R";
					argsR[2] = wd4r;
					argsR[3] = projectName;
					argsR[4] = String.valueOf(motifsDir.list().length);
					argsR[5] = String.valueOf(height);

					Runtime.getRuntime().exec(argsR);
				}
				break;
			}
		} catch (ParseException e) {
			System.out.println("Error parsing arguments: " + e.getMessage());
			formatter.printHelp("MapMotifsToProteins", options);
		}
	}

	private static Options initializeOptions() {

		Options options = new Options();

		options.addOption("p", "properties", true, "properties file");
		options.addOption("s", "step", true, "step to execute, options: 1 or 3");
		options.addOption("t", "threshold", true, "p-value threshold");
		options.addOption("c", "cutree", true, "cutree height, required for step 3");
		options.addOption("h", "help", false, "show help");

		return options;
	}



}
