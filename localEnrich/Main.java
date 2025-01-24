import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Properties;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import graph.Interaction;
import graph.Protein;
import sampling.ApproximateNormalDistribuiton;
import sampling.MotifSampling;
import sampling.ProteinAnnotations;
import utils.AnnotationCompanionFiles;
import utils.AssessEnrichment;
import utils.Calculator;
import utils.CorrelationGraphLoader;
import utils.DistanceMatrix;
import utils.MotifEnrichment;
import utils.NetworkProteins;

public class Main {

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

			/* ------ Import parameters ---------- */
			System.out.println("Loading properties file");
			Properties params = new Properties();
			params.load(new FileInputStream(cmd.getOptionValue("p")));	

			String wd = params.getProperty("working_directory");

			String networkName = params.getProperty("project_name");

			String projectName = "";

			if(Boolean.parseBoolean(params.getProperty("nullModel"))){
				projectName = networkName + "_nullModel";
			} else {
				projectName = networkName;
			}

			boolean removeOverlyConnectedProteins = Boolean.parseBoolean(params.getProperty("removeOverlyConnectedProteins"));
			int maxInteractions = Integer.parseInt(params.getProperty("maxInteractions"));

			String fastaFile = wd + params.getProperty("fastaFile");
			String mapProtToRefSeqFile = wd + params.getProperty("geneIdsFile");

			String correlationRepository = wd + params.getProperty("networkRepositoryFile");
			String proteinsInNetworkOutputFile = wd + networkName +"_listOfProteinNames.tsv";
			String protInfoFile = wd + networkName + "_proteinsInNetwork_info.tsv";

			String distanceMatrixFile = wd + networkName + "_distanceMatrix.txt"; 
			String distanceMatrix2File = wd +  networkName + "_distanceMatrix2.txt";

			if(removeOverlyConnectedProteins) {
				distanceMatrixFile = wd + networkName + "_removedOverConnectedProteins_" + maxInteractions + "_distanceMatrix.txt";
				distanceMatrix2File = wd + networkName + "_removedOverConnectedProteins_" + maxInteractions + "_distanceMatrix2.txt";
			}

			String annotationWD = params.getProperty("annotation_directory");
			String degenAnnotationPrefix = annotationWD + "/motif_enumeration/annotations/annotation_";

			String proteinAnnotationFrequencyFile = wd + projectName + "_protFreqAnnotation.tsv";

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

			String mcSamplingPrefix = wd + "mcDistributions/" + projectName + clusteringName + "_mcSamplingDistribution_";
			String normalDistributionParamsFile = wd + projectName + clusteringName +"_normalDistributionParams.tsv";

			String testedDegenMotifsOutputPrefix = wd + "motifClustering/" + projectName + clusteringName + "_MotifClustering_";
			String significanceScoresFile = wd + projectName + clusteringName + "_SignificanceScores.tsv";

			/* ----------------------------------------------------------------------------------------------------------- 
			 * Part 0 - Formatting the network - must be run every time -
			 * 
			 * + The correlation network is loaded, establishing proteins in the network and their corresponding interaction/relationships
			 * + The shortest paths are calculated and output as a distance matrix. Only the biggest connected component is kept
			 * 
			 * ----------------------------------------------------------------------------------------------------------- */
			System.out.println("\nRunning step 0 - initializing network");
			System.out.println("+ Loading interaction repository");
			ArrayList<Interaction> interactionList = CorrelationGraphLoader.loadGraphFromCorrelationNetwork(correlationRepository, fastaFile, 
					mapProtToRefSeqFile, proteinsInNetworkOutputFile, Double.parseDouble(params.getProperty("corrThreshold", "0.0")), 
					removeOverlyConnectedProteins, maxInteractions);
			System.out.println("Number of interactions:" + interactionList.size() + "\n");

			System.out.println("+ Getting list of proteins in network");
			ArrayList<Protein> proteinList = NetworkProteins.getProteinsInNetwork(interactionList);
			System.out.println("Number of Proteins: " + proteinList.size() + "\n");

			File f = new File(distanceMatrixFile);
			if(!f.exists() && !f.isDirectory()) {
				System.out.println("+ Generating distance matrix");
				DistanceMatrix.computeDistanceMatrix(interactionList, proteinList, distanceMatrixFile);
			}

			System.out.println("+ Loading distance matrix");
			double[][] distanceMatrix = DistanceMatrix.loadDistanceMatrix(distanceMatrixFile, proteinList); 

			/* Determine which proteins are disconnected*/ 
			boolean[] proteinsToKeep = Calculator.determineConnectedProteins(distanceMatrix);
			/* Update proteins to keep in the network (ie. those that contribute to the most connected component) */
			ArrayList<Protein> proteinList2 = NetworkProteins.modifyNetworkProteinsList(proteinList, proteinsToKeep, protInfoFile);
			HashSet<String> proteinSet = NetworkProteins.getProteinSet(proteinList2);

			if(proteinList.size() != proteinList2.size()) {
				//System.out.println("+ Checking for disconnected components");
				File f1 = new File(distanceMatrix2File);
				if(!f1.exists() && !f1.isDirectory()) {
					/* Update distance matrix to only include proteins that contribute to the most connected component of the graph */
					System.out.println("+ Updating distance matrix without disconnected components");
					DistanceMatrix.updateDistanceMatrix(proteinsToKeep, distanceMatrix, distanceMatrix2File);
				}

				/* Load distance matrix representing fully connected component */
				System.out.println("+ Loading updated distance matrix");
				distanceMatrix = DistanceMatrix.loadDistanceMatrix(distanceMatrix2File, proteinList2);
			}
			System.out.println("--- Completed step 0 ---");

			int lowerBound = Integer.parseInt(params.getProperty("lowerBoundToSample", "3"));
			int upperBound = Integer.parseInt(params.getProperty("upperBoundToSample", "2000"));
			int numOfSamplings = Integer.parseInt(params.getProperty("numberOfSamplings", "100000"));

			String protFreqFilePrefix = annotationWD + "/ProtAnnotationFreq_" + projectName + "_n" + lowerBound + "_" + upperBound + "/" 
					+"/" + projectName + "_n" + lowerBound + "_" + upperBound + "_proteinAnnotationFreq_" ;
			String annotationCompanionFilePrefix = annotationWD + "/CompanionFiles_" + projectName + "_n" + lowerBound + "_" + upperBound + "/" 
					+ projectName + "_n" + lowerBound + "_" + upperBound + "_motifAnnotationsCompanionFile_";

			/* ------- Get the step number  ---------------*/
			int step = 0;
			int fileNumber = 0;

			if(cmd.hasOption("step")) {
				step = Integer.parseInt(cmd.getOptionValue("step"));
			}
			System.out.println("Next step: " + step + "\n");
			switch(step) {
			case 1: 
				/* ----------------------------------------------------------------------------------------------------------- 
				 * Part 1 - Create companion files for Monte Carlo Sampling and assessing the clustering of motifs
				 *
				 * This section aims to identify any discrepancies between the annotation file and the final network. 
				 * It serves to confirm the size of annotations for reference during the motif clustering step. 
				 * ----------------------------------------------------------------------------------------------------------- */
				System.out.println("\n\nRunning step 1 - creating accessory files");

				/* required ARG */
				if (!cmd.hasOption("a")) {
					throw new MissingOptionException("The '-a' annotation file number is required for this step");
				}

				/* file number to test */
				fileNumber = Integer.parseInt(cmd.getOptionValue("a"));

				/* required directories */
				String f3 = annotationWD + "/CompanionFiles_" + projectName + "_n" + lowerBound + "_" + upperBound + "/";
				String f4 = annotationWD + "/ProtAnnotationFreq_" + projectName + "_n" + lowerBound + "_" + upperBound + "/"; 
				createDir(f3);
				createDir(f4);

				System.out.println("+ creating companion file: " + fileNumber);
				System.out.println("Companion files are stored under : " + annotationWD + "/CompanionFiles_" + projectName + "_n" + lowerBound + "_" + upperBound + "/");

				AnnotationCompanionFiles.assessAnnotationFile(degenAnnotationPrefix, annotationCompanionFilePrefix, proteinSet, 
						fileNumber, fileNumber, lowerBound, upperBound);

				System.out.println("+ enumerating protein annotation frequency file: " + fileNumber);
				System.out.println("Annotation frequency files are stored under : " + annotationWD + "/ProtAnnotationFreq_" + projectName + "_n" + lowerBound + "_" + upperBound + "/");

				ProteinAnnotations freq = new ProteinAnnotations(lowerBound, upperBound, proteinSet);
				freq.calculateProteinAnnotationFrequency(degenAnnotationPrefix, annotationCompanionFilePrefix, fileNumber, fileNumber, protFreqFilePrefix);

				System.out.println("--- Completed step 1 ---");

				break;
			case 2:
				/* ----------------------------------------------------------------------------------------------------------- 
				 * Part 2 - Combine protein annotation frequency 
				 * ----------------------------------------------------------------------------------------------------------- */
				System.out.println("\nRunning step 2 - combining annotation frequencies");
				freq = new ProteinAnnotations(lowerBound, upperBound, proteinSet);
				File directory4 = new File(annotationWD + "/ProtAnnotationFreq_" + projectName + "_n" + lowerBound + "_" + upperBound + "/"); 

				System.out.println("+ Combining protein annotation frequency data");
				System.out.println("Creating annotation frequency file :" +  proteinAnnotationFrequencyFile);

				freq.combineProteinFrequencyData(protFreqFilePrefix, directory4, proteinAnnotationFrequencyFile);

				System.out.println("--- Completed step 2 ---");
				break;
			case 3:
				/* ----------------------------------------------------------------------------------------------------------- 
				 * Part 3 - Perform Monte Carlo Sampling 
				 * ----------------------------------------------------------------------------------------------------------- */
				System.out.println("\nRunning step 3 - Monte Carlo Sampling procedure");
				createDir(wd + "/mcDistributions");

				/* required ARG */
				if (!cmd.hasOption("d")) {
					throw new MissingOptionException("The '-d' distribution sample size is required for this step");
				}
				int sampleSize = Integer.parseInt(cmd.getOptionValue("d"));

				System.out.println("+ Performing Monte Carlo Sampling Procedure : " + sampleSize);
				System.out.println("Monte Carlo Sampling distributions are stored under : " + wd + "/mcDistributions/");

				MotifSampling sampling = new MotifSampling(proteinAnnotationFrequencyFile, proteinList2, distanceMatrix, clusteringMeasure, percentThreshold); // 2 - Initialize sampling
				sampling.computeMultipleDistributions(sampleSize, sampleSize, numOfSamplings, mcSamplingPrefix); // 3 - Perform sampling for n proteins

				System.out.println("--- Completed step 3 ---");
				break;
			case 4:
				/* ----------------------------------------------------------------------------------------------------------- 
				 * Part 4 - Calculate the mean and standard deviation for each cluster size
				 * ----------------------------------------------------------------------------------------------------------- */
				System.out.println("\nRunning step 4 - Calculate normal distribution parameters");
				System.out.println("Creating normal distribution parameters file : "+ normalDistributionParamsFile);

				ApproximateNormalDistribuiton.getNormalDistributionParams(mcSamplingPrefix, lowerBound, upperBound, numOfSamplings, normalDistributionParamsFile);
				System.out.println("--- Completed step 4 ---");
				break;
			case 5:
				/* ----------------------------------------------------------------------------------------------------------- 
				 * Part 5 - Calculate the clustering and clustering significance of each motif
				 * ----------------------------------------------------------------------------------------------------------- */

				System.out.println("\nRunning step 5 - Assess motif clustering and signifcance");
				createDir(wd + "/motifClustering");

				/* required ARG */
				if (!cmd.hasOption("a")) {
					throw new MissingOptionException("The '-a' annotation file number is required for this step");
				}

				/* file number to test */
				fileNumber = Integer.parseInt(cmd.getOptionValue("a"));

				System.out.println(" + perform motif clustering: " + fileNumber);
				System.out.println("Clustering details are stored under: " + wd + "/motifClustering");

				MotifEnrichment m = new MotifEnrichment(distanceMatrix, proteinList2, normalDistributionParamsFile, lowerBound, upperBound,clusteringMeasure, percentThreshold);
				m.testMotifClustering(degenAnnotationPrefix, annotationCompanionFilePrefix, testedDegenMotifsOutputPrefix, fileNumber, fileNumber);

				System.out.println("--- Completed step 5 ---");
				break;
			case 6: 
				/* ----------------------------------------------------------------------------------------------------------- 
				 * Part 6 - Print all significance scores for FDR estimation
				 * ----------------------------------------------------------------------------------------------------------- */
				System.out.println("\nRunning step 6 - combine significance scores for FDR estimation");
				File dir = new File(wd + "motifClustering/");
				AssessEnrichment.assessSignificanceScores(testedDegenMotifsOutputPrefix, dir, significanceScoresFile);
				System.out.println("--- Completed step 6 ---");
				break;
			}
		} catch (ParseException e) {
			System.out.println("Error parsing arguments: " + e.getMessage());
			formatter.printHelp("MapMotifsToProteins", options);
		}

	}


	public static void printRefSeqIdsInNetwork(String outputFile, ArrayList<Protein> proteinList) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(Protein prot: proteinList) {

				if(prot.getProteinId() != null) {

					for(String rnaId: prot.getProteinId()) {
						out.write(rnaId + "\n");
						out.flush();
					}
				}

			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void printProtAndRefSeqIdsInNetwork(String outputFile, ArrayList<Protein> proteinList) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(Protein prot: proteinList) {

				if(prot.getProteinId() != null) {
					out.write(prot.getProteinName() + "\t");

					for(String rnaId: prot.getProteinId()) {
						out.write(rnaId + "|");
					}
					out.write("\n");
					out.flush();
				}

			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
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
		options.addOption("s", "step", true, "step to execute, range [1-6]");
		options.addOption("a", "annotation", true, "annotation file number, required for steps [2,3,6], range [0-999]");
		options.addOption("d", "distribution", true, "sample size for Monte Carlo distribution, required for step 4, recommended range [20-2000]");
		options.addOption("h", "help", false, "show help");

		return options;
	}


}
