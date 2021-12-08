import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Properties;

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
import utils.RandomizeProteinAnnotations;

public class Main {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		System.out.println("**Loading parameters file** \n");
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));		

		String wd = params.getProperty("working_directory");
		String networkName = params.getProperty("network_name");

		String projectName = "";

		if(Boolean.parseBoolean(params.getProperty("nullModel"))){
			projectName = networkName + "_nullModel";
		} else {
			projectName = networkName;
		}

		boolean removeOverlyConnectedProteins = Boolean.parseBoolean(params.getProperty("removeOverlyConnectedProteins"));
		int maxInteractions = Integer.parseInt(params.getProperty("maxInteractions"));

		String fastaFile = wd + params.getProperty("fastaFile");
		String mapProtToRefSeqFile = wd + params.getProperty("mapGeneSymbolsToRefSeqIds");

		String correlationRepository = wd + params.getProperty("networkRepositoryFile");
		String proteinsInNetworkOutputFile = wd + networkName +"_listOfProteinNames.tsv";
		String protInfoFile = wd + networkName + "_proteinsInNetwork_info.tsv";

		String distanceMatrixFile = wd + networkName + "_distanceMatrix.txt"; 
		String distanceMatrix2File = wd +  networkName + "_distanceMatrix2.txt";

		if(removeOverlyConnectedProteins) {
			distanceMatrixFile = wd + networkName + "_removedOverConnectedProteins_" + maxInteractions + "_distanceMatrix.txt";
			distanceMatrix2File = wd + networkName + "_removedOverConnectedProteins_" + maxInteractions + "_distanceMatrix2.txt";
		}

		String degenAnnotationPrefix = params.getProperty("degenAnnotationPrefix") + "corrNetTop2_degenMotifMappedToProteinsInNetwork_";
		if(Boolean.parseBoolean(params.getProperty("nullModel"))){
			degenAnnotationPrefix = params.getProperty("degenAnnotationPrefix") + "annotation_nullModel_";
		}

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

		String testedDegenMotifsOutputPrefix = wd + "motifClustering/" + projectName + clusteringName + "_testedDegenMotifClustering_";
		String significanceScoresFile = wd + projectName + clusteringName + "_listOfCalculatedSignificanceScores.tsv";

		File directory = new File(wd + "/mcDistributions"); 
		if (! directory.exists()){
			System.out.println("creating directory: mcDistributions/");
			directory.mkdir();
		}

		File directory2 = new File(wd + "/motifClustering");
		if (! directory2.exists()){
			System.out.println("creating directory: motifClustering/ \n");
			directory2.mkdir();
		}

		System.out.println("**Loading interaction repository**");
		ArrayList<Interaction> interactionList = CorrelationGraphLoader.loadGraphFromCorrelationNetwork(correlationRepository, fastaFile, 
				mapProtToRefSeqFile, proteinsInNetworkOutputFile, Double.parseDouble(params.getProperty("corrThreshold")), 
				removeOverlyConnectedProteins, maxInteractions);
		System.out.println("Number of interactions:" + interactionList.size() + "\n");

		System.out.println("**Getting list of proteins in network**");
		ArrayList<Protein> proteinList = NetworkProteins.getProteinsInNetwork(interactionList);
		System.out.println("Number of Proteins: " + proteinList.size() + "\n");
		//printProtAndRefSeqIdsInNetwork(rnaIdListFile, proteinList);
		//printRefSeqIdsInNetwork(rnaIdListFile, proteinList);

		//CheckDegreeDistributions.assessDegreeDistribution(proteinList2, interactionList, (wd + projectName + "_degreesInNetwork.tsv"));

		File f = new File(distanceMatrixFile);
		if(!f.exists() && !f.isDirectory()) {
			System.out.println("**Generating distance matrix**");
			DistanceMatrix.computeDistanceMatrix(interactionList, proteinList, distanceMatrixFile);
		}

		/* Perform motif enumeration around here */

		System.out.println("**Loading distance matrix**");
		double[][] distanceMatrix = DistanceMatrix.loadDistanceMatrix(distanceMatrixFile, proteinList); 

		/* Determine which proteins are disconnected*/ 
		boolean[] proteinsToKeep = Calculator.determineConnectedProteins(distanceMatrix);
		/* Update proteins to keep in the network (ie. those that contribute to the most connected component) */
		ArrayList<Protein> proteinList2 = NetworkProteins.modifyNetworkProteinsList(proteinList, proteinsToKeep, protInfoFile);
		HashSet<String> proteinSet = NetworkProteins.getProteinSet(proteinList2);

		if(proteinList.size() != proteinList2.size()) {
			System.out.println("**Checking for disconnected components**");
			File f1 = new File(distanceMatrix2File);
			if(!f1.exists() && !f1.isDirectory()) {
				/* Update distance matrix to only include proteins that contribute to the most connected component of the graph */
				System.out.println("**Updating distance matrix**");
				DistanceMatrix.updateDistanceMatrix(proteinsToKeep, distanceMatrix, distanceMatrix2File);
			}

			/* Load distance matrix representing fully connected component */
			System.out.println("**Loading updated distance matrix**\n");
			distanceMatrix = DistanceMatrix.loadDistanceMatrix(distanceMatrix2File, proteinList2);
		} 
		
		/* Print network for MCL algorithm */
		if(Boolean.parseBoolean(params.getProperty("mcl"))) {
			String networkOutputFile = wd + networkName + "_mclNetwork.txt";
			
			System.out.println("**Formatting network for MCL analysis**");
			opt.FormatNetworkForMCL.formatMCLnetwork(proteinList2, interactionList, networkOutputFile);
		}
		
		if(Boolean.parseBoolean(params.getProperty("generateNullModel"))) {
			/* Output protein list and 3'UTRs */ 
			String nullMapFile = wd + projectName + "_nullModel_ProteinToRefseqIds.tsv";
			File f2 = new File(nullMapFile);
			if(!f2.exists() && !f2.isDirectory()) {
				System.out.println("**Randomizing protein-3'UTR associaions**");
				RandomizeProteinAnnotations.generateNullModel(proteinList2, nullMapFile);
			}
			
			String refSeqIdToMotifsFile = wd + "enumeratedMotifsPerRefSeqId.tsv";
			String outputProteinFile = wd + projectName + "_proteinToMotifs.tsv";
			File f3 = new File(outputProteinFile);
			if(!f3.exists() && !f3.isDirectory()) {
				System.out.println("**Generating map of proteins to motifs for annotation files**");
				RandomizeProteinAnnotations.generateProteinToMotifMap(nullMapFile, refSeqIdToMotifsFile, outputProteinFile);
			}
		}

		int lowerBound = Integer.parseInt(params.getProperty("lowerBoundToSample", "3"));
		int upperBound = Integer.parseInt(params.getProperty("upperBoundToSample", "2000"));
		int numOfSamplings = Integer.parseInt(params.getProperty("numberOfSamplings"));

		/* Check annotation files and create companion file */

		File directory3 = new File(params.getProperty("degenAnnotationPrefix") + "/CompanionFiles_" + projectName + "_n" + lowerBound + "_" + upperBound + "/"); 

		String annotationCompanionFilePrefix = params.getProperty("degenAnnotationPrefix") + "/CompanionFiles_" + projectName + "_n" + lowerBound + "_" + upperBound + "/" 
				+ projectName + "_n" + lowerBound + "_" + upperBound + "_motifAnnotationsCompanionFile_";

		if (! directory3.exists()){
			System.out.println("creating directory: /CompanionFiles_" + projectName + "_n" + lowerBound + "_" + upperBound);
			directory3.mkdir();
		}

		//if(directory3.list().length < Integer.parseInt(params.getProperty("numDegenMotifFiles"))) {
		if(Boolean.parseBoolean(params.getProperty("generateCompanionFiles"))) {

			System.out.println("Creating companion files");
			AnnotationCompanionFiles.assessAnnotationFile(degenAnnotationPrefix, annotationCompanionFilePrefix, proteinSet, 
					Integer.parseInt(args[1]), Integer.parseInt(args[2]), lowerBound, upperBound);
		}

		File directory4 = new File(params.getProperty("degenAnnotationPrefix") + "/ProtAnnotationFreq_" + projectName + "_n" + lowerBound + "_" + upperBound + "/"); 
		String protFreqFilePrefix = params.getProperty("degenAnnotationPrefix") + "/ProtAnnotationFreq_" + projectName + "_n" + lowerBound + "_" + upperBound + "/" 
				+"/" + projectName + "_n" + lowerBound + "_" + upperBound + "_proteinAnnotationFreq_" ;

		if (! directory4.exists()){
			System.out.println("creating directory: /ProtAnnotationFreq_" + projectName + "_n" + lowerBound + "_" + upperBound);
			directory4.mkdir();
		}

		if(Boolean.parseBoolean(params.getProperty("calculateProteinAnnotationFreq"))) {	
			// For MC sampling 
			// 1 - Make list: protein = #motifs (degen + non degen) from full annotation list	>> Do this once

			System.out.println("**Enumerating protein annotation frequency files**");
			ProteinAnnotations freq = new ProteinAnnotations(lowerBound, upperBound, proteinSet);
			freq.calculateProteinAnnotationFrequency(degenAnnotationPrefix, annotationCompanionFilePrefix, Integer.parseInt(args[1]), Integer.parseInt(args[2]), protFreqFilePrefix);

			if(directory4.list().length == Integer.parseInt(params.getProperty("numDegenMotifFiles"))) {
				System.out.println("Combining protein annotation frequency data /n");
				freq.combineProteinFrequencyData(protFreqFilePrefix, Integer.parseInt(params.getProperty("numDegenMotifFiles")), proteinAnnotationFrequencyFile);
			}
		}

		/* Perform Monte Carlo Sampling procedure */
		if(Boolean.parseBoolean(params.getProperty("performMCprocedure"))) {
			System.out.println("**Performing Monte Carlo Sampling Procedure**");
			MotifSampling sampling = new MotifSampling(proteinAnnotationFrequencyFile, proteinList2, distanceMatrix, clusteringMeasure, percentThreshold); // 2 - Initialize sampling
			sampling.computeMultipleDistributions(Integer.parseInt(args[1]), Integer.parseInt(args[2]), numOfSamplings, mcSamplingPrefix); // 3 - Perform sampling for n proteins
		}

		if(Boolean.parseBoolean(params.getProperty("calculateNormalDistributionParams"))) {
			ApproximateNormalDistribuiton.getNormalDistributionParams(mcSamplingPrefix, lowerBound, upperBound, numOfSamplings, normalDistributionParamsFile);
		}

		/* Load and test significance annotations */
		if(Boolean.parseBoolean(params.getProperty("testMotifs"))) {
			System.out.println("**Assessing motif clustering**");
			MotifEnrichment m = new MotifEnrichment(distanceMatrix, proteinList2, normalDistributionParamsFile, lowerBound, upperBound,clusteringMeasure, percentThreshold);
			m.testMotifClustering(degenAnnotationPrefix, annotationCompanionFilePrefix, testedDegenMotifsOutputPrefix, Integer.parseInt(args[1]), Integer.parseInt(args[2]));
		}

		/* Look at the overall distribution of significance scores once all annotations have been tested */
		if(Boolean.parseBoolean(params.getProperty("assessSignificanceScores"))) {
			System.out.println("**Assessing significance scores**");
			AssessEnrichment.assessSignificanceScores(testedDegenMotifsOutputPrefix, Integer.parseInt(params.getProperty("numDegenMotifFiles")), significanceScoresFile);
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
}
