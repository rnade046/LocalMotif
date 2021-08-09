import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;

import graph.Interaction;
import graph.Protein;
import opt.CheckDegreeDistributions;
import sampling.ApproximateNormalDistribuiton;
import sampling.MotifSampling;
import sampling.ProteinAnnotations;
import utils.AssessEnrichment;
import utils.Calculator;
import utils.CorrelationGraphLoader;
import utils.DistanceMatrix;
import utils.MotifEnrichment;
import utils.NetworkProteins;

public class Main {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		
		System.out.println("**Loading parameters file** \n");
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));		

		String wd = params.getProperty("working_directory");
		String projectName = params.getProperty("project_name");

		String fastaFile = wd + params.getProperty("fastaFile");
		String mapProtToRefSeqFile = wd + params.getProperty("mapGeneSymbolsToRefSeqIds");

		String correlationRepository = wd + params.getProperty("networkRepositoryFile");
		String proteinsInNetworkOutputFile = wd + projectName +"_listOfProteinNames.tsv";
		String protInfoFile = wd + projectName + "_proteinsInNetwork_info.tsv";
		
		String distanceMatrixFile = wd + projectName + "_distanceMatrix.txt";
		String distanceMatrix2File = wd +  projectName + "_distanceMatrix2.txt";

		String degenAnnotationPrefix = params.getProperty("degenAnnotationPrefix");

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
		ArrayList<Interaction> interactionList = CorrelationGraphLoader.loadGraphFromCorrelationNetwork(correlationRepository, fastaFile, mapProtToRefSeqFile, proteinsInNetworkOutputFile, Double.parseDouble(params.getProperty("corrThreshold")));
		System.out.println("Number of interactions:" + interactionList.size() + "\n");

		System.out.println("**Getting list of proteins in network**");
		ArrayList<Protein> proteinList = NetworkProteins.getProteinsInNetwork(interactionList);
		System.out.println("Number of Proteins: " + proteinList.size() + "\n");
		//printProtAndRefSeqIdsInNetwork(rnaIdListFile, proteinList);
		//printRefSeqIdsInNetwork(rnaIdListFile, proteinList);

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
		CheckDegreeDistributions.assessDegreeDistribution(proteinList2, interactionList, (wd + projectName + "_degreesInNetwork.tsv"));
		
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

		int lowerBound = Integer.parseInt(params.getProperty("lowerBoundToSample", "3"));
		int upperBound = Integer.parseInt(params.getProperty("upperBoundToSample", "2000"));
		int numOfSamplings = Integer.parseInt(params.getProperty("numberOfSamplings"));

		if(Boolean.parseBoolean(params.getProperty("calculateProteinAnnotationFreq"))) {	
			// For MC sampling 
			// 1 - Make list: protein = #motifs (degen + non degen) from full annotation list	>> Do this once
			System.out.println("**Enumerating protein annotation frequency file**");
			ProteinAnnotations freq = new ProteinAnnotations(lowerBound, upperBound);
			freq.calculateProteinAnnotationFrequency(degenAnnotationPrefix, Integer.parseInt(params.getProperty("numDegenMotifFiles")), proteinAnnotationFrequencyFile);
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
			m.testMotifClustering(degenAnnotationPrefix, testedDegenMotifsOutputPrefix, Integer.parseInt(args[1]), Integer.parseInt(args[2]));
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
