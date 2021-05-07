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
import sampling.ApproximateNormalDistribuiton;
import sampling.MotifSampling;
import sampling.ProteinAnnotations;
import utils.Calculator;
import utils.CorrelationGraphLoader;
import utils.DistanceMatrix;
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

		String distanceMatrixFile = wd + projectName + "_distanceMatrix.txt";
		String distanceMatrix2File = wd +  projectName + "_distanceMatrix2.txt";
		
		String annotationFile = wd + params.getProperty("motifAnnotationFile");
		String degenAnnotationPrefix = wd + params.getProperty("degenAnnotationPrefix");

		String proteinAnnotationFrequencyFile = wd + projectName + "_protFreqAnnotation.tsv";
		String mcSamplingPrefix = wd + projectName + "_mcSamplingDistribution_";
		String normalDistributionParamsFile = wd + projectName + "_normalDistributionParams.tsv";

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
		ArrayList<Protein> proteinList2 = NetworkProteins.modifyNetworkProteinsList(proteinList, proteinsToKeep);

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

		if(Boolean.parseBoolean(params.getProperty("calculateProteinAnnotationFreq"))) {	
			// For MC sampling 
			// 1 - Make list: protein = #motifs (degen + non degen) from full annotation list	>> Do this once
			System.out.println("**Enumerating protein annotation frequency file**");
			ProteinAnnotations freq = new ProteinAnnotations(Integer.parseInt(params.getProperty("lowerBoundToSample", "3")), Integer.parseInt(params.getProperty("upperBoundToSample", "2000")));
			freq.calculateProteinAnnotationFrequency(annotationFile, degenAnnotationPrefix, Integer.parseInt(params.getProperty("numDegenMotifFiles")), proteinAnnotationFrequencyFile);
		}
		
		
		int lowerBound = Integer.parseInt(args[1]);
		int upperBound = Integer.parseInt(args[2]);
		
		int numOfSamplings = Integer.parseInt(params.getProperty("numberOfSamplings"));
		
		if(Boolean.parseBoolean(params.getProperty("performMCprocedure"))) {
			
			
			/* This will be done in job arrays on CC */
			System.out.println("**Performing Monte Carlo Sampling Procedure**");
			// 2 - Initialize sampling
			MotifSampling sampling = new MotifSampling(proteinAnnotationFrequencyFile, proteinList2, distanceMatrix);
			// 3 - Perform sampling for n proteins
			sampling.computeMultipleDistributions(lowerBound, upperBound, numOfSamplings, mcSamplingPrefix);
		}

		if(Boolean.parseBoolean(params.getProperty("calculateNormalDistributionParams"))) {
			ApproximateNormalDistribuiton.getNormalDistributionParams(mcSamplingPrefix, lowerBound, upperBound, numOfSamplings, normalDistributionParamsFile);
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
