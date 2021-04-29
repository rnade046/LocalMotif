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
import sampling.MotifSampling;
import sampling.ProteinAnnotations;
import utils.CorrelationGraphLoader;
import utils.DistanceMatrix;
import utils.NetworkProteins;


public class Main {

	public static void main(String[] args) throws FileNotFoundException, IOException {
		
		System.out.println("**Loading parameters file** \n");
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));		
		
		String wd = params.getProperty("working_directory");

		String fastaFile = wd + params.getProperty("fastaFile");
		String mapProtToRefSeqFile = wd + params.getProperty("mapGeneSymbolsToRefSeqIds");

		String correlationRepository = wd + params.getProperty("networkRepositoryFile");
		String proteinsInNetworkOutputFile = wd + "CorrelationNet_listOfProteinNames.tsv";

		String distanceMatrixFile = wd + params.getProperty("distanceMatrixFile");

		String annotationFile = wd + params.getProperty("motifAnnotationFile");
		String degenAnnotationPrefix = wd + params.getProperty("degenAnnotationPrefix");
		
		String proteinAnnotationFrequencyFile = wd + params.getProperty("proteinAnnotationFrequencyFile");
		String mcSamplingPrefix = wd + params.getProperty("mcSamplingPrefix");

		System.out.println("**Loading interaction repository**");
		ArrayList<Interaction> interactionList = CorrelationGraphLoader.loadGraphFromCorrelationNetwork(correlationRepository, fastaFile, mapProtToRefSeqFile, proteinsInNetworkOutputFile, 0.5);
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


		if(Boolean.parseBoolean(params.getProperty("calculateProteinAnnotationFreq"))) {	
			// For MC sampling 
			// 1 - Make list: protein = #motifs (degen + non degen) from full annotation list	>> Do this once
			File f1 = new File(proteinAnnotationFrequencyFile);
			if(!f1.exists() && !f1.isDirectory()) {
				System.out.println("Enumerating protein annotation frequency file");
				ProteinAnnotations.enumerateProteinAnnotationFrequency(annotationFile, proteinAnnotationFrequencyFile);
			}
		}
		if(Boolean.parseBoolean(params.getProperty("performMCprocedure"))) {
			/* This will be done in job arrays on CC */
			System.out.println("**Performing Monte Carlo Sampling Procedure**");
			// 2 - Initialize sampling
			MotifSampling sampling = new MotifSampling(proteinAnnotationFrequencyFile, proteinList, distanceMatrix);
			// 3 - Perform sampling for n proteins
			sampling.computeMultipleDistributions(3, 3, 10000, mcSamplingPrefix);
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
