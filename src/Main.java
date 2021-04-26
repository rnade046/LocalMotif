import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.Timestamp;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import graph.Interaction;
import graph.Protein;
import motifs.MotifDegeneration;
import sampling.MotifSampling;
import sampling.ProteinAnnotations;
import utils.Calculator;
import utils.CorrelationGraphLoader;
import utils.DistanceMatrix;
import utils.SaintGraphLoader;
import utils.Loader;
import utils.NetworkProteins;


public class Main {

	public static void main(String[] args) {

		String wd = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\";

		String cellMapFile = wd + "input_files\\saint-latest-HumanCellMap.txt";

		String mapProtToRefSeqFile = wd + "motif_enumeration\\BiomaRt_MappingRefSeqIdsToGeneSymbol_corrNet.tsv";
		
		String correlationRepository = wd + "input_files\\correlation_v2.tsv";
		String proteinsInNetwork = wd + "CorrelationNet_listOfProteinNames.tsv";
		String distanceMatrixFile = wd + "IO_files\\distanceMatrix_cellMap_CorrNet.txt";
		
		String fastaFile = wd + "input_files\\human_3UTRsequences.txt";
		String annotationFile = wd + "IO_files\\motifMapped.tsv";
		String proteinAnnotationFrequencyFile = wd + "IO_files\\test-protFreqAnnotation.tsv";
		String mcSamplingFile = wd + "IO_files\\test-MonteCarloSamplingFile_n3_s100000.tsv";
		
		System.out.println("**Loading interaction repository**");
		ArrayList<Interaction> interactionList = CorrelationGraphLoader.loadGraphFromCorrelationNetwork(correlationRepository, fastaFile, mapProtToRefSeqFile, proteinsInNetwork, 0.5);
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
		/*
		double[][] distanceMatrix = DistanceMatrix.loadDistanceMatrix(distanceMatrixFile, proteinList); 
		
		// For MC sampling 
		// 1 - Make list: protein = #motifs (degen + non degen) from full annotation list
		File f1 = new File(proteinAnnotationFrequencyFile);
		if(!f1.exists() && !f1.isDirectory()) {
			System.out.println("Enumerating protein annotation frequency file");
			ProteinAnnotations.enumerateProteinAnnotationFrequency(annotationFile, proteinAnnotationFrequencyFile);
		}
		// 2 - Initialize sampling
		MotifSampling sampling = new MotifSampling(proteinAnnotationFrequencyFile, proteinList, distanceMatrix);
		// 3 - Perform sampling for n proteins
		sampling.computeMultipleDistributions(3, 3, 10000, mcSamplingFile);*/
		
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
