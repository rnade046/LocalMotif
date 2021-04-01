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
import utils.Calculator;
import utils.DistanceMatrix;
import utils.GraphLoader;
import utils.Loader;
import utils.NetworkProteins;


public class Main {

	public static void main(String[] args) throws InterruptedException, ExecutionException {

		String wd = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\";

		String cellMapFile = wd + "input_files\\saint-latest-HumanCellMap.txt";
		String baitCellMapMappingFile = wd + "input_files\\baitsMappingFile-HumanCellMap.txt";

		String preyMappingFile = wd + "input_files\\preysMappingFile-BioMart.tsv";
		String baitMappingFile = wd + "input_files\\baitsMappingFile-BioMart.tsv";

		String distanceMatrixFile = wd + "IO_files\\cellMap_distanceMatrix.txt";

		String fastaFile = wd + "input_files\\human_3UTRsequences.txt";
		String rnaIdListFile = wd + "IO_files\\refSeqRNAids-HumanCellMap.tsv";
		ArrayList<Interaction> interactionList =  GraphLoader.loadInteractionRepository(cellMapFile, baitCellMapMappingFile, baitMappingFile, preyMappingFile, 0.01);
		System.out.println("Number of interactions: " + interactionList.size());

		ArrayList<Protein> proteinList = NetworkProteins.getProteinsInNetwork(interactionList);
		System.out.println("Number of Proteins: " + proteinList.size());

		printRefSeqIdsInNetwork(rnaIdListFile, proteinList);
		/*
		printProteinsInNetwork(proteinListFile, proteinList);
		File f = new File(distanceMatrixFile);
		if(!f.exists() && !f.isDirectory()) {
			DistanceMatrix.computeDistanceMatrix(interactionList, proteinList, distanceMatrixFile);
		}
		double[][] distanceMatrix = Loader.loadDistanceMatrix(distanceMatrixFile, proteinList); 

		// TESTING motif enumeration

		ArrayList<Protein> proteinList = new ArrayList<>();
		proteinList.add(new Protein("PLEC", new ArrayList<String>(Arrays.asList("NM_000445"))));
		//SproteinList.add(new Protein("PX6", new ArrayList<String>(Arrays.asList("NM_001316313", "NM_000287"))));

		MotifEnumerator enumerate = new MotifEnumerator(fastaFile, proteinList, 8);

		Timestamp timestamp = new Timestamp(System.currentTimeMillis());
        System.out.println(timestamp + " [start]");
		HashMap<String, HashSet<String>> enumeratedMotifs = enumerate.enumerateMotifs(7);
		timestamp = new Timestamp(System.currentTimeMillis());
        System.out.println(timestamp + " [end]"); */ 
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
}
