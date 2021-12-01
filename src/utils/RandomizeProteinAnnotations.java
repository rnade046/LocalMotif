package utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import graph.Protein;

public class RandomizeProteinAnnotations {

	public static void generateNullModel(ArrayList<Protein> proteinList, String nullMapFile) {
		
		/* Get list of proteins with RefSeqIds */
		HashMap<String, List<String>> annotatedProteinsMap = getAnnotatedProteins(proteinList);
		ArrayList<String> annotatedProteinList = new ArrayList<>(annotatedProteinsMap.keySet());
		
		/* Randomize proteins and their 3'UTRs */
		int[] proteinOrder = randomizeProteinsOrder(annotatedProteinsMap.size());
		
		/* Generate null map of proteins to refSeq IDs*/
		generateMapProteinToRefSeqIDs(annotatedProteinsMap, annotatedProteinList, proteinOrder, nullMapFile);
	}
	
	private static HashMap<String, List<String>> getAnnotatedProteins(ArrayList<Protein> proteinList){
		
		HashMap<String,	List<String>> proteinMap = new HashMap<>();
		
		for(Protein protein: proteinList) {
			if(protein.getProteinId() != null) {
				proteinMap.put(protein.getProteinName(), protein.getProteinId());
			}
		}
		return proteinMap;
	}
	
	private static int[] randomizeProteinsOrder(int numProteins) {
	
		/* Initialize array */
		int[] proteins = new int[numProteins];
		for(int i=0; i<numProteins; i++) {
			proteins[i]=i;
		}
		
		/* Randomize array indexes */
		Random ran = new Random();
		for(int i=0; i<numProteins*1000; i++) {
			
			// Select two indexes
			int idx1 = ran.nextInt(numProteins);
			int idx2 = ran.nextInt(numProteins);
			
			// swap the indexes
			proteins[idx1] = idx2;
			proteins[idx2] = idx1; 
		}
		
		return proteins;
	}
	
	private static void generateMapProteinToRefSeqIDs(HashMap<String, List<String>> annotatedProteinMap, ArrayList<String> annotatedProteinsList, int[] proteinOrder, String nullMapFile) {
		
		try {
			
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(nullMapFile)));
			
			for(int i=0; i<proteinOrder.length; i++) {
				
				String protein1 = annotatedProteinsList.get(i); // protein name 
				String protein2 = annotatedProteinsList.get(proteinOrder[i]); // protein whose annotations to switch with
				
				List<String> refSeqIds = annotatedProteinMap.get(protein2); // new RefSeqIDs
 				
				out.write(protein1 + "\t");
				for(String id: refSeqIds) {
					out.write(id + "|");
				}
				out.write("\n");
			}
			
			out.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		
	}
}
