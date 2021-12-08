package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
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
		
		int maxShuffle = ((numProteins*(numProteins-1))/2)*100;
		for(int i=0; i<maxShuffle; i++) {
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

	/**
	 * Generate map of protein = motif1|2|3 for annotation files */ 
	public static void generateProteinToMotifMap(String nullMapFile, String refSeqToMotifsFile, String outputProteinFile) {

		/* Load proteins to randomized refSeqIds */
		HashMap<String, String> proteinToRefSeqIdsMap = loadProtToRefSeq(nullMapFile); // key = protein, values = refSeqIds
		System.out.println("Loaded protein to refSeqIds: " + proteinToRefSeqIdsMap.size());
		/* Load refSeqIds to motifs */
		HashMap<String, List<String>> motifToRefSeqIDsMap = loadMotifsToRefSeqMap(refSeqToMotifsFile); // key = motif, values = refSeqIds
		System.out.println("Loaded refSeqID to motifs: " + motifToRefSeqIDsMap.size());

		/* Map proteins to motifs */
		System.out.println("Mapping proteins to motifs");
		mapProteinsToMotifs(proteinToRefSeqIdsMap, motifToRefSeqIDsMap, outputProteinFile);
	}

	
	
	private static HashMap<String, String> loadProtToRefSeq(String mapFile){

		HashMap<String, String> map = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(mapFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				String prot = line.split("\t")[0];
				
				if(line.split("\t").length > 1 ) {
				String[] values = line.split("\t")[1].split("\\|");
				for(String id : values) {
					map.put(id, prot);
				}
				
				}
				line = input.readLine();
			}

			input.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return map;
	}
	
	/**
	 * Load motifs and their associated RefSeq Ids contained in text file into HashMap
	 * Map should contain 65,466 motifs. 
	 * 
	 * @param inputFile				String - file path; file format: RefSeqId \t Motif1|Motif2|..|Motifx
	 * @return motifsToRefSeqMap	HashMap<String, List<String> - MotifX = RefSeqId1, RefSeqId2, ..., RefSeqIdX
	 */
	private static HashMap<String, List<String>> loadMotifsToRefSeqMap(String inputFile){

		HashMap<String, List<String>> motifsToRefSeqMap = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			/* Ex. line: RefSeqID \t Motif1|Motif2|Motif3|..|MotifX */
			String line = input.readLine(); // no header

			while(line!=null) {

				if(line.split("\t").length >= 2) {
					String refSeqId = line.split("\t")[0];				//RefSeqId
					String[] motifs = line.split("\t")[1].split("\\|");	//Motifs associated to refSeqId

					/* map motifX = list of RefSeqIds*/
					for(String m: motifs) {
						if(motifsToRefSeqMap.containsKey(m)) {
							motifsToRefSeqMap.get(m).add(refSeqId);		// add current RefSeq Id to existing list
						} else {
							List<String> refSeqIdList = new ArrayList<>();	// create new list
							refSeqIdList.add(refSeqId);
							motifsToRefSeqMap.put(m, refSeqIdList); 	// add entry for new motif
						}
					}
				} else { 
					//System.out.println("error: " + line.split("\t")[0]);
				}
				line = input.readLine();
			}
			input.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return motifsToRefSeqMap;
	}

	private static void mapProteinsToMotifs(HashMap<String, String> proteinToRefSeqIdsMap, HashMap<String, List<String>> refSeqToMotifsMap, String outputFile) {

		try {

			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(Entry<String, List<String>> motifEntry : refSeqToMotifsMap.entrySet()) { // motif = refSeqID1|2|n
				List<String> refSeqIds = motifEntry.getValue();

				HashSet<String> proteinSet = new HashSet<>();

				for(String id: refSeqIds) {
					if(proteinToRefSeqIdsMap.containsKey(id)) { // contains refseq id
						proteinSet.add(proteinToRefSeqIdsMap.get(id)); // get protein name
					}
				}
			
				if(!proteinSet.isEmpty()) {
				out.write(motifEntry.getKey() + "\t"); // motif
				
				for(String protein: proteinSet) {
					out.write(protein + "|");
				}
				out.write("\n");
				
				out.flush();
			
				}
			}

		out.close();
	} catch(IOException e) {
		e.printStackTrace();
	}

}
}
