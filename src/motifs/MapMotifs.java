package motifs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class MapMotifs {

	public static void mapMotifsToRefSeqIds(String motifsToRefSeqIDFile, String degenMotifsFile, String motifToDegenMotifsFile, String degenMotifsToRefSeqIdsFile) {

		/* Load (65,466) motifs and their associated RefSeqIds */
		HashMap<String, List<String>> motifsToRefSeqMap = loadMotifsToRefSeqMap(motifsToRefSeqIDFile);

		/* Load degenerative motifs that need to be indexed */
		List<String> degenMotifList = loadDegenerativeMotifsToMap(degenMotifsFile);

		/* Map degenMotifs to all associated motifs */
		HashMap<String, HashSet<String>> degenMotifToRefSeqIdMap = mapDegenerateMotifToRefSeqIds(motifToDegenMotifsFile, degenMotifList,  motifsToRefSeqMap);
		/* print degenMotif annotation */
		printMapOfDegenMotifsToRefSeqIds(degenMotifsToRefSeqIdsFile, degenMotifToRefSeqIdMap, motifsToRefSeqMap);
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

	/**
	 * Load degenerate motifs to index as a list<String> 
	 * 
	 * @param inputFile	String - file path that contains list of degenerate motifs
	 * @return degenMotifList List<String> - list of degenerate motifs
	 */
	private static List<String> loadDegenerativeMotifsToMap(String inputFile){
		List<String> degenMotifList = new ArrayList<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			while(line != null) {

				degenMotifList.add(line);
				line = input.readLine();
			}

			input.close();
		} catch(IOException e) {
			e.printStackTrace();
		}

		return degenMotifList;
	}
	
	/**
	 * Map degenerate motifs to it's corresponding non degenerate motifs. 
	 * Takes the list of degenerate motif that we are looking to map, runs through the file containing motif = list of degenMotifs, &
	 * identifies all non degenerate motif the degenMotif is associated to
	 * 
	 * @param mapMotifToDegenMotifsFile	String - file path, format: motif \t degenMotif1|degenMotif2|...|degenMotifX
	 * @param degenMotifList 			List<String> - list of degenerate motifs to map
	 * 
	 * @return degenMotifsToNonDegenMotifMap	HashMap<String, HashSet<String>> - map of degenMotif = RefSeqID1|RefSeqID2|...|RefSeqIDX
	 */
	private static HashMap<String, HashSet<String>> mapDegenerateMotifToRefSeqIds(String mapMotifToDegenMotifsFile, List<String> degenMotifList, HashMap<String, List<String>>  motifsToRefSeqMap){
		HashMap<String, HashSet<String>> degenMotifsToRefSeqIDsMap = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(mapMotifToDegenMotifsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			/* Ex. line = motifX \t DegenMotif1|DegenMotif2|...|DegenMotifX */
			String line = input.readLine(); // no header

			while(line != null) {

				String motif = line.split("\t")[0];	// motif
				HashSet<String> degenMotifsSet = new HashSet<String>(Arrays.asList(line.split("\t")[1].split("\\|"))); // List of degenerate motifs associated to "motif"

				/* for each degen motif that we're trying to map, check if it is associated to the current "motif", 
				 * if so add the current "motif's" RefSeqIds to the degenMotifs list of associated motifs*/
				for(String degenMotif: degenMotifList) {
					if(degenMotifsSet.contains(degenMotif)) {
						if(degenMotifsToRefSeqIDsMap.containsKey(degenMotif)) {
							degenMotifsToRefSeqIDsMap.get(degenMotif).addAll(motifsToRefSeqMap.get(motif));
						} else {
							HashSet<String> associatedRefSeqIDList = new HashSet<>();
							associatedRefSeqIDList.addAll(motifsToRefSeqMap.get(motif));
							degenMotifsToRefSeqIDsMap.put(degenMotif, associatedRefSeqIDList);
						}
					}
				}
				line = input.readLine();
			}

			input.close();
		} catch(IOException e) {
			e.printStackTrace();
		}

		return degenMotifsToRefSeqIDsMap;
	}

	
	private static void printMapOfDegenMotifsToRefSeqIds(String outputFile, HashMap<String, HashSet<String>> degenMotifsToAssociatedRefSeqMap, HashMap<String, List<String>> motifsToRefSeqIds) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(String degenMotif:  degenMotifsToAssociatedRefSeqMap.keySet()) {

				if(degenMotifsToAssociatedRefSeqMap.get(degenMotif).size() >= 1) {

					out.write(degenMotif + "\t");

						for(String refSeqId:  degenMotifsToAssociatedRefSeqMap.get(degenMotif)) {
								out.write(refSeqId + "|");
						}
					}
					out.write("\n");
				}
			

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
