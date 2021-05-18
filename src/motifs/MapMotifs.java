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
import java.util.Map;

public class MapMotifs {

	public static void mapDegenMotifsToRefSeqIds(String listDegenMotifsToTestInputFile, String motifToDegenMotifsFile, 
			int numMotifFiles, String motifsToRefSeqIDFile, String degenMotifsToRefSeqIdsFile, String protToRefSeqFile) {


		/* Load degenerative motifs that need to be indexed */
		List<String> degenMotifList = loadDegenerativeMotifsToMap(listDegenMotifsToTestInputFile);
		System.out.println("Loaded degen motifs to test: " + degenMotifList.size() + "\n");

		System.out.println("Iterating through degen motif map: ");
		HashMap<String, List<String>> degenMotifToRefSeqIdMap = mapDegenerateMotifToRefSeqIds(motifToDegenMotifsFile, numMotifFiles, degenMotifList);
		
		/* Load (65,466) motifs and their associated RefSeqIds */
		HashMap<String, List<String>> motifsToRefSeqMap = loadMotifsToRefSeqMap(motifsToRefSeqIDFile);
		System.out.println("Loaded motifs to refseq map: " + motifsToRefSeqMap.size());		
		HashMap<String, String> refSeqToProtMap = loadRefSeqToProteinMap(protToRefSeqFile);
		System.out.println("Loaded refSeq to protein map: " + refSeqToProtMap.size() + "\n");
		
		System.out.println("Print map of Degen Motifs To RefSeq Ids:");
		printMapOfDegenMotifsToRefSeqIds(degenMotifsToRefSeqIdsFile, degenMotifToRefSeqIdMap, motifsToRefSeqMap, refSeqToProtMap);

		
//		//for(String motif: degenMotifList) {
//			
//			/* Map degenMotifs to all associated motifs */
//			
//			//System.out.println("Testing motif: " + motif);
//			
//			HashSet<String> refSeqIdSet = mapDegenerateMotifToRefSeqIdsForAGivenMotif(motifToDegenMotifsFile, motif, motifsToRefSeqMap);
//			/* print degenMotif annotation */
//			
//			printDegenMotifsToRefSeqIds(degenMotifsToRefSeqIdsFile, motif, refSeqIdSet);
//		//}
		
	}
	
	public static void mapMotifsToRefSeqIds(String motifsToRefSeqIDFile, String protToRefSeqFile, String motifMapFile) {
	
		HashMap<String, List<String>> motifsToRefSeqMap = loadMotifsToRefSeqMap(motifsToRefSeqIDFile);
		System.out.println("Loaded motifs to refseq map: " + motifsToRefSeqMap.size());		
		HashMap<String, String> refSeqToProtMap = loadRefSeqToProteinMap(protToRefSeqFile);
		System.out.println("Loaded refSeq to protein map: " + refSeqToProtMap.size() + "\n");
		
		printMapOfMotifsToRefSeqIds(motifMapFile, motifsToRefSeqMap, refSeqToProtMap);
		
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
	
	private static HashMap<String, String> loadRefSeqToProteinMap(String inputFile){
		HashMap<String, String> refSeqToProteinMap = new HashMap<>();
		
		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String line = input.readLine(); // header 
			line = input.readLine();
			
			while(line != null) {
				
				String prot = line.split("\t")[0];
				String refSeqId = line.split("\t")[1];
				
				refSeqToProteinMap.put(refSeqId, prot);
				
				line = input.readLine();
			}
			
			input.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return refSeqToProteinMap;
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
	 * @param degenMotifsToTestList 			List<String> - list of degenerate motifs to map
	 * 
	 * @return degenMotifsToNonDegenMotifMap	HashMap<String, HashSet<String>> - map of degenMotif = RefSeqID1|RefSeqID2|...|RefSeqIDX
	 */
	private static HashMap<String, List<String>> mapDegenerateMotifToRefSeqIds(String mapMotifToDegenMotifsFile, int numMotifFile, List<String> degenMotifsToTestList){
		HashMap<String, List<String>> degenMotifsToMotifsMap = new HashMap<>();

		
		for(int i=0; i<=numMotifFile; i++) {
			// mapMotifToDegenMotifsFile = mapMotifsToDegenMotifs_
			String file = mapMotifToDegenMotifsFile + i;
			
			System.out.print("File " + i + ": ");
			
			InputStream in;
			try {
				in = new FileInputStream(new File(file));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				/* Ex. line = motifX \t DegenMotif1|DegenMotif2|...|DegenMotifX */
				String line = input.readLine(); // no header
				int lineCount = 1;
				while(line != null) {

					System.out.print(lineCount +".");
					
					String motif = line.split("\t")[0];	// motif
					HashSet<String> currentDegenMotifsSet = new HashSet<String>(Arrays.asList(line.split("\t")[1].split("\\|"))); // List of degenerate motifs associated to "motif"

					/* for each degen motif that we're trying to map, check if it is associated to the current "motif", 
					 * if so add the current "motif's" RefSeqIds to the degenMotifs list of associated motifs*/
					for(String degenMotif: degenMotifsToTestList) {
						if(currentDegenMotifsSet.contains(degenMotif)) {
							if(degenMotifsToMotifsMap.containsKey(degenMotif)) {
								degenMotifsToMotifsMap.get(degenMotif).add(motif);
							} else {
								List<String> associatedMotifList = new ArrayList<>();
								associatedMotifList.add(motif);
								degenMotifsToMotifsMap.put(degenMotif, associatedMotifList);
							}
						}
					}
					
					line = input.readLine();
					lineCount++;
				}
				System.out.println("current map size: " + degenMotifsToMotifsMap.size());
				input.close();
			} catch(IOException e) {
				e.printStackTrace();
				e.getCause();
				e.getMessage();
			}
			
		}
		
		System.out.print("Done\n\n");
		return degenMotifsToMotifsMap;
	}

	@SuppressWarnings("unused")
	private static HashSet<String> mapDegenerateMotifToRefSeqIdsForAGivenMotif(String mapMotifToDegenMotifsFile, String degenMotif, HashMap<String, List<String>>  motifsToRefSeqMap){
		HashSet<String> refSeqIDSet = new HashSet<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(mapMotifToDegenMotifsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			/* Ex. line = motifX \t DegenMotif1|DegenMotif2|...|DegenMotifX */
			String line = input.readLine(); // no header
			int lineCount = 1;
			System.out.print(lineCount + ".");
			while(line != null) {

				System.out.print(lineCount + ".");
//				if(lineCount%100000 == 0) {
//					System.out.print(lineCount + ".");
//				} 
//
//				if(lineCount%1000000 == 0) {
//					System.out.println();
//				}
				String motif = line.split("\t")[0];	// motif
				HashSet<String> degenMotifsSet = new HashSet<String>(Arrays.asList(line.split("\t")[1].split("\\|"))); // List of degenerate motifs associated to "motif"

				/* for each degen motif that we're trying to map, check if it is associated to the current "motif", 
				 * if so add the current "motif's" RefSeqIds to the degenMotifs list of associated motifs*/
				if(degenMotifsSet.contains(degenMotif)) {
					if(motifsToRefSeqMap.containsKey(motif)){
						refSeqIDSet.addAll(motifsToRefSeqMap.get(motif));
					}
				}

				line = input.readLine();
				lineCount++;

			}
			input.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		System.out.print("Done\n");
		return refSeqIDSet;
	}

	private static void printMapOfDegenMotifsToRefSeqIds(String outputFile, HashMap<String, List<String>> degenMotifsToAssociatedMotifsMap, 
		HashMap<String, List<String>> motifsToRefSeqIds, HashMap<String, String> refSeqToProtMap) {
		int maxNumberOfProteins = 0;
		int minNumberOfProteins = Integer.MAX_VALUE;
		int motifCount = 0;
		
		System.out.println("Number of degen motifs to iterate: " + degenMotifsToAssociatedMotifsMap.size());
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(Map.Entry<String, List<String>> entry : degenMotifsToAssociatedMotifsMap.entrySet()) {
				
				String degenMotif = entry.getKey();
				List<String> associatedMotifs = entry.getValue();
				
				motifCount++;
				if(motifCount%100 == 0) {
					System.out.print(motifCount + ".");
				}
				
				if(motifCount%1000 == 0) {
					System.out.println();
				}
				
				System.out.println(degenMotif + "|associated to non degenMotifs: " + associatedMotifs.size());
				if(degenMotifsToAssociatedMotifsMap.get(degenMotif).size() >= 1) {
					
					out.write(degenMotif + "\t");
					
					HashSet<String> refSeqSet = getRefSeqIdSet(associatedMotifs, motifsToRefSeqIds);
					HashSet<String> proteinSet = getProteinAssociatedToDegenMotif(refSeqSet, refSeqToProtMap);
					if(proteinSet.size() > maxNumberOfProteins) {
						maxNumberOfProteins = proteinSet.size();
					}
					if(proteinSet.size() < minNumberOfProteins) {
						minNumberOfProteins = proteinSet.size();
					}
					
					out.write(refSeqSet.size() + "\t" + proteinSet.size()  +"\t");
					
					for(String refSeqId: refSeqSet) {
						out.write(refSeqId + "|");
					}
					out.write("\t");
					for(String prot: proteinSet) {
						out.write(prot + "|");
					}
				}
				out.write("\n");
			}
			System.out.println("Done\n");
			
			System.out.println("Max #proteins associated to motifs: " + maxNumberOfProteins);
			System.out.println("Min #proteins associated to motifs: " + minNumberOfProteins);
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	private static HashSet<String> getRefSeqIdSet(List<String> listOfMotifs, HashMap<String, List<String>> motifsToRefSeqMap){
		HashSet<String> refSeqIdSet = new HashSet<>();
		
		for(String motif: listOfMotifs) {
			if(motifsToRefSeqMap.containsKey(motif)) {
				refSeqIdSet.addAll(motifsToRefSeqMap.get(motif));	
			}
		}		
		return refSeqIdSet;
	}
	
	private static HashSet<String> getProteinAssociatedToDegenMotif(HashSet<String> refSeqIds, HashMap<String, String> refSeqToProtMap){
		HashSet<String> proteinSet = new HashSet<>();
	
		for(String refSeqId : refSeqIds) {
			if(refSeqToProtMap.containsKey(refSeqId)) {
				proteinSet.add(refSeqToProtMap.get(refSeqId));	
			}
		}			
		
		return proteinSet;
	}
	
	private static void printMapOfMotifsToRefSeqIds(String outputFile, HashMap<String, List<String>> motifsToRefSeqIds, HashMap<String, String> refSeqToProtMap) {
			int maxNumberOfProteins = 0;
			int minNumberOfProteins = Integer.MAX_VALUE;
			int motifCount = 0;
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
				
				for(String motif:  motifsToRefSeqIds.keySet()) {
					motifCount++;
					if(motifCount%100 == 0) {
						System.out.print(motifCount + ".");
					}
					
					if(motifCount%1000 == 0) {
						System.out.println();
					}
											
						out.write(motif + "\t");
						
						HashSet<String> refSeqSet = new HashSet<>(motifsToRefSeqIds.get(motif));
						HashSet<String> proteinSet = getProteinAssociatedToDegenMotif(refSeqSet, refSeqToProtMap);
						
						if(proteinSet.size() > maxNumberOfProteins) {
							maxNumberOfProteins = proteinSet.size();
						}
						
						if(proteinSet.size() < minNumberOfProteins) {
							minNumberOfProteins = proteinSet.size();
						}
						
						out.write(refSeqSet.size() + "\t" + proteinSet.size() + "\t");
						
						for(String refSeqId: refSeqSet) {
							out.write(refSeqId + "|");
						}
						out.write("\t");
						
						for(String prot: proteinSet) {
							out.write(prot + "|");
						}
					out.write("\n");
				}
				System.out.println("Done\n");
				
				System.out.println("Max #proteins associated to motifs: " + maxNumberOfProteins);
				System.out.println("Min #proteins associated to motifs: " + minNumberOfProteins);
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

		}

	@SuppressWarnings("unused")
	private static void printDegenMotifsToRefSeqIds(String outputFile, String degenMotif, HashSet<String> refSeqIDSet) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile), true));


					out.write(degenMotif + "\t");

					for(String refSeqId: refSeqIDSet) {
						out.write(refSeqId + "|");
					}
				
				out.write("\n");
			


			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
