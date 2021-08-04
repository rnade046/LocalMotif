package ClusteredMotifs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

public class MotifFamily {


	public static void assessMotifFamilies(String motifFamilyFilePrefix, int numberOfFamilies,
			String clusteringFilePrefix, int numberOfClusteringFiles, String refSeqToMotifFile, 
			String proteinToRefSeqIDsFile, String outputFilePrefix) {

		
		/* Load protein : refseqId list */
		HashMap<String, HashSet<String>> proteinToRefSeqIDsMap = loadProteinsToRefSeqIdsMap(proteinToRefSeqIDsFile);
		
		/* Load refSeqID : motif list */
		HashMap<String, HashSet<String>> refSeqIdToMotifsMap = loadRefSeqIDToMotifsMap(refSeqToMotifFile);
		
		for(int i=1; i<=numberOfFamilies; i++) {

			String motifFamilyFile = motifFamilyFilePrefix + i + ".tsv";

			/* Load motifs in motif family identified during hierarchical clustering step */ 
			HashSet<String> motifSet = loadMotifsInFamily(motifFamilyFile);

			/* Obtain significant clustering information for every motif */
			HashMap<String, Double> motifSignificantMap = getMotifSignificance(motifSet, clusteringFilePrefix, numberOfClusteringFiles);

			/* Identify representative motif from list */
			String representativeMotif = getRepresentativeMotif(motifSignificantMap);

			HashSet<String> possibleMotifSet = getPossibleInstancesOfMotif(representativeMotif);

			ArrayList<String> motifInstances = getInstancesOfMotifs(possibleMotifSet, proteinToRefSeqIDsMap, refSeqIdToMotifsMap);

			String outputFile = outputFilePrefix + i + ".tsv";
			printMotifs(motifInstances, outputFile);
		}
	}

	/**
	 * Load the list of proteins in the network and their contributing refSeq Ids. 
	 * @param proteinInfoFile	String - proteinName \t ID1|ID2|..|IDn
	 * 
	 * @return loadProteinsToRefSeqIdsMap	HashMap<String, HashSet<String>> - map {protein: Set<IDs>}
	 */
	private static HashMap<String, HashSet<String>> loadProteinsToRefSeqIdsMap(String proteinInfoFile){

		HashMap<String, HashSet<String>> proteinToRefSeqIdsMap = new HashMap<>();

		try {
			InputStream in = new FileInputStream(new File(proteinInfoFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				String[] col = line.split("\t");
				proteinToRefSeqIdsMap.put(col[0], new HashSet<String>(Arrays.asList(col[1].split("\\|")))); // col[0] = protein, col[1] = ID1|ID2|ID3

				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinToRefSeqIdsMap;
	}

	
	/**
	 * Load list of RefSeqIDs and their contributing motifs. 
	 * 
	 * @param refSeqIDToMotifsFile
	 * @return
	 */
	private static HashMap<String, HashSet<String>> loadRefSeqIDToMotifsMap(String refSeqIDToMotifsFile){

		HashMap<String, HashSet<String>> refSeqIdToMotifsMap = new HashMap<>();

		try {
			InputStream in = new FileInputStream(new File(refSeqIDToMotifsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				String[] col = line.split("\t");
				refSeqIdToMotifsMap.put(col[0], new HashSet<String>(Arrays.asList(col[1].split("\\|")))); // col[0] = refSeqId, col[1] = motif1|motif2|motif3

				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return refSeqIdToMotifsMap;

	}

	/**
	 * Load motifs in motif family. Store in set. 
	 * @param motifFamilyFile
	 * @return
	 */
	private static HashSet<String> loadMotifsInFamily(String motifFamilyFile){
		HashSet<String> motifSet = new HashSet<>();

		try {
			InputStream in = new FileInputStream(new File(motifFamilyFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {

				motifSet.add(line);
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifSet;
	}

	private static HashMap<String, Double> getMotifSignificance(HashSet<String> motifSet, String clusteringFilePrefix, int numberOfClusteringFiles){

		HashMap<String, Double> motifSignificanceMap = new HashMap<>();

		/* Get p-value for all motifs in this family */ /* rank in ascending order; take smallest one; if tie take motif with least amount of degen characters */ 
		outerloop:
			for(int i=0; i<numberOfClusteringFiles; i++) {

				String clusteringFile = clusteringFilePrefix + i;
				try {
					InputStream in = new FileInputStream(new File(clusteringFilePrefix));
					BufferedReader input = new BufferedReader(new InputStreamReader(in));

					String line = input.readLine();

					while(line != null) {

						String motif = line.split("\t")[0];

						if(motifSet.contains(motif)) {
							motifSignificanceMap.put(motif, Double.parseDouble(line.split("\t")[3])); // [3] = p-value

							if(motifSignificanceMap.size() == motifSet.size()) {
								break outerloop;
							}
						}

						line = input.readLine();
					}

					input.close();
				} catch (IOException e) {
					e.printStackTrace();
				}

			}

		return motifSignificanceMap;
	}

	private static String getRepresentativeMotif(HashMap<String, Double> motifSignificanceMap) {

		String repMotif = "";

		List<Entry<String, Double>> list = new LinkedList<Entry<String, Double>>(motifSignificanceMap.entrySet());  
		list.sort(Entry.comparingByValue());

		double minPval = list.get(0).getValue();
		boolean repeatPval = true; 
		int countIdx = 1;

		/* Check if more than one motif has the same p-value */
		while(repeatPval) {
			if(list.get(countIdx).getValue() == minPval) {
				countIdx++;
			} else {
				repeatPval = false;
			}
		}

		/* Compare motifs with identical p-value*/
		if(countIdx > 1) {

			/* Identify motifs with the same low p-value */
			String[] motifs = new String[countIdx];
			for(int i = 0; i < countIdx; i++) {
				motifs[i] = list.get(i).getKey();
			}

			/* Check for the least degen motif */
			HashSet<Character> degenCharacterSet = new HashSet<>(Arrays.asList('R', 'D', 'H', 'V', 'B', '*'));
			int maxDegenCount = 8;
			for(String motif: motifs) {

				int degenCount = 0;
				for(int k=0; k < motif.length() ; k++) {

					/* keep count if added character is a degenerate character */
					if(degenCharacterSet.contains(motif.charAt(k))) {
						degenCount++;
					}
				}

				if(degenCount < maxDegenCount) {
					repMotif = motif;
					maxDegenCount = degenCount;
				}
			}

		} else { 
			repMotif = list.get(0).getKey();
		}

		return repMotif;
	}

	private static HashSet<String> getPossibleInstancesOfMotif(String repMotif){

		HashMap<Character, Character[]> charMap = defineDegenCharacterMap();
		int solutions = calculateSolutions(repMotif, charMap);

		HashSet<String> possibleMotifSet = getAllMotifs(repMotif, solutions, charMap);


		return possibleMotifSet;
	}

	/**
	 * Define the RNA nucleotide mapping to their possible substitutions (ie. degenerate characters)
	 * Note: for now this is fix - this could be passed as an input parameter in the future to enable 
	 * flexible motif representation for users
	 * Maps RNA nucleotides {A, U, C, G} to possible substitutions (themselves + degenerate characters) 
	 * 
	 * @return charMap	 HashMap<Character, Character[]> - map {nucleotide : list of possible characters}
	 */
	private static HashMap<Character, Character[]> defineDegenCharacterMap(){
		HashMap<Character, Character[]> charMap = new HashMap<>();

		Character[] a = {'A'};
		Character[] t = {'T'};
		Character[] g = {'G'};
		Character[] c = {'C'};
		Character[] r = {'A', 'G'};
		Character[] y = {'T', 'C'};
		Character[] d = {'A', 'T', 'G'};
		Character[] b = {'T', 'G', 'C'};
		Character[] h = {'A', 'T', 'C'};
		Character[] v = {'A', 'G', 'C'};
		Character[] x = {'A', 'T', 'G', 'C'};

		charMap.put('A', a);
		charMap.put('T', t);
		charMap.put('G', g);
		charMap.put('C', c);
		charMap.put('R', r);
		charMap.put('Y', y);
		charMap.put('D', d);
		charMap.put('B', b);
		charMap.put('H', h);
		charMap.put('V', v);
		charMap.put('*', x);

		return charMap;
	}

	/**
	 * Determine the number of possible solutions that will be generated for the given motif length
	 * This function works simply because there's the same amount of possible letters for each nucleotide
	 * Therefore the number of solutions will be the same. 
	 * 
	 * @param motifLength	int - number of nucleotides considered for a motif (motif size)
	 * @param charMap		HashMap<Character, Character[]> - map of {nucleotides : list of substitutions}
	 * 
	 * @return solutions	int - number of solutions for a given motif of length 8
	 */
	private static int calculateSolutions(String motif, HashMap<Character, Character[]> charMap) {
		int solutions = 1;
		for(int i = 0; i < motif.length(); i++) {
			solutions *= charMap.get(motif.charAt(i)).length;
		}
		return solutions;
	} 


	/**
	 * For a given motif produce all possible motifs (ie. degenerate motifs) with given criteria.
	 * Current criteria includes motifs with a max of 7 degenerate characters. 
	 * Inspired by Cartesian Products: https://stackoverflow.com/questions/9591561/java-cartesian-product-of-a-list-of-lists  
	 * 
	 * @param motif				String - sequence motifs 
	 * @return degenerateMotifs ArrayList<String> - all possible degenerate sequence motifs 
	 */
	private static HashSet<String> getAllMotifs(String motif, int solutions, HashMap<Character, Character[]> charMap) {
		HashSet<String> motifs = new HashSet<>();

		/* generate all solutions */ 
		for(int i = 0; i < solutions; i++) {
			int j = 1;
			String seq = "";
			/* generate a given motif */
			for(int k=0; k < motif.length() ; k++) {
				Character[] set = charMap.get(motif.charAt(k));

				char charToAppend = set[(i/j)%set.length]; // magic 
				seq += charToAppend;
				j *= set.length;
			}

			motifs.add(seq);

		}	
		return motifs;
	}

	private static ArrayList<String> getInstancesOfMotifs(HashSet<String> possibleMotifs, HashMap<String, HashSet<String>> proteinToRefSeqIdsMap,
			HashMap<String, HashSet<String>> refSeqIdToMotifsMap) {

		/* Search fasta sequence; find instance of motif (ie. convert degen motif to non degen) ; print PWM to file */
		ArrayList<String> instanceOfMotifList = new ArrayList<>();

		/* Search for motif in a given protein; only one instance of a motif will be kept if found 
		 * multiple times in sequence or in different refseqIds (variants)*/
		for(Entry<String, HashSet<String>> proteinEntry : proteinToRefSeqIdsMap.entrySet()) {
			HashSet<String> currentInstancesOfMotif = new HashSet<>();
			
			for(String id: proteinEntry.getValue()) {
				for(String motif: refSeqIdToMotifsMap.get(id)) {
					if(possibleMotifs.contains(motif)) {
						currentInstancesOfMotif.add(motif);
					}
				}
				
			}
			instanceOfMotifList.addAll(currentInstancesOfMotif);
		}

		return instanceOfMotifList;
	}

	private static void printMotifs(ArrayList<String> motifInstances, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(String motif: motifInstances) {
				out.write(motif + "\n");
				out.flush();
			}


			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
