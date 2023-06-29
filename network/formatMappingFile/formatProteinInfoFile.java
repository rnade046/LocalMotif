package formatMappingFile;

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
import java.util.Map.Entry;

public class formatProteinInfoFile {

	public static void main(String[] args) {

		String fastaFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/input_files/human_3UTRsequences.txt";
		String biomartFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motif_enumeration/BiomaRt_MappingRefSeqIdsToGeneSymbol_corrNet.tsv";
		String correlationNetworkFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/hcm/corrNet2-400_formattedNetwork.tsv";

		String outputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/hcm/corrNet2-400_protein-refSeq-info.tsv";

		/* Load map Protein = RefSeqIds from BioMart mapping - filter by RefSeqIds for which we have a sequence */
		HashMap<String, ArrayList<String>> mapRefSeqIds = getRefSeqIdsForwhichWeHaveASequence(biomartFile, fastaFile);

		/* Update mapping file to reflect only proteins in original correlation network */
		mapRefSeqIds = updateMappingForProteinsInNetwork(correlationNetworkFile, mapRefSeqIds);

		/* output formatted mapping file (network-info) */
		printMappingFile(outputFile, mapRefSeqIds);
	}

	/**
	 * Get the list of refSeqIds associated to a given protein. RefSeqIds need to have a corresponding sequence in the fasta file
	 * 
	 * @param biomartMappingFile		String - file path output from BiomaRt outlining map of HGNC symbol to RefSeqId
	 * @param refSeqSet					String - Fasta file with sequences
	 * @return mapRefSeqIds				HashMap<String, ArrayList<String>> - map of protein as HGNC symbol = list of RefSeqIds
	 */
	private static HashMap<String, ArrayList<String>> getRefSeqIdsForwhichWeHaveASequence(String biomartMappingFile, String fastaFile){

		HashMap<String, ArrayList<String>> mapRefSeqIds = new HashMap<>();

		HashSet<String> refSeqSet = generateRefSeqSet(fastaFile);

		try {
			InputStream in = new FileInputStream(new File(biomartMappingFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();
			while(line!=null) {

				String protein = line.split("\t")[0];
				String refSeqID = line.split("\t")[1];

				/* Add RefSeqId only if it has corresponding sequence in fasta file */
				if(refSeqSet.contains(refSeqID)) {

					/* Add refSeqId to existing list if protein has been seen b4, 
					 * or create new listing if protein hasn't been seen b4 */ 
					if(mapRefSeqIds.containsKey(protein)) {
						mapRefSeqIds.get(protein).add(refSeqID);
					} else {
						ArrayList<String> refSeqIdList = new ArrayList<>();
						refSeqIdList.add(refSeqID);
						mapRefSeqIds.put(protein, refSeqIdList);
					}
				}
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return mapRefSeqIds;
	}

	/**
	 * Generate an index for the FASTA file; identifying the line number for the different refseq identifiers
	 * 
	 * @param fastaFile String - file path for the FASTA sequences
	 * @return indexOfFastaFile HashMap<String, Integer> - map of {refseqId : line count} 
	 */
	private static HashSet<String> generateRefSeqSet(String fastaFile){
		HashSet<String> refSeqSet = new HashSet<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(fastaFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				/* store the line index of the start of a new sequence */
				if(line.startsWith(">")) {
					String[] col = line.split("_|\\.");
					String refSeqId = col[2] + "_" + col[3]; // col[2] = type of ID (eg. NM) ; col[3] = number ID (#####)

					refSeqSet.add(refSeqId);
				}
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return refSeqSet;
	}

	private static HashMap<String, ArrayList<String>> updateMappingForProteinsInNetwork(String networkFile, HashMap<String, ArrayList<String>> mapProtToRefSeqIds){

		HashMap<String, ArrayList<String>> updatedMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(networkFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {
				String[] col = line.split("\t");

				/* check interactor 1 */
				if(!updatedMap.containsKey(col[0])){
					if(mapProtToRefSeqIds.containsKey(col[0])) {
						updatedMap.put(col[0], mapProtToRefSeqIds.get(col[0]));
					} else {
						updatedMap.put(col[0], new ArrayList<String>());
					}

				}

				/* check interactor 2 */
				if(!updatedMap.containsKey(col[1])){
					if (mapProtToRefSeqIds.containsKey(col[1])) {
						updatedMap.put(col[1], mapProtToRefSeqIds.get(col[1]));
					} else {
						updatedMap.put(col[1], new ArrayList<String>());
					}
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return updatedMap;
	}

	private static void printMappingFile(String outputFile, HashMap<String, ArrayList<String>> mapProtToRefSeqIds) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(Entry<String, ArrayList<String>> mapping : mapProtToRefSeqIds.entrySet()) {

				out.write(mapping.getKey() + "\t"); // protein name

				if(!mapping.getValue().isEmpty()) {
					for(String id: mapping.getValue()) { // list refSeqIds 
						out.write(id + "|");
					}
				}
				out.write("\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
