

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
import java.util.List;

public class GetSequencesInNetworkForMEME {

	public static void main(String[] args) {

		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";
		String proteinInfoFile = wd + "corrNetTop2-400_proteinsInNetwork_info.tsv";
		String fastaFile = wd + "input_files/human_3UTRsequences.txt";

		String mclOutputFile = wd + "benchmark/";
		
		String outputFilePrefix = wd + "benchmark/fasta-mcli2-dec2024/corrNet2-400_sequences_fasta_MEME_mclCluster_";

		/* Load MCL output file, line-by-line (one line = one cluster) */ 
		System.out.println("**Get sequences for MEME analysis**");
		try {
			InputStream in = new FileInputStream(new File(mclOutputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // no header
			int lineCount = 1;
			while(line !=null) {
				System.out.println(lineCount);
				HashSet<String> proteins = new HashSet<>(Arrays.asList(line.split("\t"))); // set of proteins in MCL cluster
				System.out.println("cluster contains x proteins : " + proteins.size());
				if(proteins.size() >= 3) {
					/* Load refSeqIds associated to proteins in cluster */
					List<HashSet<String>> refSeqIDs = getRefSeqIdsCorrespondingToProteinsFromNetworkInfoFile(proteins, proteinInfoFile);
					System.out.println("Found corresponding refseqIds for x proteins : " + refSeqIDs.size());
					
					// Determine sequences to print
					System.out.println("Selecting refSeqIds to use:");
					HashMap<String, Boolean	> idsToPrintMap = selectSequencesToPrint(fastaFile, refSeqIDs);

					if(idsToPrintMap.size() >= 3) {
						// Iterate through FASTA file, printing lines if corresponds to refSeqId in network
						String outputFile = outputFilePrefix + lineCount+".txt";
						
						System.out.println("Printing sequences: " + idsToPrintMap.size());
						printSequencesInNetwork(fastaFile, idsToPrintMap, outputFile);
						System.out.println("Done");
					}
				}
				line = input.readLine();
				lineCount++;
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/**
	 * Obtain the refSeqIDs corresponding to the different proteins in the given MCL cluster
	 * These IDs are used in the LESMoN approach.
	 * 
	 * @param proteinsInCluster Set<String> - Set of proteins
	 * @param networkInfoFile	String - File containing proteinName \t RefSeqId1|2|..|n
	 * @return refSeqIDs		List<String[]> - List of [RefSeqId1, 2, .., n] 
	 */
	private static List<HashSet<String>> getRefSeqIdsCorrespondingToProteinsFromNetworkInfoFile(HashSet<String> proteinsInCluster, String networkInfoFile){

		List<HashSet<String>> refSeqIDList= new ArrayList<HashSet<String>>();

		try {
			InputStream in = new FileInputStream(new File(networkInfoFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			while(line !=null) {

				if(proteinsInCluster.contains(line.split("\t")[0])) { // [0] = proteinName

					if(line.split("\t").length > 1) {  // check if protein has corresponding refSeqIds (possible there are none)
						String[] refSeqIds = line.split("\t")[1].split("\\|");
						refSeqIDList.add(new HashSet<String>(Arrays.asList(refSeqIds)));
					}
				}
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return refSeqIDList;
	}

	public static HashMap<String, Boolean> selectSequencesToPrint(String fastaFile, List<HashSet<String>> refSeqIds){
		HashMap<String, Boolean> idsToPrintSet = new HashMap<>();
		
		int idCount = 0;
		for(HashSet<String> idSet : refSeqIds) { // IDs correspond to a single protein
			
			if(idSet.size() > 1) {
				/* Check if sequence exists corresponding to Id; and get sequence length */
				int maxSequenceLength=0;
				String idToPrint = "";

				try {
					InputStream in = new FileInputStream(new File(fastaFile));
					BufferedReader input = new BufferedReader(new InputStreamReader(in));

					String line = input.readLine();
					int sequenceCount = 0; 

					while(line != null || sequenceCount < idSet.size()) {

						if(line.startsWith(">")) {
							String id = line.split("\\_")[2] + "_" + line.split("\\_|\\.")[3];
							if(idSet.contains(id)) {
								sequenceCount ++; // sequence has been found

								/* calculate sequence length */ 
								String[] positions = line.split(":| |-");
								int seqLength = Integer.parseInt(positions[3]) - Integer.parseInt(positions[2]) + 1;

								if(seqLength > maxSequenceLength) {
									maxSequenceLength = seqLength;
									idToPrint = id;
								}
							}
						}
						line = input.readLine();	
					}
					input.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
				idsToPrintSet.put(idToPrint, true);
			} else {
				for(String id: idSet) {
					idsToPrintSet.put(id, true);
				}	
			}
			idCount++;
			System.out.print(idCount + ".");
			
		}
		return idsToPrintSet;
	}

	public static void printSequencesInNetwork(String fastaFile, HashMap<String, Boolean> refSeqIdSet, String outputFile) {

		try {
			InputStream in = new FileInputStream(new File(fastaFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			String line = input.readLine();
			boolean print = false; 

			while(line !=null) {

				/* check if refseq Id is in network*/
				if(line.startsWith(">")) {
					String id = line.split("\\_")[2] + "_" + line.split("\\_|\\.")[3];
					if(refSeqIdSet.containsKey(id)) {
						if(refSeqIdSet.get(id)) {
							print = true;
							refSeqIdSet.put(id, false);
						}
					} else {
						print = false;
					}
				}

				if(print) {
					out.write(line + "\n");
					out.flush();
				}

				line = input.readLine();
			}
			out.close();
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
