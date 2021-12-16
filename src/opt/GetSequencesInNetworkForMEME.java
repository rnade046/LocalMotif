package opt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class GetSequencesInNetworkForMEME {

	public static void main(String[] args) {

		String wd = "C://Users//Rachel//Documents//LESMoNlocal//analysis//";
		String proteinInfoFile = wd + "corrNetTop2-400_proteinsInNetwork_info.tsv";
		String fastaFile = wd + "input_files//human_3UTRsequences.txt";

		String mclOutputFile = wd + "benchmark//mclOutput_i2.txt";

		// Load MCL output file, line-by-line (one line = one cluster)
		System.out.println("**Get sequences for MEME analysis**");
		try {
			InputStream in = new FileInputStream(new File(mclOutputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // no header
			int lineCount = 1;
			while(line !=null) {
				System.out.println(lineCount);
				HashSet<String> proteins = new HashSet<>(Arrays.asList(line.split("\t")));
				if(proteins.size() >= 3) {
					// Load refSeqIds associated to proteins in cluster
					HashMap<String, Boolean> refSeqIdSet = loadRefSeqIdsFromNetworkInfoFile(proteins, proteinInfoFile);

					// Iterate through FASTA file, printing lines if corresponds to refSeqId in network
					String outputFile = wd + "benchmark//fasta//corrNet2-400_sequences_fasta_MEME_mclCluster_"+ lineCount+".txt";

					printSequencesInNetwork(fastaFile, refSeqIdSet, outputFile);
				}
				line = input.readLine();
				lineCount++;
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static HashMap<String, Boolean> loadRefSeqIdsFromNetworkInfoFile(HashSet<String> proteinsInCluster, String networkInfoFile){

		HashMap<String, Boolean> refSeqIDSet = new HashMap<>();

		try {
			InputStream in = new FileInputStream(new File(networkInfoFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			while(line !=null) {

				if(proteinsInCluster.contains(line.split("\t")[0])) {
					if(line.split("\t").length > 1) {
						String[] refSeqIds = line.split("\t")[1].split("\\|");
						for(String id : refSeqIds) {
							refSeqIDSet.put(id, true);
						}

					}
				}
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return refSeqIDSet;
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
