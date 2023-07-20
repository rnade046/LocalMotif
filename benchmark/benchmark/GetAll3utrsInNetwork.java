package benchmark;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashSet;

public class GetAll3utrsInNetwork {

	public static void main(String[] args) {

		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";

		String proteinInfoFile = wd + "corrNetTop2-400_proteinsInNetwork_info.tsv";
		String utrFastaFile = wd + "input_files/human_3UTRsequences.txt";

		String outputFile = wd + "benchmark/corrNet2-400_3utr_allSeqsInNetwork.fasta";

		/* get all refSeqIds associated to proteins in the network */
		HashSet<String> refSeqIdSet = getRefSeqIDsInNetwork(proteinInfoFile);
		System.out.println("Number of ids: " + refSeqIdSet.size());
		
		try {
			InputStream in = new FileInputStream(new File(utrFastaFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			String line = input.readLine();
			boolean print = false; 
			int seqCount = 0;

			while(line !=null) {

				/* check if RefSeq Id is in network*/
				if(line.startsWith(">")) {
					String id = line.split("\\_")[2] + "_" + line.split("\\_|\\.")[3];
					if(refSeqIdSet.contains(id)) {
						print = true;
						seqCount++;
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
			
			System.out.println("Number of sequences: " + seqCount);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static HashSet<String> getRefSeqIDsInNetwork(String proteinInfoFile){

		HashSet<String> refSeqIdSet = new HashSet<>();

		BufferedReader in;
		try {
			in = new BufferedReader(new InputStreamReader(new FileInputStream(proteinInfoFile)));

			String line = in.readLine();

			while(line!=null) {

				if(line.split("\t").length > 1) {
					refSeqIdSet.addAll(Arrays.asList(line.split("\t")[1].split("\\|")));
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return refSeqIdSet;
	}
}
