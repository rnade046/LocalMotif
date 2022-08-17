package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;


public class SequenceLength {

	public static void main(String[] args) {
		
		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";
		String fastaFile = wd + "/input_files/human_3UTRsequences.txt";
		
		String refSeqIdFile = wd + "corrNetTop2_proteinsInNetwork_info.tsv";
		String outptutFile = wd + "nucleotide-composition/corrNetTop2_sequenceLengths_allSeqsConsidered.tsv";
		
		System.out.println("* Assess - all sequences *");
		calculateNucleotideComp(refSeqIdFile, fastaFile, outptutFile);
		
		String refSeqIdLongestSeqsFile = wd + "/MotifPosition/corrNetTop2_longestSequence.tsv";
		String outputLongestSeqsFile = wd + "nucleotide-composition/corrNetTop2_sequenceLengths_LongestSeqs.tsv";
		
		System.out.println("* Assess - longest sequences *");
		calculateNucleotideComp(refSeqIdLongestSeqsFile, fastaFile, outputLongestSeqsFile);
	}

	private static void calculateNucleotideComp(String refSeqIdsFile, String fastaFile, String outputFile) {

		/* load RefseqIds */
		HashSet<String> idSet = loadRefSeqIds(refSeqIdsFile);

		/* determine nucleotide composition from FASTA file */
		List<Integer> seqLengths = countSequenceLength(fastaFile, idSet);
		
		/* output composition to file */
		printCounts(seqLengths, outputFile);
		
	}

	private static HashSet<String> loadRefSeqIds(String refSeqIdFile){

		HashSet<String> idSet = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(refSeqIdFile))));

			String line = in.readLine();

			while(line != null) {
				
				if(line.split("\t").length > 1 ) {
					String[] ids = line.split("\t")[1].split("\\|"); // [1] = list of IDs {id1|id2|..|idx} 
					idSet.addAll(Arrays.asList(ids));
				}

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return idSet;
	}
	
	private static List<Integer> countSequenceLength(String fastaFile, HashSet<String> ids){
		
		List<Integer> seqLengths = new ArrayList<>();
		
		/* read through FASTA file */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));			

			String line = in.readLine();

			boolean readSeq = false; 
			String seq = "";

			while(line!=null) {

				if(readSeq & !line.startsWith(">")) {
					seq += line;
				}

				/* check if sequence is associated to a protein in network */
				if(line.startsWith(">")) {
					readSeq = false;

					/* if sequence has been loaded - assess nucleotide composition  */
					if(!seq.isEmpty()) {
						seqLengths.add(seq.length());
					}

					//  >hg38_ncbiRefSeq_NM_001276352.2 range=chr1:67092165-67093579 5'pad=0 3'pad=0 strand=- repeatMasking=none					
					String id = line.split("[\\_\\s++\\.]")[2] + "_"+ line.split("[\\_\\s++\\.]")[3];

					if(ids.contains(id)) {
						readSeq = true;
						seq = "";
					}
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return seqLengths;
	}
	
	private static void printCounts(List<Integer> seqCounts, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(int i=0; i<seqCounts.size(); i++) {
				out.write(seqCounts.get(i) + "\n");
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
