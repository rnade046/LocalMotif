package randomSeq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;

public class RandomizeSequences {


	public static void generateRandomizedFasta(String refSeqIdFile, String fastaFile, String randomizedFastaFile) {

		/* Load IDs associated to fasta Seq to randomize */
		HashSet<String> refSeqIdSet = loadRefSeqIDs(refSeqIdFile);
		HashMap<String, Integer> idxOfFasta = generateIndexOfFastaFile(fastaFile);
		System.out.println("Number of seq to randomize: " +refSeqIdSet.size());
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(randomizedFastaFile)));
			int countId = 0;

			for(String refSeqId: refSeqIdSet) {			
				countId++;

				if(countId%10 == 0) {
					System.out.print(countId + ".");
				} 
				if(countId%100 == 0) {
					System.out.println();
				}
				/* Get formatted sequence as String */
				if(idxOfFasta.containsKey(refSeqId)) {

					String seq = getRNAsequence(fastaFile, idxOfFasta.get(refSeqId));
					String randomSeq = randomizeSequence(seq);

					out.write(">RANDOM_" + refSeqId + "\n");
					out.write(randomSeq + "\n\n");
					out.flush();
				} 
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Load refseq ids that match to a protein in the network.
	 * 
	 * @param refSeqIdsFile	String - file path: list of refseq ids (1 per line)
	 * @return refSeqList	HashSet<String> set of refSeqIds
	 */
	private static HashSet<String> loadRefSeqIDs(String refSeqIdsFile){
		HashSet<String> refSeqList = new HashSet<String>();

		InputStream in;
		try {
			in = new FileInputStream(new File(refSeqIdsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			while(line!=null) {
				refSeqList.add(line);
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return refSeqList;
	}

	/**
	 * Generate an index for the FASTA file; identifying the line number for the different refseq identifiers
	 * 
	 * @param fastaFile String - file path for the FASTA sequences
	 * @return indexOfFastaFile HashMap<String, Integer> - map of {refseqId : line count} 
	 */
	private static HashMap<String, Integer> generateIndexOfFastaFile(String fastaFile){
		HashMap<String, Integer> indexOfFastaFile = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(fastaFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			int lineCount = 1;

			while(line != null) {
				/* store the line index of the start of a new sequence */
				if(line.startsWith(">")) {
					String[] col = line.split("_|\\.");
					String refSeqId = col[2] + "_" + col[3]; // col[2] = type of ID (eg. NM) ; col[3] = number ID (#####)

					indexOfFastaFile.put(refSeqId, lineCount);
				}

				line = input.readLine();
				lineCount++;
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return indexOfFastaFile;
	}

	/**
	 * Obtains full RNA seqeunce from the fasta file given the line number corresponding to the sequence header. 
	 * Will store all information following sequence header until the next sequence header 
	 * 
	 * @param fastaFile		String - file path to fasta file containing all rna sequences
	 * @param lineNumber	int - corresponds to the line number from which to start reading the file
	 * 
	 * @return sequence		String - full rna sequence 
	 */
	private static String getRNAsequence(String fastaFile, int lineNumber) {
		String sequence = "";

		InputStream in;
		try {
			in = new FileInputStream(new File(fastaFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // first line

			/* read all lines until the line number specified */
			for(int i=1; i<=lineNumber; i++) {
				line = input.readLine();
			}

			/* store line info until we reach the end of the file or the next sequence header */
			while(line !=null && !line.startsWith(">")) {
				sequence += line;
				line = input.readLine();
			}
			input.close();

		}catch (IOException e) {
			e.printStackTrace();
		}
		return sequence;
	}

	/**
	 * For a given sequence, generate a it's random sequence using non overlapping window
	 * of 10 character 
	 * 
	 * @param seq	String - sequence to randomize
	 * @return randomSeq 	String - randomized sequence
	 */
	private static String randomizeSequence(String seq) {

		String randomSeq = "";

		/* iterate over sequence with a non overlaping window of size 10 */
		for(int i=0; i<seq.length(); i+=10) {

			String substring = "";

			/* substring of length 10, unless remainder sequences contains less then 10 character */ 
			if((seq.length() - i) >= 10) {
				substring = seq.substring(i, (i+10)); // returns substring of length 10
			} else {
				substring = seq.substring(i); // returns remainder of sequence (length < 10) 
			}

			if(substring.length() >= 2) {

				int max = substring.length();
				int min = 0;

				/* pairwise switch 1000 times two characters */ 
				for(int j=0; j<1000; j++) {
					String currentSubstring = "";

					/* initialize indexed of characters to switch, cannot be the same index or map to the same characters */
					int idx1, idx2 = 0;
					do {
						idx1 = (int) ((Math.random() * (max - min)) + min);
						idx2 = (int) ((Math.random() * (max - min)) + min);
					} while(idx1 == idx2);

					/* generate new substring, switching characters at selected indexes */
					for(int k=0; k<substring.length(); k++) {
						if(k == idx1) {
							currentSubstring += substring.charAt(idx2);
						} else if (k == idx2) {
							currentSubstring += substring.charAt(idx1);
						} else {
							currentSubstring += substring.charAt(k);
						}
					}
					substring = currentSubstring; // re-initialize substring
				}
			}

			randomSeq += substring; // add substring to seq  
		}
		return randomSeq;
	}
}