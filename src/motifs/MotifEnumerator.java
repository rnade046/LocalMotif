package motifs;

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

public class MotifEnumerator {
	private String fastaFile;
	private int motifLength;
	private HashMap<String, Integer> fastaIdxMap;
	private HashSet<String> refSeqIds;
	
	public MotifEnumerator(String _fastaFile, int _motifLength, HashSet<String> _refSeqIds) {
		this.fastaFile = _fastaFile;
		this.motifLength = _motifLength;
		this.refSeqIds = _refSeqIds;
		this.fastaIdxMap = generateIndexOfFastaFile(_fastaFile);
	}

	public void generateMotifList(String mapMotifsToRefSeqIdsFile, String motifsToTestFile){
		
		HashSet<String> motifSet = new HashSet<>(); // set to contain all motifs
		int countMissingSeq = 0;	// count of missing refseq ids in fasta 
		int countAllMotifs = 0;
		
		int idCount = 1;
		System.out.println("Testing id: ");
		
		for(String id: refSeqIds) {
			if(idCount%100 == 0) {
				System.out.println();
			}
			System.out.print(idCount + "|");
			idCount++;
			
			String formattedSeq = "";
			
			/* Get formatted sequence as String */
			if(fastaIdxMap.containsKey(id)) {
				formattedSeq = getRNAsequence(fastaIdxMap.get(id));
			} else {
				countMissingSeq++;
				continue; 	// do not test id if it's sequence is not in the fasta file
			}
			
			/* Enumerate motifs from formatted sequence */
			HashSet<String> currentMotifSet = generateMotifsFromFormattedSequence(formattedSeq);
			countAllMotifs += currentMotifSet.size();
			
			/* print RefseqId = motifs to file */ 
			printRefSeqMotifs(id, currentMotifSet, mapMotifsToRefSeqIdsFile);
			
			/* Add motifs to ongoing motif Set */
			motifSet.addAll(currentMotifSet);
			
		}
		System.out.println("Number of total motifs: " + countAllMotifs);
		System.out.println("Number of unique motifs: " + motifSet.size() + "\n");
	
		System.out.println("Number of refSeq Ids not found in Fasta file : " + countMissingSeq);
		
		// print motifs to test
		printAllMotifs(motifSet, motifsToTestFile);
		
		//return motifSet;
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
	private String getRNAsequence(int lineNumber) {
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
	 * Given a formatted sequence list all motifs of length n (as specified by the motifLength parameter) using a sliding window.
	 * 
	 * @param sequence	String - formatted sequence (no spaces or special characters)
	 * @return seqMotifSet	HashSet<String> - Set of all motifs associated to the current sequence / mRNA
	 */
	private HashSet<String> generateMotifsFromFormattedSequence(String sequence){
		HashSet<String> seqMotifSet = new HashSet<>();

		/* Use sliding window to obtain all possible motifs */ 
		for(int i=0; i<sequence.length()-motifLength; i++) {
			String motif = sequence.substring(i, i+motifLength); // substring keeps characters from start index (inclusive) to end index (exclusive)
			seqMotifSet.add(motif);
		}
		return seqMotifSet;
	}

	private void printRefSeqMotifs(String id, HashSet<String> motifSet, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile), true));
			
			out.write(id + "\t");
			for(String motif : motifSet) {
				out.write(motif + "|");
			}
			out.write("\n");
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private void printAllMotifs(HashSet<String> motifSet, String outputFile) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(String motif: motifSet) {
				out.write(motif + "\n");
				out.flush();
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	 /**
	 * Generate an index for the FASTA file; identifying the line number for the different refseq identifiers
	 * 
	 * @param fastaFile String - file path for the FASTA sequences
	 * @return indexOfFastaFile HashMap<String, Integer> - map of {refseqId : line count} 
	 */
	private HashMap<String, Integer> generateIndexOfFastaFile(String fastaFile){
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


}
