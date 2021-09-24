package ClusteredMotifs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;

import javax.sound.midi.Sequence;

public class PositionConservation {

	private String fasta;
	private int motifLength;
	private HashMap<String, Integer> fastaIdx;
	private HashMap<String, String[]> proteinInfoMap;

	// Constructor
	public PositionConservation(String fastaFile, String proteinInfoFile, int l) {

		this.fasta = fastaFile;
		this.motifLength = l;
		System.out.println("Generate fasta file index");
		this.fastaIdx = generateFastaFileIdx(fastaFile); // generate index for FASTA file
		System.out.println("Load protein info file");
		this.proteinInfoMap = loadProteinInfoFile(proteinInfoFile); // load protein info file
	}

	public void getMotifPositions(String extractedAnnotationsFile) {

		/* Load significant motif and its annotated proteins from extracted annotation 1 at a time */ 
		InputStream in;
		try {
			in = new FileInputStream(new File(extractedAnnotationsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {

				String[] col = line.split("\t"); // [0] = motif, [1] = proteinList
				
				/* For every protein get list of refseqIDs*/
				for(String prot: col[1].split("\\|")) {
					
					if(this.proteinInfoMap.containsKey(prot)) {
					
						String[] refSeqIDs = this.proteinInfoMap.get(prot);
						
						for(String id: refSeqIDs) {
							if(this.fastaIdx.containsKey(id)) {
								int indexOfFasta = this.fastaIdx.get(id);
								
								/* Load sequence */
								String sequence = loadRNAsequence(indexOfFasta);
								
								/* Search for motif in sequence */
								
							}
						}
						
					}
					
				
				}
				/* Search for motif in reference sequence */
				
				
				line = input.readLine();
			}
			input.close();
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
	private static HashMap<String, Integer> generateFastaFileIdx(String fastaFile){

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

	private static HashMap<String, String[]> loadProteinInfoFile(String proteinInfoFile){

		HashMap<String, String[]> proteinToRefSeqIDMap = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(proteinInfoFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {

				String[] col = line.split("\t");
				proteinToRefSeqIDMap.put(col[0], col[1].split("\\|")); // [0] = motif, [1] = list of refSeqIDs
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinToRefSeqIDMap;
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
	private String loadRNAsequence(int lineNumber) {
		String sequence = "";

		InputStream in;
		try {
			in = new FileInputStream(new File(this.fasta));
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
	
	private int searchForMotif(String motif, String sequence){
	
		int motifPosition = 0;
		HashSet<String> possibleInstances = new HashSet<>();
		
		
		int currentPosition = 0;
		for(int i=(sequence.length() - this.motifLength); i>=0; i--) {
			
			
			
			currentPosition++;
		}
		
		return motifPosition;
	}
	
	private HashSet<String> getPossibleMotifInstances(String motif, String sequence){
		
	}
}
