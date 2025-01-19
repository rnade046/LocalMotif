package annotateMotifs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class MotifDegeneration {

	private int motifLength;
	private int maxDegenThreshold;
	private HashSet<Character> degenCharacterSet;
	
	/**
	 * Constructor for motif enumeration class 
	 * 
	 * @param fastaFile_				String - fasta file containing all rna sequence
	 * @param proteinsInNetworkList_	ArrayList<Protein> - list of proteins in network
	 * @param motifLength_				Int - motif length of interst
	 */
	public MotifDegeneration(int motifLength_, int maxDegenThreshold_) {
		this.motifLength = motifLength_;
		this.maxDegenThreshold = maxDegenThreshold_;
		this.degenCharacterSet = new HashSet<>(Arrays.asList('R', 'D', 'H', 'V', 'B', '*')); // Define set of characters that are degenerate characters
	}

	public void generateAllPossibleMotifs(String outputFilePrefix) {

		System.out.println("Generating set of degenerate motifs");
		
		List<Character> nucleotides = Arrays.asList('A', 'T', 'C', 'G', 'R', 'Y', 'B', 'D', 'H', 'V', '*'); 
		
		String accum = "";
		HashMap<String, ArrayList<Integer>> container = new HashMap<String, ArrayList<Integer>>();
		
		int breaks = (int) Math.round(Math.pow(nucleotides.size(), this.motifLength))/999 ;
		Integer fileIdx = 0;

		fileIdx = permutation(nucleotides, this.motifLength, accum, container, breaks, fileIdx, outputFilePrefix);
		printDegenMotifSet(outputFilePrefix, fileIdx, container);
	}

	private Integer permutation(List<Character> NA, int k, String accumulated, HashMap<String, ArrayList<Integer>> container, int breaks, Integer fileIndex, String outputFilePrefix){

		if(k == 0)
		{
			//System.out.println(accumulated);

			ArrayList<Integer> value = new ArrayList<Integer>();
			container.put(accumulated, value);

			if(container.size()%breaks == 0) {
				printDegenMotifSet(outputFilePrefix, fileIndex, container);
				fileIndex++;
			}

			return fileIndex;
		}
		for(Character na:NA)
			fileIndex = permutation(NA, k - 1, accumulated + na, container, breaks, fileIndex, outputFilePrefix);

		return fileIndex;

	}

	private void printDegenMotifSet(String outputFile, int fileIdx, HashMap<String, ArrayList<Integer>> degenMotifs) {

		System.out.print(fileIdx + ".");
		if(fileIdx%50==0 && fileIdx !=0) {
			System.out.println();
		}

		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(outputFile + fileIdx)));


			for(String motif: degenMotifs.keySet()) {

				int degenCount = 0;
				/* generate a given motif */
				for(int k=0; k < motif.length() ; k++) {

					/* keep count if added character is a generate character */
					if(degenCharacterSet.contains(motif.charAt(k))) {
						degenCount++;
					}
				}

				/* motif must pass max number of degenerate character threshold to be considered in our approach */
				if(degenCount <= this.maxDegenThreshold) {
					out.write(motif + "\n");
					out.flush();
				}

			}
			degenMotifs.clear();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
