package precompileMotifs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.zip.DeflaterOutputStream;

public class MotifDegeneration {

	private int motifLength;
	private int maxDegenThreshold;
	private HashMap<Character, Character[]> degenCharacterMap;
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
		this.degenCharacterMap = defineDegenCharacterMap(); 											 // Define nucleotide mapping to various degenerate characters (following IUPAC standard)
		this.degenCharacterSet = new HashSet<>(Arrays.asList('R', 'D', 'H', 'V', 'B', '*')); // Define set of characters that are degenerate characters
	}

	public void enumerateNonDegenerateMotifs(String inputFilePrefix, int i, String outputFilePrefix) {

//		int fileCount = dir.list().length;
//
//		for(int i=0; i<fileCount; i++) {

			FileInputStream in;
			try {
				in = new FileInputStream(new File(inputFilePrefix + i));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String degenMotif = input.readLine();
				System.out.println("Degenerating motifs: ");
				int motifCount=1;
				while(degenMotif != null) {
					if(motifCount%1000 == 0){
						System.out.print(motifCount + ".");
					}

					int solutions = calculateSolutions(degenMotif);
					HashSet<String> nonDegenMotifSet = getAllMotifs(degenMotif, solutions);
					printMotifs(outputFilePrefix + i, degenMotif, nonDegenMotifSet);

					degenMotif = input.readLine();
					motifCount++;

					if(motifCount%10000 == 0) {
						System.out.println();
					}
				}
				System.out.print("Done" +"\n");
				input.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
//		}

	}
	/**
	 * Define the RNA nucleotide mapping to their possible substitutions (ie. degenerate characters)
	 * Note: for now this is fix - this could be passed as an input parameter in the future to enable 
	 * flexible motif representation for users
	 * Maps RNA nucleotides {A, U, C, G} to possible substitutions (themselves + degenerate characters) 
	 * 
	 * @return charMap	 HashMap<Character, Character[]> - map {nucleotide : list of possible characters}
	 */
	private HashMap<Character, Character[]> defineDegenCharacterMap(){
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
	private int calculateSolutions(String motif) {
		int solutions = 1;
		for(int i = 0; i < motifLength; i++) {
			solutions *= degenCharacterMap.get(motif.charAt(i)).length;
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
	private HashSet<String> getAllMotifs(String motif, int solutions) {
		HashSet<String> motifs = new HashSet<>();

		/* generate all solutions */ 
		for(int i = 0; i < solutions; i++) {
			int j = 1;
			String seq = "";
			/* generate a given motif */
			for(int k=0; k < motif.length() ; k++) {
				Character[] set = degenCharacterMap.get(motif.charAt(k));

				char charToAppend = set[(i/j)%set.length]; // magic 
				seq += charToAppend;
				j *= set.length;
			}

			motifs.add(seq);

		}	
		return motifs;
	}

	/**
	 * Print degenerate motifs associated to the real motif. Lines are appended to existing file
	 * No header. Example line: 
	 * Motif \t Degen1|Degen2|..|DegenX \n
	 * 
	 * @param outputFile	String - file path to generated degenerate motifs
	 * @param motif			String - motif for which degenerate motifs were generated
	 * @param degenMotifs	HashSet<String> - Set of degenerated motifs 
	 */
	private static void printMotifs(String outputFile, String motif, HashSet<String> degenMotifs) {

//		BufferedWriter out;
//		try {
//			out = new BufferedWriter(new FileWriter(new File(outputFile), true));
//
//			out.write(motif + "\t");
//
//			for(String degenMotif: degenMotifs) {
//				out.write(degenMotif + "|");
//				out.flush();
//			}
//			out.write("\n");
//			out.flush();
//			out.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		
		try {
			DeflaterOutputStream out = new DeflaterOutputStream(new FileOutputStream(new File(outputFile), true));

			String line = motif + "\t";
			out.write(line.getBytes());
			
			for(String degenMotif: degenMotifs) {
				
				line = degenMotif + "|";
				out.write(line.getBytes());
				out.flush();
			}
			line = "\n";
			out.write(line.getBytes());
			out.flush();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void generateAllPossibleMotifs(String outputFilePrefix) {
		System.out.println("Generating set of degenerate motifs");
		List<Character> nucleotides = Arrays.asList('A', 'T', 'C', 'G', 'R', 'Y', 'B', 'D', 'H', 'V', '*'); 
		String accum = "";
		HashMap<String, ArrayList<Integer>> container = new HashMap<String, ArrayList<Integer>>();
		int breaks = (int) Math.round(Math.pow(nucleotides.size(), this.motifLength))/999 ;
		Integer fileIdx = 0;

		fileIdx = permutation(nucleotides, this.motifLength, accum, container, breaks, fileIdx, outputFilePrefix);
		//	System.out.println("\nPrinting motifs");
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
