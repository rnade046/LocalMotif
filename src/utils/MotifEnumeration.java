package utils;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import graph.*;

public class MotifEnumeration {

	public static void enumerateMotifs (String fastaFile, ArrayList<Protein> proteinsInNetworkList, int motifLength) {
		
		/* Define nucleotide mapping to various degenerate characters (following IUPAC standard) */
		HashMap<Character, Character[]> degenCharacterMap = defineIUPAC();
		
		/* Define set of characters that are degenerate characters */
		HashSet<Character> degenCharacterSet = new HashSet<>(Arrays.asList('R', 'D', 'H', 'V', 'B', '*'));
		
		/* Calculate the number of motif solutions */
		int solutions = calculateSolutions(motifLength, degenCharacterMap);
		
		/* Generate index for the fasta file */
		HashMap<String, Integer> fastaIndex = generateIndexOfFastaFile(fastaFile);
	}
	
	/**
	 * Define the RNA nucleotide mapping to their possible substitutions (ie. degenerate characters)
	 * Note: for now this is fix - this could be passed as an input parameter in the future to enable 
	 * flexible motif representation for users
	 * Maps RNA nucleotides {A, U, C, G} to possible substitutions (themselves + degenerate characters) 
	 * 
	 * @return charMap	 HashMap<Character, Character[]> - map {nucleotide : list of possible characters}
	 */
	private static HashMap<Character, Character[]> defineIUPAC(){
		HashMap<Character, Character[]> charMap = new HashMap<>();
		
		Character[] a = {'A', 'R', 'D', 'H', 'V', '*'};
		Character[] u = {'U', 'Y', 'B', 'D', 'H', '*'};
		Character[] g = {'G', 'R', 'B', 'D', 'V', '*'};
		Character[] c = {'C', 'Y', 'B', 'H', 'V', '*'};
		
		charMap.put('A', a);
		charMap.put('U', u);
		charMap.put('G', g);
		charMap.put('C', c);
		
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
	private static int calculateSolutions(int motifLength, HashMap<Character, Character[]> charMap) {
		int solutions = 1;
		for(int i = 0; i < motifLength; i++) {
	    	solutions *= charMap.get('A').length;
	    }
		return solutions;
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
	
	private static HashMap<String, HashSet<String>> generateMotifsForGivenProtein(String refSeqID){
		HashMap<String, HashSet<String>> seqMotifMap = new HashMap<>();
		
		String rnaSequence = "";
		
		return seqMotifMap;
	}
	
	private static String getRNAsequence(String fastaFile, HashMap<String, Integer> indexOfFastaMap) {
		String sequence = "";
		
		
		return sequence;
	}
}
