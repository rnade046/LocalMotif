package motifs;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import graph.*;

public class MotifDegeneration {

	private String fastaFile;
	private ArrayList<Protein> proteinList;
	private int motifLength;
	private int solutions;
	private HashMap<Character, Character[]> degenCharacterMap;
	private HashMap<String, Integer> fastaIdxMap;
	private HashSet<Character> degenCharacterSet;
	/**
	 * Constructor for motif enumeration class 
	 * 
	 * @param fastaFile_				String - fasta file containing all rna sequence
	 * @param proteinsInNetworkList_	ArrayList<Protein> - list of proteins in network
	 * @param motifLength_				Int - motif length of interst
	 */
	public MotifDegeneration(String fastaFile_, ArrayList<Protein> proteinsInNetworkList_, int motifLength_) {
		this.fastaFile = fastaFile_;
		this.proteinList = proteinsInNetworkList_;
		this.motifLength = motifLength_;

		this.degenCharacterMap = defineIUPAC(); 											 // Define nucleotide mapping to various degenerate characters (following IUPAC standard)
		this.degenCharacterSet = new HashSet<>(Arrays.asList('R', 'D', 'H', 'V', 'B', '*')); // Define set of characters that are degenerate characters
		this.solutions = calculateSolutions(); 												 // Calculate the number of motif solutions

		this.fastaIdxMap = generateIndexOfFastaFile(); // Generate index for the fasta file
	}

	/**
	 * Enumerate motifs for all 3'UTR sequences associated to a protein within the BioID network
	 * 
	 * @param maxDegenThreshold		Integer - max number of degenerate characters allowed in the motifs
	 * @return generatedMotifs		Map<String, HashSet> - map of {RefSeqID : list of motif} 
	 * 
	 * @throws InterruptedException
	 * @throws ExecutionException
	 */
	public HashMap<String, HashSet<String>> enumerateMotifsInParallel (int maxDegenThreshold) throws InterruptedException, ExecutionException {

		HashMap<String, HashSet<String>> generatedMotifs = new HashMap<>();
		List<Callable<HashMap<String, HashSet<String>>>> tasks = new ArrayList<Callable<HashMap<String, HashSet<String>>>>();

		for(Protein prot : proteinList) {
			List<String> refSeqList = prot.getProteinId();
			for(String refSeqID: refSeqList) {
				Callable<HashMap<String, HashSet<String>>> c = new Callable<HashMap<String, HashSet<String>>>() {
					@Override
					public HashMap<String, HashSet<String>> call() throws Exception {
						return generateMotifsForGivenProtein(refSeqID, maxDegenThreshold);
					}
				};
				tasks.add(c);

			}

		}
		ExecutorService EXEC = Executors.newCachedThreadPool();
		try {
			long start = System.currentTimeMillis();

			List<Future<HashMap<String, HashSet<String>>>> results = EXEC.invokeAll(tasks);
			for (Future<HashMap<String, HashSet<String>>> fr : results) {
				generatedMotifs.putAll(fr.get());

			}
			long elapsed = System.currentTimeMillis() - start;
			System.out.println(String.format("Eslapsed time: %d ms", elapsed));

		}  finally {
			EXEC.shutdown();
		}
		return generatedMotifs;
	}


	public HashMap<String, HashSet<String>> enumerateMotifs (int maxDegenThreshold) throws InterruptedException, ExecutionException {

		HashMap<String, HashSet<String>> generatedMotifs = new HashMap<>();

		for(Protein prot : proteinList) {
			List<String> refSeqList = prot.getProteinId();
			for(String refSeqID: refSeqList) {

				System.out.println("Getting motifs for seq: " + refSeqID);
				Timestamp timestamp = new Timestamp(System.currentTimeMillis());
				System.out.println(timestamp + " [start]");
				HashMap<String, HashSet<String>> motifMap =  generateMotifsForGivenProteinInParallelWithFutures(refSeqID, maxDegenThreshold);
				timestamp = new Timestamp(System.currentTimeMillis());
				System.out.println(timestamp + " [end]" + "\n");
				generatedMotifs.putAll(motifMap);
			}

		}
		return generatedMotifs;
	}
	/**
	 * Define the RNA nucleotide mapping to their possible substitutions (ie. degenerate characters)
	 * Note: for now this is fix - this could be passed as an input parameter in the future to enable 
	 * flexible motif representation for users
	 * Maps RNA nucleotides {A, U, C, G} to possible substitutions (themselves + degenerate characters) 
	 * 
	 * @return charMap	 HashMap<Character, Character[]> - map {nucleotide : list of possible characters}
	 */
	private HashMap<Character, Character[]> defineIUPAC(){
		HashMap<Character, Character[]> charMap = new HashMap<>();

		Character[] a = {'A', 'R', 'D', 'H', 'V', '*'};
		Character[] u = {'U', 'Y', 'B', 'D', 'H', '*'};
		Character[] g = {'G', 'R', 'B', 'D', 'V', '*'};
		Character[] c = {'C', 'Y', 'B', 'H', 'V', '*'};

		charMap.put('A', a);
		charMap.put('T', u);
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
	private int calculateSolutions() {
		int solutions = 1;
		for(int i = 0; i < motifLength; i++) {
			solutions *= degenCharacterMap.get('A').length;
		}
		return solutions;
	} 

	/**
	 * Generate an index for the FASTA file; identifying the line number for the different refseq identifiers
	 * 
	 * @param fastaFile String - file path for the FASTA sequences
	 * @return indexOfFastaFile HashMap<String, Integer> - map of {refseqId : line count} 
	 */
	private HashMap<String, Integer> generateIndexOfFastaFile(){
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
	 * For a given RefSeq ID obtain 3'UTR sequence and generate all possible motifs (including degenerate motifs)
	 * Motifs can include up to the max number of degenerate characters set by maxDegenThreshold (inclusively)
	 * 
	 * @param refSeqID				String - RefSeqID identifying the required sequence
	 * @param maxDegenThreshold		Integer - max number of degenerate characters allowed in a given motif
	 * @return seqMotifSet			Set<String> - Set of degenerate sequence motif for given reference sequence
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 */
	private HashMap<String, HashSet<String>> generateMotifsForGivenProteinInParallel(String refSeqID, int maxDegenThreshold) throws InterruptedException, ExecutionException{
		HashMap<String, HashSet<String>> seqMap = new HashMap<>();
		HashSet<String> seqMotifSet = new HashSet<>();

		/* find line number where sequence begins in fasta file */ 
		int lineNumber = fastaIdxMap.get(refSeqID);  

		/* obtain full rna sequence from fasta file */
		String rnaSequence = getRNAsequence(lineNumber);

		List<Callable<HashSet<String>>> tasks = new ArrayList<Callable<HashSet<String>>>();
		/* Use sliding window to obtain all possible motifs */ 
		for(int i=0; i<rnaSequence.length()-motifLength; i++) {
		//for(int i=0; i<20; i++) {
			String motif = rnaSequence.substring(i, i+motifLength); // substring keeps characters from start index (inclusive) to end index (exclusive)
			Callable<HashSet<String>> c = new Callable<HashSet<String>>() {
				@Override
				public HashSet<String> call() throws Exception {
					return getAllMotifs(motif, maxDegenThreshold);
				}
			};
			tasks.add(c);
		}
		/* produce all degenerate motifs and add to hashSet */ 
		ExecutorService exec = Executors.newFixedThreadPool(10);
		try {
			long start = System.currentTimeMillis();

			List<Future<HashSet<String>>> results = exec.invokeAll(tasks);
			for (Future<HashSet<String>> fr : results) {
				seqMotifSet.addAll(fr.get()); 
				 System.out.println("Memory usage: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())+ " bytes");
			}
			long elapsed = System.currentTimeMillis() - start;
			System.out.println(String.format("Eslapsed time Motif: %d ms", elapsed));

		}  finally {
			exec.shutdown();
		}




		seqMap.put(refSeqID, seqMotifSet);
		return seqMap;
	}

	private HashMap<String, HashSet<String>> generateMotifsForGivenProteinInParallelWithFutures(String refSeqID, int maxDegenThreshold) throws InterruptedException, ExecutionException{
		HashMap<String, HashSet<String>> seqMap = new HashMap<>();
		HashSet<String> seqMotifSet = new HashSet<>();

		/* find line number where sequence begins in fasta file */ 
		int lineNumber = fastaIdxMap.get(refSeqID);  

		/* obtain full rna sequence from fasta file */
		String rnaSequence = getRNAsequence(lineNumber);
		List<Future<HashSet<String>>> futureList = new ArrayList<>();
		ExecutorService exec = Executors.newFixedThreadPool(10);
		/* Use sliding window to obtain all possible motifs */ 
		//for(int i=0; i<rnaSequence.length()-motifLength; i++) {
		for(int i=0; i<35; i++) {
			String motif = rnaSequence.substring(i, i+motifLength); // substring keeps characters from start index (inclusive) to end index (exclusive)
			Callable<HashSet<String>> c = new Callable<HashSet<String>>() {
				@Override
				public HashSet<String> call() throws Exception {
					return getAllMotifs(motif, maxDegenThreshold);
				}
			};
			Future<HashSet<String>> future = exec.submit(c);
		    futureList.add(future);
		}
		/* produce all degenerate motifs and add to hashSet */ 
		
		try {
			long start = System.currentTimeMillis();

			for (Future<HashSet<String>> fr : futureList) {
				seqMotifSet.addAll(fr.get()); 
				 System.out.println("Memory usage: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())+ " bytes");
			}
			long elapsed = System.currentTimeMillis() - start;
			System.out.println(String.format("Eslapsed time Motif: %d ms", elapsed));

		}  finally {
			exec.shutdown();
		}

		seqMap.put(refSeqID, seqMotifSet);
		return seqMap;
	}
	
	private HashMap<String, HashSet<String>> generateMotifsForGivenProtein(String refSeqID, int maxDegenThreshold){
		HashMap<String, HashSet<String>> seqMap = new HashMap<>();
		HashSet<String> seqMotifSet = new HashSet<>();

		/* find line number where sequence begins in fasta file */ 
		int lineNumber = fastaIdxMap.get(refSeqID);  

		/* obtain full rna sequence from fasta file */
		String rnaSequence = getRNAsequence(lineNumber);
		System.out.println("Motifs to test: " + (rnaSequence.length()-motifLength));
		/* Use sliding window to obtain all possible motifs */ 
		for(int i=0; i<rnaSequence.length()-motifLength; i++) {
			long start = System.currentTimeMillis();
			String motif = rnaSequence.substring(i, i+motifLength); // substring keeps characters from start index (inclusive) to end index (exclusive)
			HashSet<String> motifList = getAllMotifs(motif, maxDegenThreshold);
			/* produce all degenerate motifs and add to hashSet */ 
			seqMotifSet.addAll(motifList);
			long elapsed = System.currentTimeMillis() - start;
			System.out.println(String.format("Eslapsed time Motif " + i+ ": %d ms | motifSet size: " + seqMotifSet.size(), elapsed));
		}

		seqMap.put(refSeqID, seqMotifSet);
		return seqMap;
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
	 * For a given motif produce all possible motifs (ie. degenerate motifs) with given criteria.
	 * Current criteria includes motifs with a max of 7 degenerate characters.
	 * 
	 * @param motif				String - sequence motifs 
	 * @return degenerateMotifs ArrayList<String> - all possible degenerate sequence motifs 
	 */
	private HashSet<String> getAllMotifs(String motif, int maxDegenThreshold) {
		HashSet<String> degenerateMotifs = new HashSet<>();

		/* generate all solutions */ 
		for(int i = 0; i < solutions; i++) {
			int j = 1;
			String seq = "";
			int degenCount = 0;
			/* generate a given motif */
			for(int k=0; k < motif.length() ; k++) {
				Character[] set = degenCharacterMap.get(motif.charAt(k));

				char charToAppend = set[(i/j)%set.length]; // magic 
				seq += charToAppend;
				j *= set.length;

				/* keep count if added character is a generate character */
				if(degenCharacterSet.contains(charToAppend)) {
					degenCount++;
				}
			}

			/* motif must pass max number of degenerate character threshold to be considered in our approach */
			if(degenCount <= maxDegenThreshold) {
				degenerateMotifs.add(seq);
			}

		}	
		return degenerateMotifs;
	}
}
