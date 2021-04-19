package motifs;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class MotifDegeneration {

	private int motifLength;
	private int solutions;
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
		this.degenCharacterMap = defineIUPAC(); 											 // Define nucleotide mapping to various degenerate characters (following IUPAC standard)
		this.degenCharacterSet = new HashSet<>(Arrays.asList('R', 'D', 'H', 'V', 'B', '*')); // Define set of characters that are degenerate characters
		this.solutions = calculateSolutions(); 												 // Calculate the number of motif solutions
	}

	public void enumerateDegenerateMotifs(String inputFile, String outputFile) {
		FileInputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String motif = input.readLine();
			System.out.println("Degenerating motifs: ");
			int motifCount=1;
			while(motif != null) {
				System.out.print(motifCount + ".");
				HashSet<String> degenMotifs = getAllMotifs(motif);
				printMotifs(outputFile, motif, degenMotifs);
				
				motif = input.readLine();
				motifCount++;
				if(motifCount%10 == 0) {
					System.out.println();
				}
			}
			System.out.print("Done" +"\n");
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
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
	 * For a given motif produce all possible motifs (ie. degenerate motifs) with given criteria.
	 * Current criteria includes motifs with a max of 7 degenerate characters. 
	 * Inspired by Cartesian Products: https://stackoverflow.com/questions/9591561/java-cartesian-product-of-a-list-of-lists  
	 * 
	 * @param motif				String - sequence motifs 
	 * @return degenerateMotifs ArrayList<String> - all possible degenerate sequence motifs 
	 */
	private HashSet<String> getAllMotifs(String motif) {
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
			if(degenCount <= this.maxDegenThreshold) {
				degenerateMotifs.add(seq);
			}

		}	
		return degenerateMotifs;
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
		
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(outputFile), true));
			
			out.write(motif + "\t");
			
			for(String degenMotif: degenMotifs) {
				out.write(degenMotif + "|");
				out.flush();
			}
			out.write("\n");
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	public void generateAllPossibleMotifs(String outputFile) {
		System.out.println("Generating set of degenerate motifs");
		List<Character> nucleotides = Arrays.asList('A', 'T', 'C', 'G', 'R', 'Y', 'B', 'D', 'H', 'V', '*'); 
		String accum = "";
		HashMap<String, ArrayList<Integer>> container = new HashMap<String, ArrayList<Integer>>();
		
		permutation(nucleotides, this.motifLength, accum, container, outputFile);
	//	System.out.println("\nPrinting motifs");
		printDegenMotifSet(outputFile, container);
	}
	
	private void permutation(List<Character> NA, int k, String accumulated, HashMap<String, ArrayList<Integer>> container, String outputFile){
		
		if(k == 0)
		{
			//System.out.println(accumulated);
	
			ArrayList<Integer> value = new ArrayList<Integer>();
			container.put(accumulated, value);
			
			if(container.size()%1000000 == 0) {
				printDegenMotifSet(outputFile, container);
			}
			
			return;
		}
		for(Character na:NA)
			permutation(NA, k - 1, accumulated + na, container, outputFile);
		
	}
	
	private void printDegenMotifSet(String outputFile, HashMap<String, ArrayList<Integer>> degenMotifs) {
		
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(outputFile), true));
			
			
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
/*	private HashMap<String, HashSet<String>> generateMotifsForGivenProteinInParallel(String refSeqID, int maxDegenThreshold) throws InterruptedException, ExecutionException{
		HashMap<String, HashSet<String>> seqMap = new HashMap<>();
		HashSet<String> seqMotifSet = new HashSet<>();

		 find line number where sequence begins in fasta file  
		int lineNumber = fastaIdxMap.get(refSeqID);  

		 obtain full rna sequence from fasta file 
		String rnaSequence = getRNAsequence(lineNumber);

		List<Callable<HashSet<String>>> tasks = new ArrayList<Callable<HashSet<String>>>();
		 Use sliding window to obtain all possible motifs  
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
		 produce all degenerate motifs and add to hashSet  
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
	}*/

	/*private HashMap<String, HashSet<String>> generateMotifsForGivenProteinInParallelWithFutures(String refSeqID, int maxDegenThreshold) throws InterruptedException, ExecutionException{
		HashMap<String, HashSet<String>> seqMap = new HashMap<>();
		HashSet<String> seqMotifSet = new HashSet<>();

		 find line number where sequence begins in fasta file  
		int lineNumber = fastaIdxMap.get(refSeqID);  

		 obtain full rna sequence from fasta file 
		String rnaSequence = getRNAsequence(lineNumber);
		List<Future<HashSet<String>>> futureList = new ArrayList<>();
		ExecutorService exec = Executors.newFixedThreadPool(10);
		 Use sliding window to obtain all possible motifs  
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
		 produce all degenerate motifs and add to hashSet  
		
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
	}*/
	
/*	private HashMap<String, HashSet<String>> generateMotifsForGivenProtein(String refSeqID, int maxDegenThreshold){
		HashMap<String, HashSet<String>> seqMap = new HashMap<>();
		HashSet<String> seqMotifSet = new HashSet<>();

		 find line number where sequence begins in fasta file  
		int lineNumber = fastaIdxMap.get(refSeqID);  

		 obtain full rna sequence from fasta file 
		String rnaSequence = getRNAsequence(lineNumber);
		System.out.println("Motifs to test: " + (rnaSequence.length()-motifLength));
		 Use sliding window to obtain all possible motifs  
		for(int i=0; i<rnaSequence.length()-motifLength; i++) {
			long start = System.currentTimeMillis();
			String motif = rnaSequence.substring(i, i+motifLength); // substring keeps characters from start index (inclusive) to end index (exclusive)
			HashSet<String> motifList = getAllMotifs(motif, maxDegenThreshold);
			 produce all degenerate motifs and add to hashSet  
			seqMotifSet.addAll(motifList);
			long elapsed = System.currentTimeMillis() - start;
			System.out.println(String.format("Eslapsed time Motif " + i+ ": %d ms | motifSet size: " + seqMotifSet.size(), elapsed));
		}

		seqMap.put(refSeqID, seqMotifSet);
		return seqMap;
	}*/

	
	/**
	 * Enumerate motifs for all 3'UTR sequences associated to a protein within the BioID network
	 * 
	 * @param maxDegenThreshold		Integer - max number of degenerate characters allowed in the motifs
	 * @return generatedMotifs		Map<String, HashSet> - map of {RefSeqID : list of motif} 
	 * 
	 * @throws InterruptedException
	 * @throws ExecutionException
	 */
/*	public HashMap<String, HashSet<String>> enumerateMotifsInParallel (int maxDegenThreshold) throws InterruptedException, ExecutionException {

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
	}*/

}
