package evolutionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class DetermineGenomePositionsOfMotifs {

	public static void main(String[] args) {

		String proteinsToTest = args[0];
		String proteinInfoFile = "corrNetTop2-400_proteinsInNetwork_info.tsv";
		String companionFile = "/home/rnade046/scratch/annotationFiles_Full_FWD/CompanionFiles_corrNetTop2-400_n20_2000/corrNetTop2-400_n20_2000_motifAnnotationsCompanionFile_";
		//String companionFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/corrNetTop2-400_coreTPD_p0.4_coreProteinsByMotif.tsv";
		String fastaFile = "human_3UTRsequences.txt";
		String outputFile = "positions/genomePositions_" + args[1];
		//String outputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/evolutionConservation/positionsByProtein.tsv";

		/* search for motifs in given annotation/companion file */
		assessSequenceWithRegex(proteinsToTest, proteinInfoFile, fastaFile, companionFile, outputFile);
	}

	public static void assessSequenceWithRegex(String proteinsToTest, String proteinInfo, String fastaFile, String companionFile, String outputFile) {

		/* load proteins to test */
		HashSet<String> proteins = loadProteinsToTest(proteinsToTest);	
		System.out.println("proteins to assess : " + proteins.size());

		/* get refSeqIds for proteins to test */ 
		HashSet<String> refSeqIds = loadProteinRefSeqIdMap(proteinInfo, proteins);
		System.out.println("sequence to assess : " + refSeqIds.size());

		/* load sequences */
		List<SearchSequence> sequences = loadSequences(refSeqIds, fastaFile);
		System.out.println("\nfound sequences : " + sequences.size());

		/* generate list of order for companion files to test */
		List<Integer> fileOrder = setOrderOfFilesToTest(1000);

		/* set character map */
		HashMap<Character, String> characterMap = setCharacterMapForRegularExpression();

		/* iterate over companion files */
		System.out.println("searching for motifs");
		for(int i=0; i<fileOrder.size(); i++) {

			if(i%50 ==0) {
				System.out.println("");
			}
			System.out.print(i + ".");
			
			/* load motifs to test */
			HashSet<String> testedMotifs = SequenceUtils.loadRepresentativeMotifs(companionFile + i);
			System.out.println(testedMotifs.size());

			for(String motif: testedMotifs) {

				/* convert motif with regular expression */
				String formattedMotif = formatMotifWithRegularExpression(motif, characterMap);
				Pattern pattern = Pattern.compile(formattedMotif);		

				/* search for motif positions */
				for(SearchSequence seq : sequences) {

					List<Integer> positionsToSearch = seq.getPositionsToSearch();
					List<Integer> positionsToRemove = new ArrayList<>();
					if(!positionsToSearch.isEmpty()) {
						Matcher matcher = pattern.matcher(seq.getSequence()).region(positionsToSearch.get(0), positionsToSearch.get(positionsToSearch.size()-1));	// match pattern to sequence 

						/* check for all instances of motif in sequence */
						while(matcher.find()) {
							int pos = matcher.start();

							seq.setFoundPosition(pos);
							positionsToRemove.add(pos);
						}	
					}
					seq.removePositionToSearch(positionsToRemove);
				}
			}

			/* update list of positions to search & check if all sequence position have been found */
 			boolean keepSearching = false;
			for(SearchSequence seq : sequences) {

				if(!seq.getPositionsToSearch().isEmpty()) { // if not empty keep searching

					seq.updatePositionsToSearch();
					keepSearching = true;
				}
			}

			/* if all sequence positions have been found - break out of loop */
			if(!keepSearching) {
				break;
			}
		}

		/* assess genome positions and print to file */
		System.out.println("printing results");
		printGenomePositions(sequences, outputFile);
	}

	/**
	 * Load proteins to test from subset protein file
	 * 
	 * @param proteinsToTestFile	String - file containing list of proteins who's sequences we are testing
	 * @return proteins				HashSet<String> - list of proteins
	 */
	private static HashSet<String> loadProteinsToTest(String proteinsToTestFile){

		HashSet<String> proteins = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinsToTestFile))));

			String line = in.readLine();
			while(line!=null) {

				proteins.add(line.split("\t")[0]);
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteins;
	}

	/**
	 * Load RefSeq IDs associated to motif
	 * @param proteinInfo
	 * @return proteinMap	HashMap<String, HashSet<String>> - protein = id1|2|3..|4
	 */
	public static HashSet<String> loadProteinRefSeqIdMap(String proteinInfo, HashSet<String> proteins){

		HashSet<String> refSeqIds = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinInfo))));

			String line = in.readLine();
			while(line!=null) {

				String[] col = line.split("\t"); // [0] = protein, [1] = refSeqIds

				if(col.length > 1) {
					if(proteins.contains(col[0])) {
						refSeqIds.addAll(Arrays.asList(col[1].split("\\|")));
					}
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return refSeqIds;
	}

	private static List<SearchSequence> loadSequences(HashSet<String> refSeqIds, String fastaFile){

		List<SearchSequence> sequences = new ArrayList<>(); 
		int seqCount=0;
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));			

			String line = in.readLine();

			boolean readSeq = false; 
			String seq = "";
			String header = "";

			while(line!=null) {

				if(readSeq) {
					seq += line;
				}

				/* check if sequence is associated to a protein in network */
				if(line.startsWith(">")) {
					readSeq = false;

					/* if sequence has been loaded - check for motif in the sequence  */
					if(!seq.isEmpty()) {	
						seqCount++;
						System.out.print(seqCount + ".");
						if(seqCount%50 == 0 ) {
							System.out.println();
						}
						
						sequences.add(new SearchSequence(seq, header));
						seq = "";
					}

					//  >hg38_ncbiRefSeq_NM_001276352.2 range=chr1:67092165-67093579 5'pad=0 3'pad=0 strand=- repeatMasking=none					
					String id = line.split("[\\_\\s++\\.]")[2] + "_"+ line.split("[\\_\\s++\\.]")[3];

					if(refSeqIds.contains(id)) { // do not consider alternate chromosomes

						if(!line.contains("alt") && !line.contains("_fix")) {
							readSeq = true;
							header = line;
							seq = "";
						}
					}
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return sequences;
	}

	private static List<Integer> setOrderOfFilesToTest(int numberOfFiles){

		List<Integer> fileOrder = new ArrayList<>();
		for(int i=0; i<numberOfFiles; i++) {
			fileOrder.add(i);
		}
		Collections.shuffle(fileOrder);

		return fileOrder;
	}


	public static List<String> searchSeqForMotif(String motif, String sequence, String header) {

		List<String> positions = new ArrayList<>();

		Pattern pattern = Pattern.compile(motif);		// compile motif as REGEX
		Matcher matcher = pattern.matcher(sequence);	// match pattern to sequence 

		String chromosome = header.split(" ")[1].split("=")[1].split(":")[0];
		chromosome = chromosome+":";
		int startPos =  Integer.parseInt(header.split(" ")[1].split("=")[1].split(":")[1].split("-")[0]);
		int endPos =  Integer.parseInt(header.split(" ")[1].split("=")[1].split(":")[1].split("-")[1]);

		/* check for all instances of motif in sequence */
		while(matcher.find()) {

			/* obtain start position of motif */
			int pos = matcher.start();

			String genomicPosition = "";

			if(!header.contains("+")){

				int motifEnd = endPos-pos;
				int motifStart = motifEnd - 7;

				genomicPosition = chromosome+motifStart+"-"+motifEnd;
			}
			else{
				int motifStart = startPos+pos;
				int motifEnd = motifStart+7;
				genomicPosition = chromosome+motifStart+"-"+motifEnd;
			}

			positions.add(genomicPosition); 
		}
		return positions;
	}

	/**
	 * Set character map - hard coded with alphabet used throughout analysis 
	 * 
	 * Values follow known IUPAC regular expression e.g. R = "[AG]" 
	 * @return characterMap		HashMap<Character, String> 
	 */
	private static HashMap<Character, String> setCharacterMapForRegularExpression(){

		HashMap<Character, String> characterMap = new HashMap<>();

		characterMap.put('A', "A");
		characterMap.put('C', "C");
		characterMap.put('G', "G");
		characterMap.put('T', "T");
		characterMap.put('R', "[AG]");
		characterMap.put('Y', "[CT]");
		characterMap.put('D', "[ATG]");
		characterMap.put('B', "[TGC]");
		characterMap.put('H', "[AUC]");
		characterMap.put('V', "[AGC]");
		characterMap.put('*', ".");
		return characterMap;
	}

	/**
	 * Take input motif and translate to corresponding regular expression. 
	 * 
	 * @param motif				String - initial motif
	 * @param characterMap		HashMap<Character, String> - conversion map with regular expressions
	 * @return formattedMotif	String - motif formatted with regular expression
	 */
	private static String formatMotifWithRegularExpression(String motif, HashMap<Character, String> characterMap) {

		String formattedMotif = "";

		for(int i = 0; i < motif.length(); i++){
			formattedMotif += characterMap.get(motif.charAt(i));
		}

		return formattedMotif;
	}

	public static void printGenomePositions(List<SearchSequence> sequences, String outputFile) {


		/* determine genome positions */
		HashMap<String, HashSet<Integer>> positions = new HashMap<>(); // chromoSome = pos1,2,3,..,n
		System.out.println("combining positions");
		for(SearchSequence seq : sequences) {
			
			String chr = seq.getChromosome();
			HashSet<Integer> currentIdxs = seq.getFoundPositions();
			HashSet<Integer> currentPositions = new HashSet<>();

			/* determine genome position of index */
			for(int idx: currentIdxs) {
				if(seq.isPositive()) {
					currentPositions.add(seq.getStartPosition() + idx);
				} else {
					currentPositions.add(seq.getEndPosition() - idx);
				}
			}

			System.out.println(chr + "\t" + currentPositions.size());
			
			/* update positions per CHROMOSOME */
			if(positions.containsKey(chr)) {
				HashSet<Integer> allPositions = positions.get(chr);
				allPositions.addAll(currentPositions);

				positions.put(chr, allPositions);
			} else {
				positions.put(chr, currentPositions);
			}
		}

		/* print and sort positions; 1 line per chromosomes*/
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(outputFile)));
			System.out.println("printing");
			for(Entry<String, HashSet<Integer>> e: positions.entrySet()) {

				out.write(e.getKey() + ":\t");
				
				List<Integer> order = new ArrayList<Integer>(e.getValue());
				Collections.sort(order);
				
				System.out.println(e.getKey() + "\t" + order.size());
				
				int startIndex = order.get(0); // start index 
				int idx1 = startIndex;
				int idx2 = 0; 
				for(int i=1; i<order.size(); i++) {

					idx2 = order.get(i);

					if(idx2 != (idx1+1)) {
						out.write(startIndex + "-" + idx1 + ",");
						out.flush();

						startIndex = idx2;
					}
					
					idx1 = idx2;
				}
				
				
				
				out.write(startIndex + "-" + idx1 + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
