import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MotifFamily {

	public static void setMotifFamilyRepresentative(String motifFamilyFilePrefix, String significantMotifFile, File motifsDir, String motifsInfoFile, String height) {

		int familyCount = getMotifFamilyFiles(motifsDir, height);
 		for(int i=1; i<=familyCount; i++) {

			String motifFamilyFile = motifFamilyFilePrefix + i + ".tsv";

			/* Load motifs in motif family identified during hierarchical clustering step */ 
			HashSet<String> motifSet = loadMotifsInFamily(motifFamilyFile);

			/* Obtain significant clustering information for every motif */
			HashMap<String, Double> motifSignificantMap = getMotifSignificance(motifSet, significantMotifFile);

			/* Identify representative motif from list */
			String representativeMotif = getRepresentativeMotif(motifSignificantMap);

			printMotifInfo(motifsInfoFile, i, representativeMotif);

		}
	}

	public static void assessMotifFamilies(String motifsInfoFile, int familyNumber, String proteinToRefSeqIDsFile, 
			String motifInstancesPrefix, String outputFilePrefix, String annotatedProteinsFile, String fastaFile) {

		/* Load protein : refseqId list */
		HashMap<String, HashSet<String>> proteinToRefSeqIDsMap = loadProteinsToRefSeqIdsMap(proteinToRefSeqIDsFile);

		/* map of regular expression characters */
		HashMap<Character, String> characterMap = setCharacterMapForRegularExpression();

		System.out.println("Motif family : " + familyNumber);
		String motif = getMotifCorrespondingToFamilyNumber(motifsInfoFile, familyNumber);

		/* regular expression motif */
		String regexMotif = formatMotifWithRegularExpression(motif, characterMap);

		/* Load list of proteins annotated by representative motif */
		HashSet<String> annotatedProteinSet = loadProteinsAnnotatedByRepresentativeMotif(motif, annotatedProteinsFile);

		ArrayList<String> motifInstances = obtainMotifInstancesForRegexMotif(fastaFile, regexMotif, proteinToRefSeqIDsMap, annotatedProteinSet);
		double[][] ppm = calculatePPM(motifInstances, motifInstances.get(0).length());

		String outputFile = outputFilePrefix + familyNumber + ".tsv";
		String motifInstancesFile = motifInstancesPrefix + familyNumber;
		printMotifs(motifInstances, motifInstancesFile);
		printPPM(ppm, outputFile);
	}


	/**
	 * Load the list of proteins in the network and their contributing refSeq Ids. 
	 * @param proteinInfoFile	String - proteinName \t ID1|ID2|..|IDn
	 * 
	 * @return loadProteinsToRefSeqIdsMap	HashMap<String, HashSet<String>> - map {protein: Set<IDs>}
	 */
	private static HashMap<String, HashSet<String>> loadProteinsToRefSeqIdsMap(String proteinInfoFile){

		HashMap<String, HashSet<String>> proteinToRefSeqIdsMap = new HashMap<>();

		try {
			InputStream in = new FileInputStream(new File(proteinInfoFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				String[] col = line.split("\t");
				if(col.length == 2) {
					proteinToRefSeqIdsMap.put(col[0], new HashSet<String>(Arrays.asList(col[1].split("\\|")))); // col[0] = protein, col[1] = ID1|ID2|ID3
				}
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinToRefSeqIdsMap;
	}

	/**
	 * Load motifs in motif family. Store in set. 
	 * @param motifFamilyFile
	 * @return
	 */
	private static HashSet<String> loadMotifsInFamily(String motifFamilyFile){
		HashSet<String> motifSet = new HashSet<>();

		try {
			InputStream in = new FileInputStream(new File(motifFamilyFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {

				motifSet.add(line);
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifSet;
	}

	private static HashMap<String, Double> getMotifSignificance(HashSet<String> motifSet, String significantMotifsFile){

		HashMap<String, Double> motifSignificanceMap = new HashMap<>();

		/* Get p-value for all motifs in this family 
		 * rank in ascending order; take smallest one; 
		 * if tie take motif with least amount of degen characters */ 
		try {
			InputStream in = new FileInputStream(new File(significantMotifsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {

				String motif = line.split("\t")[0];

				if(motifSet.contains(motif)) {
					motifSignificanceMap.put(motif, Double.parseDouble(line.split("\t")[3])); // [3] = p-value
				}

				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


		return motifSignificanceMap;
	}

	private static String getRepresentativeMotif(HashMap<String, Double> motifSignificanceMap) {

		String repMotif = "";

		List<Entry<String, Double>> list = new LinkedList<Entry<String, Double>>(motifSignificanceMap.entrySet());  
		list.sort(Entry.comparingByValue());

		double minPval = list.get(0).getValue();
		boolean repeatPval = true; 
		int countIdx = 1;

		/* Check if more than one motif has the same p-value */
		if(list.size() > 1) {
			while(repeatPval) {
				if(list.get(countIdx).getValue() == minPval) {
					countIdx++;
				} else {
					repeatPval = false;
				}
			}
		}
		/* Compare motifs with identical p-value*/
		if(countIdx > 1) {

			/* Identify motifs with the same low p-value */
			String[] motifs = new String[countIdx];
			for(int i = 0; i < countIdx; i++) {
				motifs[i] = list.get(i).getKey();
			}

			/* Check for the least degen motif */
			HashSet<Character> degenCharacterSet = new HashSet<>(Arrays.asList('R', 'D', 'H', 'V', 'B', '*'));
			int maxDegenCount = 8;
			for(String motif: motifs) {

				int degenCount = 0;
				for(int k=0; k < motif.length() ; k++) {

					/* keep count if added character is a degenerate character */
					if(degenCharacterSet.contains(motif.charAt(k))) {
						degenCount++;
					}
				}

				if(degenCount < maxDegenCount) {
					repMotif = motif;
					maxDegenCount = degenCount;
				}
			}

		} else { 
			repMotif = list.get(0).getKey();
		}

		return repMotif;
	}

	private static HashSet<String> loadProteinsAnnotatedByRepresentativeMotif(String representativeMotif, String annotatedProteinsFile){

		HashSet<String> proteinSet = new HashSet<>();

		try {
			InputStream in = new FileInputStream(new File(annotatedProteinsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				String motif = line.split("\t")[0];

				/* if line corresponds to representative motif, store all it's annotated proteins to set */
				if(motif.equals(representativeMotif)) {

					String[] annotatedProteins = line.split("\t")[2].split("\\|");
					proteinSet.addAll(Arrays.asList(annotatedProteins));
				}
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinSet;
	}

	private static ArrayList<String> obtainMotifInstancesForRegexMotif(String fastaFile, String regexMotif, HashMap<String, HashSet<String>> proteinToRefSeqIdsMap, HashSet<String> annotatedProteinSet) {

		/* Search fasta sequence; find instance of motif (ie. convert degen motif to non degen) ; print PWM to file */
		ArrayList<String> instanceOfMotifList = new ArrayList<>();

		/* Search for motif in a given protein; only one instance of a motif will be kept if found 
		 * multiple times in sequence or in different refseqIds (variants) > this is accomplished by 
		 * making a hash set for every protein (therefore only one instance of any possible motif will
		 * be kept per protein (from all 3'UTR variants) */ 

		int count = 0;
		for(String protein : annotatedProteinSet) {
			count++;

			if(count%10 == 0) {
				System.out.print(count+".");
			}

			if(count%100 == 0) {
				System.out.println();
			}

			/* find instances of motifs in fasta sequences */
			HashSet<String> currentInstancesOfMotif = searchSequencesForRegexMotifForOneProtein(fastaFile, regexMotif, proteinToRefSeqIdsMap.get(protein));

			/* merge instances with those found */
			instanceOfMotifList.addAll(currentInstancesOfMotif);
		}
		return instanceOfMotifList;
	}

	private static double[][] calculatePPM(ArrayList<String> motifInstances, int motifLength){

		double[][] pfm = new double[4][motifLength];
		double[][] ppm = new double[4][motifLength];

		/* compute position frequency matrix */
		for(String motif: motifInstances) {
			for(int i=0; i<motifLength; i++) {
				switch(motif.charAt(i)) {
				case 'A': 
					pfm[0][i] += 1;
					break;
				case 'C': 
					pfm[1][i] += 1;
					break;
				case 'G':
					pfm[2][i] += 1;
					break;
				case 'T':
					pfm[3][i] += 1;
					break;
				}
			}
		}

		/* convert position frequency matrix to position probability matrix */
		for(int i=0; i<4; i++) {
			for(int j=0; j<motifLength; j++) {
				ppm[i][j] = pfm[i][j]/ (double) motifInstances.size();
			}
		}

		return ppm;
	}

	private static void printMotifs(ArrayList<String> motifInstances, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(String motif: motifInstances) {
				out.write(motif + "\n");
				out.flush();
			}


			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void printPPM(double[][] ppm, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=0; i<ppm.length; i++) {
				for(int j=0;j<ppm[i].length; j++) {
					out.write(ppm[i][j] + "\t");
				}
				out.write("\n");
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void printMotifInfo(String outputFile, int family, String repMotif) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile), true));

			if(new File(outputFile).length() == 0) {
				out.write("RepMotif\tFamilyNumber\n");
			}

			out.write(repMotif+ "\t" + family + "\n");
			out.flush();
			out.close();

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/**
	 * Take input motif and translate to corresponding regular expression.
	 * 
	 * @param motif        String - initial motif
	 * @param characterMap HashMap<Character, String> - conversion map with regular
	 *                     expressions
	 * @return formattedMotif String - motif formatted with regular expression
	 */
	private static String formatMotifWithRegularExpression(String motif, HashMap<Character, String> characterMap) {

		String formattedMotif = "";

		for (int i = 0; i < motif.length(); i++) {
			formattedMotif += characterMap.get(motif.charAt(i));
		}

		return formattedMotif;
	}

	/**
	 * Search for regex motif in the provided sequence
	 * 
	 * @param motif    String - motif formatted with regular expression
	 * @param sequence String - RNA sequence to search
	 * @return motifFound Boolean - true if motif is found in sequence
	 */
	private static HashSet<String> searchSeqForMotif(String motif, String sequence) {

		HashSet<String> motifInstances = new HashSet<>();

		Pattern pattern = Pattern.compile(motif); // compile motif as REGEX
		Matcher matcher = pattern.matcher(sequence); // match pattern to sequence

		/* check if motif is contained in the sequence */
		while (matcher.find()) {
			// Get the matched substring
			motifInstances.add(matcher.group());
		}
		return motifInstances;
	}

	/**
	 * Set character map - hard coded with alphabet used throughout analysis
	 * 
	 * Values follow known IUPAC regular expression e.g. R = "[AG]"
	 * 
	 * @return characterMap HashMap<Character, String>
	 */
	private static HashMap<Character, String> setCharacterMapForRegularExpression() {

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

	private static HashSet<String> searchSequencesForRegexMotifForOneProtein(String fastaFile, String motif, HashSet<String> refSeqIds) {

		HashSet<String> instances = new HashSet<>();
		try {
			BufferedReader in = new BufferedReader(
					new InputStreamReader(new FileInputStream(new File(fastaFile))));

			String line = in.readLine();

			boolean readSeq = false;
			String seq = "";
			String id = "";

			while (line != null) {

				if (readSeq) {
					seq += line;
				}
				/* check if sequence is associated to a protein in network */
				if (line.startsWith(">")) {
					readSeq = false;

					/* if sequence has been loaded - check for motif in the sequence */
					if (!seq.isEmpty()) {

						instances.addAll(searchSeqForMotif(motif, seq));
					}
					// >hg38_ncbiRefSeq_NM_001276352.2 range=chr1:67092165-67093579 5'pad=0 3'pad=0
					// strand=- repeatMasking=none
					id = line.split("[\\_\\s++\\.]")[2] + "_" + line.split("[\\_\\s++\\.]")[3];

					/* reinitialize parameters for next sequence */
					if (refSeqIds.contains(id)) {
						// chromosomes
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
		return instances;
	}


	private static String getMotifCorrespondingToFamilyNumber(String motifFamilyFile, int familyNumber){

		String motif = "";
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(motifFamilyFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {

				String[] col = line.split("\t"); // [0] = IUPAC motif, [1] = family number
				if(Integer.parseInt(col[1]) == familyNumber) {
					motif = col[0];
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motif;
	}

	private static int getMotifFamilyFiles(File dir, String height) {
		
		File[] files = null;
		int fileCount = 0;
		if (dir.isDirectory()) {

			FilenameFilter filter = new FilenameFilter() {
				public boolean accept(File f, String name) 
				{ 
					return name.startsWith("MotifFamily_h" + height); 
				} 
			};
			
			files = dir.listFiles(filter);
			fileCount = files.length;
		}

		return fileCount;
	}

}
