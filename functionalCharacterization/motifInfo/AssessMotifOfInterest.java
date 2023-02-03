package motifInfo;

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
import java.util.Map.Entry;

import positionConservation.SequenceUtils;

public class AssessMotifOfInterest {

	private static String proteinAtlasFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/subcellular_location.tsv";
	private static String fastaFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/input_files/human_3UTRsequences.txt";
	private static String proteinInfoFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/corrNetTop2_proteinsInNetwork_info.tsv";
	private static String annotatedProteinsFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/corrNetTop2_protFreqAnnotation.tsv";

	public static void main(String[] args) {

		String motifsOfInterest = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/motifTest.txt";
		String annotationSubset = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/corrNetTop2-400_coreTPD_p0.4_coreProteinsByMotif.tsv";

		String outputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/coreTPD_motifOfInterest_bins2_";
		String localizationFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/coreTPD_motifOfInterest_localization2_";


		//getInformationForMotifsOfInterestWithRegex(motifsOfInterest, annotationSubset, outputFile, localizationFile);
		getInformationForMotifsOfInterestWithBins(motifsOfInterest, annotationSubset, outputFile, localizationFile);

	}
	public static void getInformationForMotifsOfInterestWithBins(String motifsOfInterestFile, String annotationSubsetFile, String outputPrefixFile, String localizationFile) {

		/* Get motifs of interest  */
		HashSet<String> motifsToTest = loadMotifsToTest(motifsOfInterestFile);
		System.out.println("Motifs to test : " + motifsToTest.size());
		
		/* Get proteins in network */
		HashSet<String> annotatedProteins = loadAnnotatedProteins(annotatedProteinsFile);
		
		/* For motif of interest, obtain list of proteins associated */
		for(String motif: motifsToTest) {
			System.out.println("testing motif : " + motif);

			/* For current motif, get list of proteins*/
			HashSet<String> proteinSet = getProteinsAssociatedToMotif(motif, annotationSubsetFile, annotatedProteins);
			System.out.println("load proteins: " + proteinSet.size());	

			/* For list of proteins get their associated refSeqIds*/
			HashMap<String, HashSet<String>> refSeqdIdsMap = getListOfRefSeqIds(proteinSet);
			System.out.println("load refSeqIds");

			/* Get protein localizations */
			HashMap<String, String> proteinLocalizationMap = determineProteinLocalization(proteinSet);
			System.out.println("got localizations");

			//HashSet<String> motifPossibilities = SequenceUtils.getPossibleMotifInstances(e.getKey());

			HashSet<String> motifPossibilities = SequenceUtils.getPossibleMotifInstances(motif);

			/* Search for motif positions */
			HashMap<String, List<String>> motifPositionByIdMap = new HashMap<>();
			int idCount = 1;
			for(HashSet<String> refSeqIds: refSeqdIdsMap.values()) { // e.key = protein, e.value = refSeqIDs

				System.out.print(idCount+ ".");
				if(idCount%50 ==0) {
					System.out.println();
				}

				motifPositionByIdMap.putAll(determineMotifPositions(refSeqIds, motifPossibilities));
				idCount++; 
			}

			/* Print results */
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputPrefixFile + motif + ".tsv")));

				out.write("Motif = " + motif + " | numberOfProtein = " + proteinSet.size() + "\n");

				/*For each protein */
				int protCount = 1;
				System.out.println("proteins to asssess : " + proteinSet.size());
				for(String protein: proteinSet) {

					System.out.print(protCount + ".");
					if(protCount%50 ==0) {
						System.out.println();
					}

					/* output protein info */
					out.write(protein + "\t" + proteinLocalizationMap.get(protein) + "\t");

					HashSet<String> refSeqIdsAssociatedToProtein = refSeqdIdsMap.get(protein);
					for(String id: refSeqIdsAssociatedToProtein) {

						if(motifPositionByIdMap.containsKey(id)) {

							out.write(id + "=");

							List<String> motifPositions = motifPositionByIdMap.get(id);
							for(String pos : motifPositions) {
								out.write(pos + ",");
							}
							out.write("|");
						}
					}
					out.write("\n");
					out.flush();
					protCount++;
				}

				out.close();
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			HashMap<String, Double> localizationCount = new HashMap<>();

			for(String local : proteinLocalizationMap.values()) {

				String[] localizations = local.split("\\;");
				double increment = 1/localizations.length;

				for(String l : localizations) {
					if(localizationCount.containsKey(l)) {
						localizationCount.put(l, localizationCount.get(l)+ increment);
					} else {
						localizationCount.put(l, increment);
					}
				}	
			}


			printLocalization(localizationFile + motif + ".tsv", localizationCount);
		}

	}

	public static void getInformationForMotifsOfInterestWithRegex(String motifsOfInterestFile, String annotationSubsetFile, String outputPrefixFile, String localizationFile) {

		/* Get motifs of interest  */
		HashSet<String> motifsToTest = loadMotifsToTest(motifsOfInterestFile);
		System.out.println("Motifs to test : " + motifsToTest.size());
		
		/* Get proteins in network */
		HashSet<String> annotatedProteins = loadAnnotatedProteins(annotatedProteinsFile);

		/* For motif of interest, obtain list of proteins associated */

		/* Iterate for each motif */
		for(String motif: motifsToTest) {
			System.out.println("testing motif : " + motif);

			/* For current motif, get list of proteins*/
			HashSet<String> proteinSet = getProteinsAssociatedToMotif(motif, annotationSubsetFile, annotatedProteins);
			System.out.println("load proteins: " + proteinSet.size());	

			/* For list of proteins get their associated refSeqIds*/
			HashMap<String, HashSet<String>> refSeqdIdsMap = getListOfRefSeqIds(proteinSet);
			System.out.println("load refSeqIds");

			/* Get protein localizations */
			HashMap<String, String> proteinLocalizationMap = determineProteinLocalization(proteinSet);
			System.out.println("got localizations");

			//HashSet<String> motifPossibilities = SequenceUtils.getPossibleMotifInstances(e.getKey());

			String regexMotif = RegexSearch.formatMotifWithRegularExpression(motif, RegexSearch.setCharacterMapForRegularExpression()); 

			/* Search for motif positions */
			HashMap<String, List<String>> motifPositionByIdMap = new HashMap<>();
			int idCount = 1;
			for(HashSet<String> refSeqIds: refSeqdIdsMap.values()) { // e.key = protein, e.value = refSeqIDs

				System.out.print(idCount+ ".");
				if(idCount%50 ==0) {
					System.out.println();
				}

				motifPositionByIdMap.putAll(RegexSearch.searchForMotifPostions(regexMotif, fastaFile, refSeqIds));
				idCount++;
			}

			/* Print results */
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputPrefixFile + motif + ".tsv")));

				out.write("Motif = " + motif + " | numberOfProtein = " + proteinSet.size() + "\n");

				/*For each protein */
				int protCount = 1;
				System.out.println("proteins to asssess : " + proteinSet.size());
				for(String protein: proteinSet) {

					System.out.print(protCount + ".");
					if(protCount%50 ==0) {
						System.out.println();
					}

					/* output protein info */
					out.write(protein + "\t" + proteinLocalizationMap.get(protein) + "\t");

					HashSet<String> refSeqIdsAssociatedToProtein = refSeqdIdsMap.get(protein);
					for(String id: refSeqIdsAssociatedToProtein) {

						if(motifPositionByIdMap.containsKey(id)) {

							out.write(id + "=");

							List<String> motifPositions = motifPositionByIdMap.get(id);
							for(String pos : motifPositions) {
								out.write(pos + ",");
							}
							out.write("|");
						}
					}
					out.write("\n");
					out.flush();
					protCount++;
				}

				out.close();
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			HashMap<String, Double> localizationCount = new HashMap<>();

			for(String local : proteinLocalizationMap.values()) {

				String[] localizations = local.split("\\;");
				double increment = 1/localizations.length;

				for(String l : localizations) {
					if(localizationCount.containsKey(l)) {
						localizationCount.put(l, localizationCount.get(l)+ increment);
					} else {
						localizationCount.put(l, increment);
					}
				}	
			}


			printLocalization(localizationFile + motif + ".tsv", localizationCount);
		}

	}

	private static HashSet<String> loadMotifsToTest(String motifsToTestFile){

		HashSet<String> motifsToTest = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(motifsToTestFile))));

			String line = in.readLine(); // no-header

			while(line!=null) {

				motifsToTest.add(line.split("\t")[1]);
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifsToTest;
	}
	
	private static HashSet<String> loadAnnotatedProteins(String protFile){
		
		HashSet<String> annotatedProteins = new HashSet<>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(protFile))));

			String line = in.readLine(); // no-header

			while(line!=null) {
				
				annotatedProteins.add(line.split("\t")[0]);
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return annotatedProteins;
	}

	private static HashSet<String> getProteinsAssociatedToMotif(String motif, String proteinsAssociatedToMotifsFile, HashSet<String> annotatedProteins){

		HashSet<String> proteinsAssociatedToMotif = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinsAssociatedToMotifsFile))));

			String line = in.readLine(); // no-header

			while(line!=null) {

				String currentMotif = line.split("\t")[0];

				if(motif.equals(currentMotif)) {
					
					for(String prot : line.split("\t")[2].split("\\|")) {
						if(annotatedProteins.contains(prot)) {
							proteinsAssociatedToMotif.add(prot);
						}
					}
					break;
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinsAssociatedToMotif;
	}

	private static HashMap<String, HashSet<String>> getListOfRefSeqIds(HashSet<String> proteinSet){
		HashMap<String, HashSet<String>> refSeqIdMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinInfoFile))));

			String line = in.readLine(); // no-header

			while(line!=null && refSeqIdMap.size() < proteinSet.size()) {

				String p = line.split("\t")[0]; // [0] protein
				if(proteinSet.contains(p)) {

					if(line.split("\t").length > 1) {
						HashSet<String> refSeqIds = new HashSet<>(Arrays.asList(line.split("\t")[1].split("\\|")));
						refSeqIdMap.put(p, refSeqIds);
					} else { 
						System.out.println(p);
					}
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return refSeqIdMap;
	}

	@SuppressWarnings("unused")
	private static HashMap<String, List<String>> determineMotifPositions(HashSet<String> refSeqIds, HashSet<String> motifPossibilities){

		HashMap<String, List<String>> motifPositionsPerId = new HashMap<>();
		int motifLength = 8;
		int bins = 100;

		for(String id : refSeqIds) {
			/* get sequence from FASTA file */
			String seq = SequenceUtils.getFastaForId(id, fastaFile);

			/* format sequence into bins */
			String[] binnedSeq = SequenceUtils.formatSequenceInBins(seq, bins);
			int[] motifPositions = new int[bins];

			/* determine if motif is in bin */
			for(int i=0; i<binnedSeq.length; i++) {

				String subsequence = binnedSeq[i];

				/* search for motif in subsequence */
				for(int j=0; j<subsequence.length()-motifLength; j++) {

					/* update count when motif is found */
					if(motifPossibilities.contains(subsequence.substring(j, j+8))) {
						motifPositions[i] += 1;
					}
				}
			}
			List<String> motifInBins = new ArrayList<>();
			/* Assess bins */
			for(int i=0; i<binnedSeq.length; i++) {
				if(motifPositions[i] >= 1) {
					motifInBins.add(String.valueOf(i+1));
				}
			}
			motifPositionsPerId.put(id, motifInBins);
		}
		return motifPositionsPerId;
	}

	private static HashMap<String, String> determineProteinLocalization(HashSet<String> proteinSet) {

		HashMap<String, String> proteinLocalizationMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinAtlasFile))));

			String line = in.readLine(); // no-header

			while(line!=null && proteinLocalizationMap.size() < proteinSet.size()) {

				String protein = line.split("\t")[1];

				if(proteinSet.contains(protein)) {
					proteinLocalizationMap.put(protein, line.split("\t")[3]); // [3] = localization
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinLocalizationMap;
	}

	private static void printLocalization(String file, HashMap<String, Double> localizationCount) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(file)));

			out.write("Localization\tCount\n");
			for(Entry<String, Double> e : localizationCount.entrySet()) {

				out.write(e.getKey() + "\t" + e.getValue() + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
