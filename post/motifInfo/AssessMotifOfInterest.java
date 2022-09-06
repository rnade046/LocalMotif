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

	public static void main(String[] args) {
		
		String motifsOfInterest = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/coreTPD_motifsOfInterest.tsv";
		String annotationSubset = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/corrNetTop2-400_coreTPD_p0.4_coreProteinsByMotif.tsv";
		
		String outputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/coreTPD_motifOfInterest_";
		
		getInformationForMotifsOfInterest(motifsOfInterest, annotationSubset, outputFile);
		
	}

	public static void getInformationForMotifsOfInterest(String motifsOfInterestFile, String annotationSubsetFile, String outputPrefixFile) {

		/* Get motifs of interest  */
		HashSet<String> motifsToTest = loadMotifsToTest(motifsOfInterestFile);
		System.out.println("Motifs to test : " + motifsToTest.size());
		/* For motif of interest, obtain list of proteins associated */
		HashMap<String, HashSet<String>> proteinsToMotifsMap = getProteinsAssociatedToMotif(motifsToTest, annotationSubsetFile);
		
		/* Iterate for each motif */
		for(Entry<String, HashSet<String>> e : proteinsToMotifsMap.entrySet()) {
			System.out.println("testing motif : " + e.getKey());
			HashSet<String> motifPossibilities = SequenceUtils.getPossibleMotifInstances(e.getKey());

			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputPrefixFile + e.getKey() + ".tsv")));

				out.write("Motif = " + e.getKey() + " | numberOfProtein = " + e.getValue().size() + "\n");

				/*For each protein */
				int protCount = 1;
				System.out.println("proteins to asssess : " + e.getValue().size() );
				for(String protein: e.getValue()) {
					System.out.print(protCount + ".");
					if(protCount%50 ==0) {
						System.out.println();
					}

					/* get the list of RefSeqIDs */
					HashSet<String> refSeqIds = getListOfRefSeqIds(protein);

					/* For each RefSeqId; get FASTA sequence - determine bin position of motif */
					HashMap<String, List<String>> motifPositions = determineMotifPositions(refSeqIds, motifPossibilities);

					/* determine main localization of protein (from Human Protein Atlas) */	
					String localization = determineProteinLocalization(protein);

					/* output protein info */
					out.write(protein + "\t" + localization + "\t");

					for(Entry<String, List<String>> refSeqEntry : motifPositions.entrySet()) {
						
						if(!refSeqEntry.getValue().isEmpty()){
							
							out.write(refSeqEntry.getKey() + "=");

							for(String bin : refSeqEntry.getValue()) {
								out.write(bin + ",");
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

	private static HashMap<String, HashSet<String>> getProteinsAssociatedToMotif(HashSet<String> motifsToTest, String proteinsAssociatedToMotifsFile){

		HashMap<String, HashSet<String>> motifToProteinMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinsAssociatedToMotifsFile))));

			String line = in.readLine(); // no-header

			while(line!=null && motifToProteinMap.size() < motifsToTest.size()) {

				String motif = line.split("\t")[0];
				if(motifsToTest.contains(motif)) {

					HashSet<String> proteinSet = new HashSet<>(Arrays.asList(line.split("\t")[2].split("\\|")));
					motifToProteinMap.put(motif, proteinSet);
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifToProteinMap;
	}

	private static HashSet<String> getListOfRefSeqIds(String protein){

		HashSet<String> refSeqIds = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinInfoFile))));

			String line = in.readLine(); // no-header

			while(line!=null) {

				String p = line.split("\t")[0];
				if(p.equals(protein)) {
					refSeqIds.addAll(Arrays.asList(line.split("\t")[1].split("\\|")));
					break;
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return refSeqIds;
	}

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

	private static String determineProteinLocalization(String protein) {

		String localization = "";

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinAtlasFile))));

			String line = in.readLine(); // no-header

			while(line!=null) {

				if(line.split("\t")[1].equals(protein)) {
					localization = line.split("\t")[3];
					break;
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return localization;
	}
}
