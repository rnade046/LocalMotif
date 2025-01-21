package sampling;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class ProteinAnnotations {

	private int lowerBoundToSample;
	private int upperBoundToSample;

	private HashSet<String> proteinSet;

	public ProteinAnnotations(int _lowerBoundToSample, int _upperBoundToSample, HashSet<String> _proteinSet) {
		this.lowerBoundToSample = _lowerBoundToSample;
		this.upperBoundToSample = _upperBoundToSample;

		this.proteinSet = _proteinSet;
	}

	public void calculateProteinAnnotationFrequency(String degenMotifAnnotationPrefix, String annotationCompanionFileString,
			int index1, int index2, String outputProteinFrequencyFilePrefix) {

		/***
		 * Given an annotation list of type 
		 * {annotation; protein id list; protein symbol list}
		 * 
		 * We want to compute the frequency of every protein in the list, and
		 * we want to output these frequencies in a text file
		 * to be used for a weighted Monte Carlo Sampling approach.
		 */

		HashMap<String, Integer> proteinToFrequencyMap = new HashMap<>();

//		System.out.print("Searching degen motif annotation files:");
		//for(int i=0; i<numDegenFiles; i++) {
		for(int i=index1; i<= index2; i++) {	
//			if(i%100==0) {
//				System.out.println();
//			}
//			if(i%10==0) {
//				System.out.print(i +".");
//			}
			
			System.out.println("Searching degen motif annotation files: " + i );
			
			String degenAnnotationFile = degenMotifAnnotationPrefix + i + ".tsv";
			String annotationCompanionFile = annotationCompanionFileString + i;

			HashSet<String> motifsToTest = loadMotifsToTest(annotationCompanionFile);
			computeFrequencyOfProteins(degenAnnotationFile, motifsToTest, proteinToFrequencyMap);


			System.out.println("Printing protein annotation frequencies");
			printProteinsToOccurence(outputProteinFrequencyFilePrefix + i , proteinToFrequencyMap);
		}
	}

	/**
	 * Load annotation list and store in HashMap proteins and their number of occurrence in annotation list.
	 * @param inputFile					text file containing annotation list {annotation; protein id list; protein symbol list}
	 * @return proteinToFrequencyMap	HashMap<String, Integer> {protein id: number of occurrence}
	 */
	private void computeFrequencyOfProteins(String inputFile, HashSet<String> motifsToTest, HashMap<String, Integer> proteinToFrequencyMap){

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			//int count = 0;
			String line = in.readLine(); // no header

			while(line!=null) {

				if(motifsToTest.contains(line.split("\\t")[0])) {

					String[] protein_ids = line.split("\\t")[2].split("\\|"); // idx[2] = protein (name) list 
					ArrayList<String> proteinList = checkProteinsInNetwork(protein_ids);

					if(proteinList.size() >= lowerBoundToSample && proteinList.size() <= upperBoundToSample) {
						/* For all proteins of a given annotation; if in list update number of occurrence, otherwise initialize */
						for(int i=0; i<proteinList.size(); i++) {
							if(proteinToFrequencyMap.containsKey(proteinList.get(i))) {
								proteinToFrequencyMap.replace(proteinList.get(i), proteinToFrequencyMap.get(proteinList.get(i)) + 1);
							} else {
								proteinToFrequencyMap.put(proteinList.get(i), 1);
							}
						}
					}
				}
				line = in.readLine();
				//count++;
			}			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Output annotated proteins and their occurrences in a text file
	 * @param outputFile				text file to contain Protein : Occurrence
	 * @param proteinToOccurrenceMap	map of {protein: occurrence}
	 */
	private static void printProteinsToOccurence(String outputFile, HashMap<String, Integer> proteinToOccurrenceMap) {
		try {
			OutputStreamWriter out = new OutputStreamWriter(new FileOutputStream(new File(outputFile)));
			/* iterate through map */
			for(String protein: proteinToOccurrenceMap.keySet()) {
				out.write(protein + "\t" + proteinToOccurrenceMap.get(protein) + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private ArrayList<String> checkProteinsInNetwork(String[] proteinList) {

		ArrayList<String> finalProteinList = new ArrayList<>();

		for(String prot: proteinList) {
			if(this.proteinSet.contains(prot)) {
				finalProteinList.add(prot);
			}
		}

		return finalProteinList;
	}

	private HashSet<String> loadMotifsToTest(String annotationCompanionFile){

		HashSet<String> motifSet = new HashSet<>();

		try {
			InputStream in = new FileInputStream(new File(annotationCompanionFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				motifSet.add(line.split("\t")[0]);
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifSet;
	}

	public void combineProteinFrequencyData(String protFreqFilePrefix, File dir, String outputFile) {
		
		HashMap<String, Integer> proteinToOccurrenceMap = new HashMap<>();
		
		int numFiles = dir.list().length;
		for(int i=0; i<numFiles; i++) {
			
			/* Load individual file protein frequencies and update proteinToOccurenceMap */
			String protFreqFile = protFreqFilePrefix + i;
			try {
				InputStream in = new FileInputStream(new File(protFreqFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine();

				while(line != null) {
					String[] col = line.split("\t");
					
					if(proteinToOccurrenceMap.containsKey(col[0])) {
						proteinToOccurrenceMap.put(col[0], proteinToOccurrenceMap.get(col[0]) + Integer.parseInt(col[1]));
					} else {
						proteinToOccurrenceMap.put(col[0], Integer.parseInt(col[1]));
					}
					line = input.readLine();
				}

				input.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			printProteinsToOccurence(outputFile, proteinToOccurrenceMap);
		}
	}
	
}