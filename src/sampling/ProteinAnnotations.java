package sampling;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashMap;

public class ProteinAnnotations {

	private int lowerBoundToSample;
	private int upperBoundToSample;

	public ProteinAnnotations(int _lowerBoundToSample, int _upperBoundToSample) {
		this.lowerBoundToSample = _lowerBoundToSample;
		this.upperBoundToSample = _upperBoundToSample;
	}

	public void calculateProteinAnnotationFrequency(String degenMotifAnnotationPrefix, int numDegenFiles, String outputProteinFrequencyFile) {

		/***
		 * Given an annotation list of type 
		 * {annotation; protein id list; protein symbol list}
		 * 
		 * We want to compute the frequency of every protein in the list, and
		 * we want to output these frequencies in a text file
		 * to be used for a weighted Monte Carlo Sampling approach.
		 */

		HashMap<String, Integer> proteinToFrequencyMap = new HashMap<>();
		
		System.out.print("Searching degen motif annotation files:");
		for(int i=0; i<numDegenFiles; i++) {
			if(i%100==0) {
				System.out.println();
			}
			if(i%10==0) {
				System.out.print(i +".");
			}
			
			String degenAnnotationFile = degenMotifAnnotationPrefix + i;
			computeFrequencyOfProteins(degenAnnotationFile, proteinToFrequencyMap);
		}
		
		System.out.println("Printing protein annotation frequencies");
		printProteinsToOccurence(outputProteinFrequencyFile, proteinToFrequencyMap);

	}

	/**
	 * Load annotation list and store in HashMap proteins and their number of occurrence in annotation list.
	 * @param inputFile					text file containing annotation list {annotation; protein id list; protein symbol list}
	 * @return proteinToFrequencyMap	HashMap<String, Integer> {protein id: number of occurrence}
	 */
	private void computeFrequencyOfProteins(String inputFile, HashMap<String, Integer> proteinToFrequencyMap){

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			//int count = 0;
			String line = in.readLine(); // no header

			while(line!=null) {

				int numOfAnnotatedProteins = Integer.parseInt(line.split("\\t")[1]);
				if(numOfAnnotatedProteins >= this.lowerBoundToSample && numOfAnnotatedProteins<= this.upperBoundToSample) {


					String[] protein_ids = line.split("\\t")[2].split("\\|"); // idx[2] = protein (name) list 

					/* For all proteins of a given annotation; if in list update number of occurrence, otherwise initialize */
					for(int i=0; i<protein_ids.length; i++) {
						if(proteinToFrequencyMap.containsKey(protein_ids[i])) {
							proteinToFrequencyMap.replace(protein_ids[i], proteinToFrequencyMap.get(protein_ids[i]) + 1);
						} else {
							proteinToFrequencyMap.put(protein_ids[i], 1);
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
}