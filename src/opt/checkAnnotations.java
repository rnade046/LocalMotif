package opt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

public class checkAnnotations {

	public static void main(String[] args) {

		String inputFiles = "mapDegenMotifsToProteins_";
		String regMotifFiles = "";
		String outputFile = "";

		int numFiles = 992;
		HashMap<Integer,Integer> mapOfProteinFreq = new HashMap<>();
		
		calculateProteinFreqFromDegenMotifs(mapOfProteinFreq, inputFiles, numFiles);
		calculateProteinFreqFromRegularMotifs(mapOfProteinFreq, regMotifFiles);
		
		outputFreqOfProteinMap(mapOfProteinFreq, outputFile);

	}

	/**
	 * Iterate through all annotation files (degen motifs) to calculate the occurrence of 
	 * a motif annotation X number of proteins 
	 * 
	 * @param mapOfProteinFreq	HashMap<Integer, Integer> - to contain number of annotated proteins = occurrence 
	 * @param inputFiles		String - file path with annotation file prefix
	 * @param numFiles			int - number of degenMotifs files (corresponds to the file labeling)
	 */
	private static void calculateProteinFreqFromDegenMotifs(HashMap<Integer, Integer> mapOfProteinFreq, String inputFiles, int numFiles) {

		/* iterate over all degen motif annotation files */
		for(int i=0; i < numFiles; i++) {

			String annotationFile = inputFiles + i; // current annotation file (ex. mapDegenMotifsToProteins_200) 

			InputStream in;
			try {
				in = new FileInputStream(new File(annotationFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // no header

				while(line!=null) {

					Integer numOfProt = Integer.parseInt(line.split("\t")[2]); // number of proteins that contain this motif

					/* Update map */
					if(mapOfProteinFreq.containsKey(numOfProt)) {	
						mapOfProteinFreq.put(numOfProt, mapOfProteinFreq.get(numOfProt) + 1);
					} else {
						mapOfProteinFreq.put(numOfProt, 1);
					}

					line = input.readLine();
				}

				input.close();
			} catch(IOException e) {
				e.printStackTrace();
			}
		}
	}

	/**
	 * Iterate through all annotation files (degen motifs) to calculate the occurrence of 
	 * a motif annotation X number of proteins 
	 * 
	 * @param mapOfProteinFreq	HashMap<Integer, Integer> - to contain number of annotated proteins = occurrence 
	 * @param inputFiles		String - file path with annotation file prefix
	 * @param numFiles			int - number of degenMotifs files (corresponds to the file labeling)
	 */
	private static void calculateProteinFreqFromRegularMotifs(HashMap<Integer, Integer> mapOfProteinFreq, String inputFile) {

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // no header

			while(line!=null) {

				Integer numOfProt = Integer.parseInt(line.split("\t")[2]); // number of proteins that contain this motif

				/* Update map */
				if(mapOfProteinFreq.containsKey(numOfProt)) {	
					mapOfProteinFreq.put(numOfProt, mapOfProteinFreq.get(numOfProt) + 1);
				} else {
					mapOfProteinFreq.put(numOfProt, 1);
				}

				line = input.readLine();
			}

			input.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
	/**
	 * Print the number of annotated proteins to occurrence map
	 * 
	 * @param mapOfProteinFreq	HashMap<Integer, Integer> number of proteins annotated = occurrence in annotation file
	 * @param outputFile		String - file path
	 */
	private static void outputFreqOfProteinMap(HashMap<Integer, Integer> mapOfProteinFreq, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int numProt: mapOfProteinFreq.keySet()) {
				out.write(numProt + "\t" + mapOfProteinFreq.get(numProt) + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
