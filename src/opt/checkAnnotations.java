package opt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Properties;

public class checkAnnotations {

	public static void main(String[] args) throws FileNotFoundException, IOException {
		
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));
		
		String wd = params.getProperty("working_directory");
		String projectName = params.getProperty("project_name");
		
		String degenMotifFilesPrefix = wd + "annotationFiles/"+ projectName+  "_degenMotifMappedToProteinsInNetwork_";
		String regMotifFiles =  wd + "annotationFiles/" + params.getProperty("motifAnnotationFile");
		
		String protFreqFile = wd + projectName + "_checkAnnotations_numberOfProteinsAnnotatedOccurence.tsv";
		String unseenNumberOfAnnotatedProteinsFile = wd + projectName + "_checkAnnotations_numberOfProteinsNotToSample.tsv";
		String seenNumberOfAnnotatedProteinsFile = wd + projectName + "_checkAnnotations_numberOfProteinsToSample.tsv";
		
		int numFiles = Integer.parseInt(params.getProperty("numDegenMotifFiles"));
		HashMap<Integer,Integer> mapOfProteinFreq = new HashMap<>();
		
		System.out.println("**Checking degen motifs annotation files**");
		calculateProteinFreqFromDegenMotifs(mapOfProteinFreq, degenMotifFilesPrefix, numFiles);
		
		System.out.println("**Checking regular motifs annotation file**\n");
		calculateProteinFreqFromRegularMotifs(mapOfProteinFreq, regMotifFiles);

		System.out.println("**Outputing freq map**");
		outputFreqOfProteinMap(mapOfProteinFreq, protFreqFile);
		outputNumberOfProteinsNotSeenInAnnotationMap(mapOfProteinFreq, unseenNumberOfAnnotatedProteinsFile);
		outputNumberOfProteinsSeenInAnnotationMap(mapOfProteinFreq, seenNumberOfAnnotatedProteinsFile);
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
			
			if(i%10==0) {
				System.out.print(i + ".");
			}
			
			if(i%100==0) {
				System.out.println();
			}
			
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
		System.out.println("Done\n");
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
			out.write("NumAnnotatedProteins\tOccurrence\n");
			
			for(int numProt: mapOfProteinFreq.keySet()) {
				out.write(numProt + "\t" + mapOfProteinFreq.get(numProt) + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void outputNumberOfProteinsNotSeenInAnnotationMap(HashMap<Integer,Integer> mapOfProteinFreq, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=3; i<=2000; i++) {

				if(!mapOfProteinFreq.containsKey(i)) {
					out.write(i + "\n");
					out.flush();
				}
				
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void outputNumberOfProteinsSeenInAnnotationMap(HashMap<Integer,Integer> mapOfProteinFreq, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=3; i<=2000; i++) {

				if(mapOfProteinFreq.containsKey(i)) {
					out.write(i + "\n");
					out.flush();
				}
				
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
