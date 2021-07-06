package ClusteredMotifs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map.Entry;

public class IdentifyMotifs {

	public static void getMotifs() {




		/* Search through annotation Files to get proteins annotated by */

	}

	/**
	 * Identify the motifs that are significantly clustered at the specified p-value threshold.
	 * (1) Output motif details to text file for easy look up
	 * (2) Store motifs and their file index in hash map for future use
	 * 
	 * @param motifClusteringPrefix		String - text file prefix for motifs and their clustering measure and significance
	 * @param numOfFiles				int - number of total files to search through
	 * @param pvalThreshold				double - significance threshold
	 * @param outputFile				String - text file path to output significant motifs that pass threshold
	 * 
	 * @return	motifsMapFileIdx		Map<String, Integer> - list of significant motifs and the file index that contains them
	 */
	public static HashMap<String, Integer> getSignificantMotifs(String motifClusteringPrefix, int numOfFiles, double pvalThreshold, String outputFile) {

		HashMap<String, Integer> motifsMapFileIdx = new HashMap<>();

		try {

			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("motif\tnProts\tTPD\tPval\t#File\n"); //header

			/* search all files for significant motifs */
			for(int i=0; i<numOfFiles; i++) {

				String motifClusteringFile = motifClusteringPrefix + i;

				InputStream in = new FileInputStream(new File(motifClusteringFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine();

				while(line!=null) {

					/* motif that pass significant threshold have their info printed to file and stored in map */
					double pval = Double.parseDouble(line.split("\t")[3]);
					if(pval <= pvalThreshold) {
						motifsMapFileIdx.put(line.split("\t")[0], i);
						out.write(line + "\t" + i + "\n");
					}
					line = input.readLine();
				}
				input.close();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifsMapFileIdx;
	}

	/**
	 * Find each significant motif in the annotation file and obtain the proteins in annotates, store in map. 
	 * 
	 * @param motifMapOfFileIdxs	Map<String, Integer> - map of motifs and their corresponding file idex
	 * @param annotationFilePrefix	String - file path prefix to annotation files
	 * 
	 * @return motifMapOfannotatedProteins	Map<String, String[]> - map of motif and it's list of annotated proteins
	 */
	public static HashMap<String, String[]> getAnnotatedProteinInfo(HashMap<String, Integer> motifMapOfFileIdxs, String annotationFilePrefix){

		HashMap<String, String[]> motifMapOfAnnotatedProteins = new HashMap<>();

		/* Search for every motif */
		for(Entry<String, Integer> m: motifMapOfFileIdxs.entrySet()) {

			String annotationFile = annotationFilePrefix + m.getValue();
			String motif = m.getKey();
			try {
				InputStream in = new FileInputStream(new File(annotationFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));
				
				String line = input.readLine();
				
				while(line!=null) {
					
					/* Store annotated proteins for the found motif && break out of loop */
					if(line.split("\t")[0].equals(motif)) {
						motifMapOfAnnotatedProteins.put(motif, line.split("\t")[2].split("|"));
						break;  
					}
					line = input.readLine();
				}
				
				input.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return motifMapOfAnnotatedProteins;
	}

}
