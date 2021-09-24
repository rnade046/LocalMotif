package ClusteredMotifs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

public class IdentifyMotifs {


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
	public static void getSignificantMotifs(String motifClusteringPrefix, int numOfFiles, double pvalThreshold, String outputFile) {

		try {

			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("motif\tnProts\tTPD\tPval\t#File\n"); //header

			/* search all files for significant motifs */
			for(int i=0; i<numOfFiles; i++) {
				
				if(i%10 == 0) {
					System.out.print(i + ".");
				}
				
				if(i%100 == 0) {
					System.out.println();
				}
				
				String motifClusteringFile = motifClusteringPrefix + i;

				InputStream in = new FileInputStream(new File(motifClusteringFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine();

				while(line!=null) {

					/* motif that pass significant threshold have their info printed to file and stored in map */
					double pval = Double.parseDouble(line.split("\t")[3]);
					if(pval <= pvalThreshold) {
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
		
		System.out.println("Done");
	}
	
	
	public static HashMap<String, Integer> loadSignificantMotifs(String motifClusteringFile) {

		HashMap<String, Integer> motifsMapFileIdx = new HashMap<>();

		try {

				InputStream in = new FileInputStream(new File(motifClusteringFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // header
				line = input.readLine();

				while(line!=null) {

					motifsMapFileIdx.put(line.split("\t")[0], Integer.parseInt(line.split("\t")[4])); // [0] = motif, [4] = file number
					
					line = input.readLine();
				}
				input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		System.out.println("Done");
		return motifsMapFileIdx;
	}
	

	/**
	 * Find each significant motif in the annotation file and obtain the proteins in annotates, store in map. 
	 * 
	 * @param motifMapOfFileIdxs	Map<String, Integer> - map of motifs and their corresponding file index
	 * @param annotationFilePrefix	String - file path prefix to annotation files
	 * 
	 * @return motifMapOfannotatedProteins	Map<String, String[]> - map of motif and it's list of annotated proteins
	 */
	public static void getAnnotatedProteinInfo(HashMap<String, Integer> motifMapOfFileIdxs, String proteinAnnotationFreqFile, String outputFile, String annotationFilePrefix){
		
		HashSet<String> proteinsInNetwork = getProteinsInNetwork(proteinAnnotationFreqFile);
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			/* Search for every motif */
			int motifCount = 0;
			for(Entry<String, Integer> m: motifMapOfFileIdxs.entrySet()) {
				motifCount++;
				
				if(motifCount%10==0) {
					System.out.print(motifCount + ".");
				}
				
				if(motifCount%100==0) {
					System.out.println();
				}
				
				String annotationFile = annotationFilePrefix + m.getValue();
				String motif = m.getKey();

				InputStream in = new FileInputStream(new File(annotationFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine();

				while(line!=null) {

					/* Store annotated proteins for the found motif && break out of loop */
					if(line.split("\t")[0].equals(motif)) {
						
						/* Check if protein was considered in analysis */
						ArrayList<String> annotatedProteinsInNetwork = new ArrayList<>();
						
						for(String prot : line.split("\t")[2].split("\\|")) { // all proteins annotated by motif in annotation file
							if(proteinsInNetwork.contains(prot)) {
								annotatedProteinsInNetwork.add(prot);
							}
						}
						
						/* output {motif	#proteins	proteinList} */
						out.write(motif + "\t" + annotatedProteinsInNetwork.size() + "\t");
						
						for(String prot: proteinsInNetwork) {
							out.write(prot + "|");
						}
						out.write("\n");
						out.flush();
						break;  
					}
					line = input.readLine();
				}

				input.close();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Done");
	}

	private static HashSet<String> getProteinsInNetwork(String protFile){
	
		HashSet<String> proteinSet = new HashSet<>();
		
		try {

			InputStream in = new FileInputStream(new File(protFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String line = input.readLine();

			while(line!=null) {
			
				proteinSet.add(line.split("\t")[0]);
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		return proteinSet;
	}
	
	public static HashMap<String, String[]> loadAnnotatedProteinInfo(String annotationFile){

		HashMap<String, String[]> motifMapOfAnnotatedProteins = new HashMap<>();

		try {

			InputStream in = new FileInputStream(new File(annotationFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String line = input.readLine();

			while(line!=null) {
			
				String[] col = line.split("\t");
				motifMapOfAnnotatedProteins.put(col[0], col[2].split("\\|"));
				
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Done");
		return motifMapOfAnnotatedProteins;
	}
}
