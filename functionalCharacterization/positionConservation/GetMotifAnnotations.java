package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;

public class GetMotifAnnotations {

	public static void main(String[] args) {

		String motifFile = args[0];
		String proteinFreqFile = args[1];
		String annotationPrefix = args[2];
		int numFiles = 999;
		String outputFile = args[3];

		getAnnotations(annotationPrefix, numFiles, motifFile, proteinFreqFile, outputFile);
		
	}


	private static void getAnnotations(String annotationPrefix, int numOfFiles, String motifFile, String proteinAnnotationFreqFile, String outputFile) {
	
		HashSet<String> proteinsInNetwork = getProteinsInNetwork(proteinAnnotationFreqFile);
		HashSet<String> motifSet = loadMotifs(motifFile);
		System.out.println("Motifs to test: " + motifSet.size());
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=0; i<numOfFiles; i++) {
				
				if(i%50 == 0) {
					System.out.println();
				}
				System.out.print(i + ".");
				
				int motifCount = 0;
				InputStream in = new FileInputStream(new File(annotationPrefix + i));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine();

				while(line!=null) {

					String motif = line.split("\t")[0];
					if(motifSet.contains(motif)) {

						/* Check if protein was considered in analysis */
						ArrayList<String> annotatedProteinsInNetwork = new ArrayList<>();

						for(String prot : line.split("\t")[2].split("\\|")) { // all proteins annotated by motif in annotation file
							if(proteinsInNetwork.contains(prot)) {
								annotatedProteinsInNetwork.add(prot);
							}
						}
						
						/* output {motif	#proteins	proteinList} */
						out.write(line.split("\t")[0] + "\t" + annotatedProteinsInNetwork.size() + "\t");

						for(String prot: annotatedProteinsInNetwork) {
							out.write(prot + "|");
						}
						out.write("\n");
						out.flush();
						
						motifCount++;

						if(motifCount == motifSet.size()) {
							break;
						}
					}
					line = input.readLine();
				}
				input.close();
				
				if(motifCount == motifSet.size()) {
					break;
				}
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

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
	
	private static HashSet<String> loadMotifs(String motifFile) {
		HashSet<String> motifSet = new HashSet<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(motifFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();

			while(line!=null) {

				motifSet.add(line.split("\t")[1]); // [1] = motif
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifSet;
	}


}
