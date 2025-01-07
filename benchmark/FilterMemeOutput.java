

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map.Entry;

public class FilterMemeOutput {

	public static void main(String[] args) {

		String memeInputFilePrefix = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/benchmark/MEME_de/meme_";
		String memeFilteredFilePrefix = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/benchmark/MEME_de/memeFiltered_0.1_";

		String summaryFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/benchmark/MEME_de/meme_anr_summary_0.1.txt";
		int numOfFilesToFilter = 118;
		double threshold = 0.1;

		HashMap<Integer, Integer> significantMotifsMap = new HashMap<>();

		/* filter every MEME.txt file output by MEME analysis on MCL clusters */
		for(int i=1; i<= numOfFilesToFilter; i++) {
			System.out.println("File " + i);

			/* check if file exits */
			File f = new File(memeInputFilePrefix + i + ".txt");
			if(f.exists() && !f.isDirectory()) { 

				int motifsTotal = 0;
				int motifsKept = 0;

				try {
					BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(memeInputFilePrefix + i + ".txt"))));
					BufferedWriter out = new BufferedWriter(new FileWriter(new File(memeFilteredFilePrefix + i + ".txt")));

					/* Header */
					out.write("MEME version 4\n\n" + "ALPHABET= ACGT\n\n" + 
							"Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"); // header

					String line = in.readLine();

					boolean printMotif = false;
					String motif="";
					int starCount = 0;
					while(line != null) {

						if(line.startsWith("*")) {

							starCount++;

							switch(starCount) {
							case 1: 			// start of header 
								line = in.readLine(); // next line
								if(line.startsWith("MOTIF")) {
									motifsTotal++;
									/* check E-val*/
									if(Double.parseDouble(line.split("\\s+")[14]) < threshold) {
										motif = line.split("\\s+")[2];
										printMotif = true;
										motifsKept++;
										
										starCount++;
									}
								}
								break;
							case 3: 			// end of sections
								starCount = 0; 
								break;
							}
						}

						/* print lines */
						if(printMotif) {
							boolean ppmFound = false;

							while(!ppmFound) {

								if(line.contains("position-specific probability matrix")) {
									ppmFound = true;

									out.write("MOTIF " + motif + "\n");
									line = in.readLine();
								}
								line = in.readLine();
							}

							while (!line.contains("---")) {
								out.write(line + "\n");
								out.flush();
								line = in.readLine();
							}
							
							out.write("\n");
							printMotif=false; 
						}
						line = in.readLine();
					}
					System.out.println("Motifs kept : " + motifsKept + " / " + motifsTotal);

					significantMotifsMap.put(i, motifsKept);
					out.close();
					in.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			} else {
				System.out.println("No file");
				significantMotifsMap.put(i, 0);
			}
		}
		printMEMEmotifSummary(summaryFile, significantMotifsMap);

	}

	private static void printMEMEmotifSummary(String outputFile, HashMap<Integer, Integer> significantMotifMap) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("Cluster\tSignificantMotifs\n"); 

			for(Entry<Integer, Integer> e : significantMotifMap.entrySet()) {
				if(e.getValue() != 0) {
					out.write(e.getKey() + "\t" + e.getValue() + "\n");
				}
			
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
