package benchmark;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class FilterMemeOutput {

	public static void main(String[] args) {

		String memeInputFilePrefix = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\benchmark\\meme-zoops-results-2\\meme_cluster_";
		String memeFilteredFilePrefix = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\benchmark\\meme-zoops-results-2\\memeFiltered_cluster_";
		int numOfFilesToFilter = 104;
		double threshold = 0.05;

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

					String line = in.readLine();
					String stars = line;

					boolean print = true;

					boolean header = false;

					int starCount=0;
					while(line != null) {

						if(line.startsWith("*")) {
							starCount++;

							switch(starCount) {
							case 1: 			// start of header 
								header = true;
								print = true;

								line = in.readLine(); // next line

								if(line.startsWith("MOTIF")) {
									motifsTotal++;

									/* check e-val*/
									if(Double.parseDouble(line.split("\\s+")[14]) > threshold) {
										print = false;
									} else {
										motifsKept++;
									}
								}
								break;
							case 2 : 			// end of header
								header = false;
								break;
							case 3: 			// end of sections
								starCount = 0; 
								break;
							}
						}

						/* print lines */
						if(print) {
							if(header) {
								out.write(stars + "\n");
							}
							out.write(line +"\n");
						}
						line = in.readLine();
					}
					System.out.println("Motifs kept : " + motifsKept + " / " + motifsTotal);
					out.close();
					in.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			} else {
				System.out.println("No file");
			}
			
		}

	}

}
