package opt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class GetMotifsThatPassSignificanceThreshold {

	public static void main(String[] args) {

		// args[0] = file path
		// args[1] = threshold

		double threshold = Double.parseDouble(args[1]);
//		double threshold = 0.05;
		
		String inputFolder = args[0];
//		String inputFolder = "";
		String outputFile = "significantMotis_at_p" + threshold + ".tsv";

		/* iterate over files in motifClustering/ */
		File dir = new File(inputFolder);
		File[] directoryListing = dir.listFiles();
		
		int countFiles = 0;
		int countMotifs = 0;
		if (directoryListing != null) {
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
				for (File f : directoryListing) {
					System.out.println("searching file: " + f);
					countFiles++;
					BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(f)));

					String line = in.readLine();
					while (line!=null) {
						// line 1: AYAGCCTA        169     2807.637        0.9484101824541571
						// col[0] = motif; col[4] = p-val
						String[] col=line.split("\t");
						if(Double.parseDouble(col[3]) <= threshold) {
							out.write(col[0] + "\t" + col[3]);
							out.flush();
						}
						line = in.readLine();
						countMotifs++;
					}
					in.close();
				}
				out.close();
				
				System.out.println("totalFiles: " + countFiles);
				System.out.println("totalMotifs: " + countMotifs);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
}
