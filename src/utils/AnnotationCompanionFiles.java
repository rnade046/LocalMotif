package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;

public class AnnotationCompanionFiles {

	/**
	 * Load annotation file; check proteins associated to annotations 1-by-1; 
	 * if the actual number of proteins lies between the lower and upper bound, 
	 * print the motif, the number of originally listed proteins and the actual number of proteins (in this network)
	 * 
	 * @param annotationPrefix
	 * @param companionFilePrefix
	 * @param proteinSet
	 * @param numFiles
	 * @param lowerBound
	 * @param upperBound
	 */
	public static void assessAnnotationFile(String annotationPrefix, String companionFilePrefix, HashSet<String> proteinSet, int numFiles, int lowerBound, int upperBound) {

		for(int i=0; i<numFiles; i++) {

			try {
				InputStream in = new FileInputStream(new File(annotationPrefix + i));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				BufferedWriter out = new BufferedWriter(new FileWriter(new File(companionFilePrefix + i)));

				String line = input.readLine(); // no header

				while(line!= null) {

					String[] col = line.split("\t");
					String motif = col[0];
					String[] proteins = col[2].split("\\|");

					/* Count proteins in network */
					int protCount = 0;
					for(String protein: proteins) {
						if(proteinSet.contains(protein)) {
							protCount++;
						}
					}

					/* print motif info if it respects upper and lower bound */
					if(protCount >= lowerBound && protCount <= upperBound) {
						out.write(motif + "\t" + proteins.length + "\t" + protCount + "\n");
						out.flush();
					}

					line = input.readLine();
				}
				input.close();
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

		}
	}
}
