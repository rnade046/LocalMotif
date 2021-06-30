package opt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class FindSignificantMotifs {

	public static void main(String[] args) {

		String motifClusteringPrefix = "/home/rnade046/projects/rrg-mlaval/rnade046/motifs_Full_FWD/motifClustering/corrNetTop2_testedDegenMotifClustering_";
		int numOfFiles = 999;

		double pvalThreshold = 1.88371083687543E-08;
		
		String outputFile = "/home/rnade046/projects/rrg-mlaval/rnade046/motifs_Full_FWD/corrNetTop2_signficantMotifs_p" + pvalThreshold +".tsv";

		getSignificantMotifs(motifClusteringPrefix, numOfFiles, pvalThreshold, outputFile);

	}

	public static void getSignificantMotifs(String motifClusteringPrefix, int numOfFiles, double pvalThreshold, String outputFile) {


		try {
			
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(int i=0; i<numOfFiles; i++) {

				String motifClusteringFile = motifClusteringPrefix + i;

				InputStream in = new FileInputStream(new File(motifClusteringFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine();
				
				while(line!=null) {
					
					double pval = Double.parseDouble(line.split("\t")[3]);
					if(pval <= pvalThreshold) {
						out.write(line + "\n");
					}
					
					line = input.readLine();
				}

				input.close();
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
