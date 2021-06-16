package fdr;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;

public class MainFDR {

	public static void main(String[] args) throws FileNotFoundException, IOException {
		
		System.out.println("**Loading parameters file** \n");
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));		

		String wd = params.getProperty("working_directory");
		String projectName = params.getProperty("project_name");

		String motifsPrefix = wd + "motifClustering/" + projectName + "_testedDegenMotifClustering_";
		String nullMotifsPrefix = wd + "motifClustering/" + projectName + "_testedDegenMotifClustering_";

		String motifs_significanceScoresFile = wd + projectName + "_listOfCalculatedSignificanceScores.tsv";
		String nullModel_significanceScoresFile = wd + projectName + "_listOfCalculatedSignificanceScores.tsv";

		/* Compute FDRs between motifs and null model + monotonic transformation */
		FdrCalculator fdrCalc = new FdrCalculator(motifs_significanceScoresFile, nullModel_significanceScoresFile);
		ArrayList<FalseDiscoveryRate> fdr = fdrCalc.computeFdr();
		
		
		
	}

}
