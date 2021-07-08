package fdr;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Properties;

public class MainFDR {

	public static void main(String[] args) throws FileNotFoundException, IOException {
		
		System.out.println("**Loading parameters file** \n");
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));		

		String wd = params.getProperty("working_directory");
		String projectName = params.getProperty("project_name");

		//String motifsPrefix = wd + "motifClustering/" + projectName + "_testedDegenMotifClustering_";
		//String nullMotifsPrefix = wd + "motifClustering/" + projectName + "_testedDegenMotifClustering_";

		int clusteringMeasure = Integer.parseInt(params.getProperty("clusteringMeasure", "0"));
		double percentThreshold = Double.parseDouble(params.getProperty("percentThreshold", "0.2"));

		String clusteringName = "";

		switch(clusteringMeasure) {
		case 0: clusteringName = "";
		break;
		case 1: clusteringName = "_TPPD" + percentThreshold;
		break;
		case 2: clusteringName = "_coreTPD" + percentThreshold;
		break;
		}
		
		String motifs_significanceScoresFile = wd + "fdr/" + projectName + clusteringName + "_listOfCalculatedSignificanceScores.tsv";
		String nullModel_significanceScoresFile = wd + "fdr/" + projectName + "_nullModel" + clusteringName + "_listOfCalculatedSignificanceScores.tsv";
		
		SimpleDateFormat sdf1 = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss");
		Timestamp timestamp = new Timestamp(System.currentTimeMillis());
		String fdrOutput = wd + "fdr/" + sdf1.format(timestamp) + "_" + projectName + clusteringName + "_FDRatThresholds_monotonicTransformation.tsv";
		//String significantMotifs = wd + pro
		
		/* Get list of significance scores */
		
		/* Compute FDRs between motifs and null model + monotonic transformation */
		System.out.println("**Initializing FDR calculation**");
		FdrCalculator fdrCalc = new FdrCalculator(motifs_significanceScoresFile, nullModel_significanceScoresFile);
		System.out.println("**Computing FDR calculation**");
		ArrayList<FalseDiscoveryRate> fdr = fdrCalc.computeFdr();
		
		/* Print information */
		ExportFdrDetails.exportFDR(fdr, fdrOutput);
	}

}
