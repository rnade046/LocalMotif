import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.util.ArrayList;

public class MainFDR {

	public static void main(String[] args) {
		
		/* input arguments */
		String wd = args[0];
		String motifs_significanceScoresFile = wd + args[1];
		String nullModel_significanceScoresFile =  wd + args[2];
		
		/* output file */
		SimpleDateFormat sdf1 = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss");
		Timestamp timestamp = new Timestamp(System.currentTimeMillis());
		String fdrOutput = wd + sdf1.format(timestamp) + "_FDRatThresholds_monotonicTransformation.tsv";
		
		/* Compute FDRs between motifs and null model + monotonic transformation */
		System.out.println("Initializing FDR calculation");
		FdrCalculator fdrCalc = new FdrCalculator(motifs_significanceScoresFile, nullModel_significanceScoresFile);
		System.out.println("Computing FDR at various p-value thresholds");
		ArrayList<FalseDiscoveryRate> fdr = fdrCalc.computeFdr();
		
		/* Print information */
		System.out.println("FDR calculations are stored: " + fdrOutput);
		ExportFdrDetails.exportFDR(fdr, fdrOutput);
		System.out.println("--- Completed FDR estimation ---");
	}

}
