package summarize;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class SummarizeSignificantMotifs {

	public static void main(String args[]) {

		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";

		String motifFamilyFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv";
		String significantMotifFile = wd + "corrNetTop2-400_coreTPD_p0.4_signficantMotifs_p6.94181641095822E-7.tsv";
		String proteinAnnotationFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteinsByMotif.tsv";

		String fdrInfoFile = wd + "/fdr/corrNet2-400/coreTPD/2022.03.08.11.22.43_corrNetTop2-400_coreTPD_p0.4_FDRatThresholds_monotonicTransformation_s100.000.tsv";
		String strandSpecificityFile = wd + "/MotifPosition/coreTPD0.4/StrandSpecificity_corrNet2-400_coreTPD_p0.4_3UTR_allMotifs.tsv";

		String outputFile = wd + "/fdr/corrNet2-400_coreTPD_0.4_s100.000_motifSummary.tsv";

		/* load motif family info map {motif = number} */
		System.out.println("** load motif families **");
		HashMap<String, Integer> motifMap = loadMotifFamilies(motifFamilyFile);

		/* load FDR info as ascending list {p-value, FDR} */
		System.out.println("** load FDR info **");
		List<Double[]> fdrInfo = loadFDRinfo(fdrInfoFile);

		/* load strand specificity */
		System.out.println("** load strand specificity **");
		HashMap<String, Double[]> strandSpecificityMap = loadStrandSpecificity(strandSpecificityFile);

		/* load protein list */
		System.out.println("** load proteins **");
		HashMap<String, String> proteinMap = loadProteinList(proteinAnnotationFile);

		/* load significant motifs - store as Motif element - combine info */
		System.out.println("** load significant motifs **");
		List<Motif> significantMotifs = loadSignificantMotifs(significantMotifFile, motifMap, fdrInfo, strandSpecificityMap, proteinMap);

		/* print list */
		System.out.println("** print motif summary **");
		printMotifSummary(significantMotifs, outputFile);

	}

	private static HashMap<String, Integer> loadMotifFamilies(String motifFamilyFile){

		HashMap<String, Integer> motifMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(motifFamilyFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {

				String[] col = line.split("\t");
				motifMap.put(col[0], Integer.parseInt(col[1]));

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifMap;
	}

	private static List<Double[]> loadFDRinfo(String fdrInfoFile){

		List<Double[]> fdrInfo = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fdrInfoFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {

				String[] col = line.split("\t");
				fdrInfo.add(new Double[] {Double.parseDouble(col[1]), Double.parseDouble(col[0])});

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return fdrInfo;
	}

	private static HashMap<String, Double[]> loadStrandSpecificity(String ssFile){

		HashMap<String, Double[]> ssMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(ssFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {

				String[] col = line.split("\t");
				ssMap.put(col[0], new Double[] {Double.parseDouble(col[3]), Double.parseDouble(col[5])});

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return ssMap;

	}

	private static HashMap<String, String> loadProteinList(String proteinAnnotationFile){

		HashMap<String, String> proteinMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinAnnotationFile))));

			String line = in.readLine(); // no header

			while(line != null) {

				String[] col = line.split("\t");
				proteinMap.put(col[0], col[2]);
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinMap;	
	}

	private static List<Motif> loadSignificantMotifs(String significantMotifsFile, HashMap<String, Integer> motifFamily, List<Double[]> fdrInfo, HashMap<String, Double[]> ssMap,
			HashMap<String, String> proteinMap){

		List<Motif> significantMotifs = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(significantMotifsFile))));

			String line = in.readLine(); // header 
			line = in.readLine();

			int motifCount = 1; 

			while(line != null) {

				System.out.print(motifCount + ".");

				if(motifCount%50 == 0) {
					System.out.println();
				}

				String[] col = line.split("\t");
				significantMotifs.add(new Motif(col, motifFamily, fdrInfo, ssMap, proteinMap));

				line = in.readLine();
				motifCount++;
			}
			in.close();
			System.out.println();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return significantMotifs;
	} 

	private static void printMotifSummary(List<Motif> significantMotifs, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			/* header */
			out.write("Motif\tFamily#\t#Proteins\tClusteringMeasure\tp-value\tFDR\tStrandSpecificity\tStrandSpecificity(adj.pval)\tProteinList\n");

			/* table */
			for(Motif m: significantMotifs) {
				out.write(m.getMotif() + "\t" + m.getFamily() + "\t" + m.getNumberOfProteins()+ "\t" + m.getClusteringMeasure() + "\t" 
						+ m.getPval() + "\t" + m.getFDR() + "\t" + m.getStrandSpecificity() + "\t" + m.getStrandSpecificitySignificance()+ "\t" + m.getProteinList() +"\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


	}
}
