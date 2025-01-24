package summarizeMotifs;
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

import ClusteredMotifs.Family;

public class SummarizeSignificantMotifs {

	public static void summarizeMotifs(ArrayList<Family> motifs, String fdrInfoFile, String proteinAnnotationFile, String significantMotifFile, String outputFile) {

		/* load FDR info as ascending list {p-value, FDR} */
		System.out.println("** load FDR info **");
		List<Double[]> fdrInfo = loadFDRinfo(fdrInfoFile);

		/* load protein list */
		System.out.println("** load proteins **");
		HashMap<String, String> proteinMap = loadProteinList(proteinAnnotationFile);

		/* load significant motifs - store as Motif element - combine info */
		System.out.println("** load significant motifs **");
		List<ClusteredMotif> significantMotifs = loadSignificantMotifs(significantMotifFile, motifs, fdrInfo, proteinMap);

		/* print list */
		System.out.println("** print motif summary **");
		printMotifSummary(significantMotifs, outputFile);
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

	private static List<ClusteredMotif> loadSignificantMotifs(String significantMotifsFile, ArrayList<Family> motifs, List<Double[]> fdrInfo, HashMap<String, String> proteinMap){

		List<ClusteredMotif> significantMotifs = new ArrayList<>();

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
				significantMotifs.add(new ClusteredMotif(col, motifs, fdrInfo, proteinMap));

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

	private static void printMotifSummary(List<ClusteredMotif> significantMotifs, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			/* header */
			out.write("Motif\tFamilyNumber\tRepresentative\tNumberOfProteins\tClusteringMeasure\tp-value\tFDR\tProteinList\n");
			
			/* table */
			for(ClusteredMotif m: significantMotifs) {
				
				out.write(m.getMotif() + "\t" + m.getFamily() +"\t" + m.isRepresentative() + "\t" + m.getNumberOfProteins()+ "\t" + m.getClusteringMeasure() + "\t" 
						+ m.getPval() + "\t" + m.getFDR() + "\t" + m.getProteinList() +"\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


	}
}
