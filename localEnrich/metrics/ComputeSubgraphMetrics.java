package metrics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


public class ComputeSubgraphMetrics {

	public static void main(String[] args) {
		
		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";
		String motifFamiliesFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv";
		String annotationFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteinsByMotif.tsv";
		String proteinIdxFile = wd + "corrNetTop2-400_proteinsInNetwork_info.tsv";
		
		String dmFile = wd + "corrNetTop2-400_removedOverConnectedProteins_400_distanceMatrix2.txt";
		String outputFile = wd + "network-topology/diameter/coreTPD0.4_s100.000_motifFamilies_diameter.tsv";
		
		/* load significant motifs */
		HashMap<String, Integer> motifFamilies = loadMotifFamilies(motifFamiliesFile);
		
		/* obtain list of proteins annotated by motifs */
		HashMap<String, String[]> proteinMap = loadProteinList(annotationFile, motifFamilies);
		
		/* obtain indexes of these proteins in shortest path matrix */
		HashMap<String, Integer> proteinIdxMap = loadProteinIdxes(proteinIdxFile);
		
		/* load distance matrix */
		double[][] dm = loadDistanceMatrix(dmFile, proteinIdxMap.size());
		
		/* calculate diameter = longest shortest path for subgraph */
		System.out.println("Motifs");
		List<Motif> motifList = new ArrayList<>();
		for(String motif: motifFamilies.keySet()) {
		
			motifList.add(new Motif(motif, motifFamilies.get(motif), proteinMap.get(motif), proteinIdxMap, dm));
			System.out.println(motifList.size());
		}
		/* calculate degree centrality for each node of the subgraph = number of connections in the subgraph
		 * need to input original network */
		
		printInfo(outputFile, motifList);
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
	
	private static HashMap<String, String[]> loadProteinList(String proteinAnnotationFile, HashMap<String, Integer> motifMap){

		HashMap<String, String[]> proteinMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinAnnotationFile))));

			String line = in.readLine(); // no header

			while(line != null) {

				String[] col = line.split("\t");
				if(motifMap.containsKey(col[0])) {
					proteinMap.put(col[0], col[2].split("\\|"));
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinMap;	
	}
	
	private static HashMap<String, Integer> loadProteinIdxes(String proteinOrderFile){

		HashMap<String, Integer> proteinMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinOrderFile))));

			String line = in.readLine(); // no header
			int lineCount = 0;
			
			while(line != null) {

				proteinMap.put(line.split("\t")[0], lineCount);
				line = in.readLine();
				lineCount++;
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinMap;	
	}
	
	public static double[][] loadDistanceMatrix(String distance_matrixFile, int numProts) {
		/* Import distance matrix from text file, it's dimensions are based on the size of the 
		 * proteinNetwork List initially used to build the distance Matrix file */

		double[][] distanceMatrix = new double[numProts][numProts]; // Initialize distance matrix

		try {

			InputStream in = new FileInputStream(new File(distance_matrixFile));				
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // read first line

			int x = 0; // initialize row counter

			while (line != null) {

				String[] col = line.split("\t"); // split columns
				int y = 0; // initialize column counter (resets at the end of every row)

				for (String str : col) { // str is the index to go through all element of col

					distanceMatrix[x][y] = Double.parseDouble(str);;// set the value (distance) at the appropriate coordinates
					y++; // adds one to the value of y (change column)
				}
				x++; // adds one to the vale of x (change row)
				line = input.readLine(); // read next line

			}
			input.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		return distanceMatrix;
	} // end import distance matrix
	
	private static void printInfo(String outputFile, List<Motif> motifList) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			/* header */
			out.write("Motif\tFamily#\t#Proteins\tDiameter\n");

			/* table */
			for(Motif m: motifList) {
				out.write(m.getMotif() + "\t" + m.getFamilyNumber() + "\t" + m.getNumberOfProteins()+ "\t" + m.getDiameter() +"\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		
	}
}
