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
import java.util.Map.Entry;

public class IdentifyCoreProteins {

	public static void getCoreProteins(String annotationSubsetFile, double percentThreshold, String proteinInfoFile, String distanceMatrixFile, String corePorteinsFile) {
		
		/* Load protein info file (the proteins in this file are in the order of the distance matrix) */
		ArrayList<String> proteinsInNetwork = loadProteinsInNetwork(proteinInfoFile);
		
		/* Load distance matrix file */
		double[][] distanceMatrix = loadDistanceMatrix(distanceMatrixFile, proteinsInNetwork);
		
		/* Load annotation subset line by line */ 
		try {

			InputStream in = new FileInputStream(new File(annotationSubsetFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // no header

			while(line!=null) {

				String[] col = line.split("\t");
				String[] allProteins = col[2].split("\\|");
				
				/* Identify number of core nodes */
				int nodeThreshold = (int) Math.ceil(allProteins.length*percentThreshold); // determine number of core proteins
				
				/* get protein indexes in distance matrix */ 
				ArrayList<Integer> annotatedProtsIdxs = getAnnotatedProteinIndexes(proteinsInNetwork, allProteins);
				
				ArrayList<Integer> coreProteinsIdxs = identifyCoreProteins(annotatedProtsIdxs, distanceMatrix, nodeThreshold); 	// Identify core proteins
				
				/* print info */ 
				printCoreProteinsInfo(corePorteinsFile, col[0], nodeThreshold, coreProteinsIdxs, proteinsInNetwork);
				
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static ArrayList<String> loadProteinsInNetwork(String inputFile){
		
		ArrayList<String> proteinsInNetworkList = new ArrayList<>();
		
		try {
			InputStream in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // no header

			while(line!=null) {

				proteinsInNetworkList.add(line.split("\t")[0]);
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinsInNetworkList;
	}
	
	private static double[][] loadDistanceMatrix(String distance_matrixFile, ArrayList<String> ProteinList) {
		/* Import distance matrix from text file, it's dimensions are based on the size of the 
		 * proteinNetwork List initially used to build the distance Matrix file */

		double[][] distanceMatrix = new double[ProteinList.size()][ProteinList.size()]; // Initialize distance matrix

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
	
	private static ArrayList<Integer> getAnnotatedProteinIndexes(ArrayList<String> proteinsInNetwork, String[] annotatedProteins){
		
		ArrayList<Integer> idxsOfAnnotatedProteins = new ArrayList<>();
		
		for(String protein: annotatedProteins) {
			
			for(int i=0; i<proteinsInNetwork.size(); i++) {
				if(proteinsInNetwork.get(i).equals(protein)) {
					idxsOfAnnotatedProteins.add(i);
					break;
				}
			}
		}
		return idxsOfAnnotatedProteins;
	}
	
	/**
	 * Identify core proteins from a list of proteins in the network. Core proteins are the top % nodes with the lowest sum of pairwise distance 
	 * with all other annotated proteins. The list of core proteins indexes in the network is returned.
	 * 
	 * @param proteinIndexes	List<Integer> index of annotated proteins in the network
	 * @param distanceMatrix	double[][]	distance matrix of shortest paths
	 * @param nodeThreshold		int threshold for protein core
	 * 
	 * @return coreProteinsList	List<Integer> indexes of proteins in the core (most clustered annotated proteins)
	 */
	private static ArrayList<Integer> identifyCoreProteins(ArrayList<Integer> proteinIndexes, double[][] distanceMatrix, int nodeThreshold){
		ArrayList<Integer> coreProteinsList = new ArrayList<>(); 

		/* Calculate the sum of interaction weights for the protein to all proteins annotated */
		HashMap<Integer, Double> mapOfSumInteractionWeights = new HashMap<>();
		for(int proteinIdx : proteinIndexes) {
			double sumOfWeights = 0;

			for(int j: proteinIndexes) {
				sumOfWeights += distanceMatrix[proteinIdx][j];			
			}

			mapOfSumInteractionWeights.put(proteinIdx, sumOfWeights);
		}
		/* Sort map by value */
		
		List<Entry<Integer, Double>> rankedProteinslist = new ArrayList<>(mapOfSumInteractionWeights.entrySet());
		rankedProteinslist.sort(Entry.comparingByValue());
		
		/* Identify top % nodes */ 
		for(int i=0; i<nodeThreshold; i++) {
			coreProteinsList.add(rankedProteinslist.get(i).getKey());
		}
		
		return coreProteinsList;
	}
	
	private static void printCoreProteinsInfo(String outputFile, String motif, int numProts, ArrayList<Integer> annotatedProtsIdxs, ArrayList<String> proteinsInNetwork) {
		
		try {
			
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile), true));

			out.write(motif + "\t" + numProts + "\t");
			for(int protIdx : annotatedProtsIdxs) {
				out.write(proteinsInNetwork.get(protIdx) + "|");
				out.flush();
			}
			out.write("\n");
			
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
}

