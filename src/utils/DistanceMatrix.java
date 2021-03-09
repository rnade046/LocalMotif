package utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

import graph.Interaction;
import graph.Protein;

public class DistanceMatrix {

	public static void computeDistanceMatrix(ArrayList<Interaction> interactionList, ArrayList<Protein> proteinList, String outputFile) {
		/* Compute distance matrix using the list of proteins in the network and the list of known interactions in the network.
		 * This distance matrix assumes a distance of 1 if proteins interact. Outputs file : ./IO_files/DistanceMatrix.txt */
		/* Generate a weighted distance matrix based on the fold change of proteins in the network  */
		
		double[][] distance_matrix = new double[proteinList.size()][proteinList.size()];
		
		for (int i=0; i<proteinList.size(); i++) {
			
			HashMap<String, Double> proteinsInteractingWithProt1Map = new HashMap<String, Double>(); // initialize Map to contain names of interacting proteins and their weight
			String prot1 = proteinList.get(i).getProteinName(); // get protein object
			
			/* Get all interacting proteins and their associated weight with prot1 > store in HashMap<String, Double> interacting protein : weight */ 
			for (int h = 0; h < interactionList.size(); h++) {
				Interaction ppi = interactionList.get(h); // get interaction
				
				if (ppi.getProtein1().equals(prot1)) {
					proteinsInteractingWithProt1Map.put(ppi.getProtein2(), ppi.getWeight());
				} else if (ppi.getProtein2().equals(prot1)) {
					proteinsInteractingWithProt1Map.put(ppi.getProtein1(), ppi.getWeight());
				}
			}

			/* Initialize distance matrix row for prot1*/
			for (int j = 0; j < distance_matrix.length; j++) {
			
				if (i == j) { // if it's the same proteins
					distance_matrix[i][j] = 0;
					continue;
				}

				if (proteinsInteractingWithProt1Map.containsKey(proteinList.get(j).getProteinName())) { // if protein i and j are connected
					distance_matrix[i][j] = proteinsInteractingWithProt1Map.get(proteinList.get(j).getProteinName()); // set strength
				} else { // otherwise set to max_value
					distance_matrix[i][j] = Double.MAX_VALUE;
				}
			}
		}
		
		/* Implement the floyd-warshall algorithm */
		for (int k = 0; k < distance_matrix.length; k++) {
			
			if(k%100 == 0) {
				System.out.println("k = " + k);
			}
			
			for (int i = 0; i < distance_matrix.length; i++) {
				for (int j = 0; j < distance_matrix.length; j++) {
					if (distance_matrix[i][j] > (distance_matrix[i][k] + distance_matrix[k][j])) {
						distance_matrix[i][j] = distance_matrix[i][k] + distance_matrix[k][j];
					}
				}
			}
		}

		//Print distance matrix
		try {
			BufferedWriter o = new BufferedWriter(new FileWriter(outputFile));

			for (int i = 0; i < distance_matrix.length; i++) {
				for (int j = 0; j < distance_matrix.length; j++) {
					o.write(distance_matrix[i][j] + "\t");
					o.flush();
				}
				o.write("\n");
				o.flush();
			}
			o.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public static void updateDistanceMatrix(boolean[] proteinsToKeep, double[][] distance_matrix, String filePath) {

		int protCount = 0;

		// Determine number of proteins that are connected to the network
		for(int i=0; i<proteinsToKeep.length; i++) {
			if(proteinsToKeep[i]) {
				protCount++;
			} 
		}

		//System.out.println("Number proteins to keep : " + protCount);

		double[][] updatedDistanceMatrix = new double[protCount][protCount];

		int row = 0;
		for(int i=0; i<distance_matrix.length; i++) {
			int col = 0;

			if(proteinsToKeep[i]) {
				for (int j=0; j<distance_matrix.length; j++){
					if(proteinsToKeep[i] == true && proteinsToKeep[j] == true) {
						//System.out.println("row=" +row +"; col=" +col +"; i=" +i +"; j=" +j);
						updatedDistanceMatrix[row][col] = distance_matrix[i][j];
						col++;
					}
				}
				row++;
			}

		}


		try {
			// Define new text file

			BufferedWriter out = new BufferedWriter(new FileWriter(new File(filePath)));

			for (int i = 0; i < updatedDistanceMatrix.length; i++) {
				for (int j = 0; j < updatedDistanceMatrix.length; j++) {

					out.write(updatedDistanceMatrix[i][j] + "\t");
				}
				out.write("\n");
				out.flush();
			}

			out.close();

		} catch (Exception e) {
			e.printStackTrace();
		}	
	}
}
