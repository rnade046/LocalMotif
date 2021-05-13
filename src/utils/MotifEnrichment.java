package utils;

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

import org.apache.commons.math3.distribution.NormalDistribution;

import graph.Protein;

public class MotifEnrichment {

	private double[][] distanceMatrix;
	private HashMap<String, Integer> indexOfProteinsInNetwork;
	private HashMap<Integer, double[]> normalDistributionParams;
	private int lowerBound;
	private int upperBound;

	public MotifEnrichment(double[][] _distance_matrix, ArrayList<Protein> proteinsInNetwork, String normalDistributionParamsFile, int _lowerBound, int _upperBound) {

		this.distanceMatrix = _distance_matrix;
		this.indexOfProteinsInNetwork = getIndexOfProteins(proteinsInNetwork);

		this.normalDistributionParams = loadNormalDistributionParams(normalDistributionParamsFile);
		this.lowerBound = _lowerBound;
		this.upperBound = _upperBound;
	}

	public void testMotifClustering(String annotationFilePrefix, String outputPrefix, int currentLowerBoundToTest, int currentUpperBoundToTest) {

		for(int i=currentLowerBoundToTest; i<=currentUpperBoundToTest; i++) {
			
			/* Load annotations 1 line at a time; 
			 * 1st check: if num of proteins >= 3 or <= 2000
			 * 2nd check: ensure proteins associated to annotation are in the network 
			 * 3rd check: confirm 1st check */ 
			String annotationFile = annotationFilePrefix + i;
			String annotationOutputFile = outputPrefix + i;
			
			try {
				InputStream in = new FileInputStream(new File(annotationFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // no header
				int motifCount = 1;
				while(line!=null) {

					if(motifCount%1000 == 0) {
						System.out.print(motifCount+".");
					}
					
					if(motifCount%10000 == 0) {
						System.out.println();
					}
					
					String col[] = line.split("\t");
					
					/* check the number of proteins annotated by motifs is within our set bounds of n proteins to be tested */
					int numProt = Integer.parseInt(col[2]);
					if(numProt>=lowerBound && numProt<=upperBound) {

						/* get list of proteins in Network */
						ArrayList<String> proteinInNetworkAssociatedToMotif = getListOfProteinsInNetwork(col[4].split("\\|"));
						
						/* check againt that the number of proteins annotated by motifs is within our set bounds of n proteins to be tested */
						if(proteinInNetworkAssociatedToMotif.size()>=lowerBound && proteinInNetworkAssociatedToMotif.size() <= upperBound) {
							
							/* compute TPD */
							ArrayList<Integer> indexOfProteinsInNetworkAnnotatedByMotif = getIndexOfProteinsInNetworkAnnotatedByMotif(proteinInNetworkAssociatedToMotif);
							double tpd = Calculator.computeTPD(distanceMatrix, indexOfProteinsInNetworkAnnotatedByMotif);
							
							/* assess clustering significance */
							double[] params = normalDistributionParams.get(proteinInNetworkAssociatedToMotif.size());
							NormalDistribution nd = new NormalDistribution(params[0], params[1]);
							
							double p_val = nd.probability(0, tpd);
							
							/* print motif; numProts; listProteins; tpd; p-val */
							printAnnotationDetails(annotationOutputFile, col[0], proteinInNetworkAssociatedToMotif, tpd, p_val);
						}
					}

					line = input.readLine();
					motifCount++;
				}

				input.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}


	/**
	 * Given the list of proteins in the network, generate an a map of it's index (ie. the protein = 1)
	 * where 1 corresponds to it's index in the distance matrix
	 * @param proteinsInNetworkList
	 * @return indexOfProteinsInNetwork		Map<String, Integer> - {protein name = index in distance matrix}
	 */
	private static HashMap<String, Integer> getIndexOfProteins(ArrayList<Protein> proteinsInNetworkList){
		HashMap<String, Integer> indexOfProteinsInNetwork = new HashMap<>();

		for(int i=0; i<proteinsInNetworkList.size(); i++) {
			indexOfProteinsInNetwork.put(proteinsInNetworkList.get(i).getProteinName(), i);
		}

		return indexOfProteinsInNetwork;
	}

	/**
	 * Load normal distributions calculated from Monte Carlo Distribution from file. Store in map to access by key = number of proteins.
	 * 
	 * @param normalDistributionParamsFile	String - file path to calculated parameters
	 * @return ndParams 	HashMap<Integer, double[]> - map number of proteins = [mean, stdev]
	 */
	private static HashMap<Integer, double[]> loadNormalDistributionParams(String normalDistributionParamsFile){
		HashMap<Integer, double[]> ndParams = new HashMap<>();

		try {
			InputStream in = new FileInputStream(new File(normalDistributionParamsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();

			while(line!=null) {
				String[] col = line.split("\t"); // [0] = number of proteins, [1] = mean, [2] = standard deviation

				double[] params = new double[] {Double.parseDouble(col[1]), Double.parseDouble(col[2])};
				ndParams.put(Integer.parseInt(col[0]), params);

				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return ndParams;
	}
	
	/**
	 * Check that proteins associated to motifs are indeed found in our network. Return the updated list of 
	 * proteins found in network associated to the motif whose clustering is being tested. 
	 * 
	 * @param proteinsAssociatedToMotif		String[] - list of protein associated to the motif (from annotation file)
	 * @return proteinsInNetworkAssociatedToMotifList	List<String> - list of proteins found in network and associated to the motif 
	 */
	private ArrayList<String> getListOfProteinsInNetwork(String[] proteinsAssociatedToMotif){

		ArrayList<String> proteinsInNetworkAssociatedToMotifList = new ArrayList<>();

		for(String prot: proteinsAssociatedToMotif) {
			if(this.indexOfProteinsInNetwork.containsKey(prot)) {
				proteinsInNetworkAssociatedToMotifList.add(prot);
			}
		}

		return proteinsInNetworkAssociatedToMotifList;
	}
	
	/**
	 * Get the indexes (associated to the distance matrix) for all the proteins in the network annotated by the tested motif.
	 * @param listOfProteinsInNetworkAnnotatedByMotif
	 * @return 
	 */
	private ArrayList<Integer> getIndexOfProteinsInNetworkAnnotatedByMotif(ArrayList<String> listOfProteinsInNetworkAnnotatedByMotif){
		
		ArrayList<Integer> listOfIndexOfProteinsInNetworkAnnotatedByMotif = new ArrayList<>();
		
		for(String prot: listOfProteinsInNetworkAnnotatedByMotif) {
			listOfIndexOfProteinsInNetworkAnnotatedByMotif.add(indexOfProteinsInNetwork.get(prot));
		}
		
		return listOfIndexOfProteinsInNetworkAnnotatedByMotif;
	}
	
	private void printAnnotationDetails(String annotationOutputFile, String motif, ArrayList<String> proteinsInNetworkAnnotatedByMotif, double tpd, double pValue) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(annotationOutputFile), true));
			
			out.write(motif + "\t" + proteinsInNetworkAnnotatedByMotif.size() + "\t");
			
			for(String prot:proteinsInNetworkAnnotatedByMotif) {
				out.write(prot + "|");
			}
			
			out.write("\t" + tpd + "\t" + pValue + "\n");
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

	
}
