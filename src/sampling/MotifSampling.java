package sampling;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ThreadLocalRandom;

import graph.Protein;
import utils.TopPercentPairwiseDistance;

public class MotifSampling {
	/**
	 * Weighted motif sampling inspired by: https://stackoverflow.com/questions/1761626/weighted-random-numbers
	 * Here, we generate a list of cumulative weights to sample from. 
	 * We generate a random number between 0 and the sum of cumulative weights. 
	 * The selected protein will be the one where its assigned cumulative weight is greater or equal to the random number.
	 */
	double[][] distanceMatrix;
	ArrayList<Protein> proteinsInNetworkList;
	HashMap<String, Integer> proteinToOccurrenceMap; 	// Map containing annotated proteins and their occurrence of annotation
	long[] cumulativeSumOfWeights;						// array containing the cumulative sum of weights for random weighted selection
	long maxCumulativeWeight;
	HashSet<Integer> protsNotToSample;
	
	int clusteringMeasure;
	double percentThreshold;

	public MotifSampling(String inputFile, ArrayList<Protein> protList, double[][] dm, int clustering_measure, double percent_threshold) {
		distanceMatrix = dm;
		proteinsInNetworkList = protList;
		
		clusteringMeasure = clustering_measure;
		percentThreshold = percent_threshold;
		
		proteinToOccurrenceMap = loadProteinOccurrenceList(inputFile);
		cumulativeSumOfWeights = computeCumulativeSumOfWeights();

		maxCumulativeWeight = cumulativeSumOfWeights[cumulativeSumOfWeights.length-1];
		protsNotToSample = listProteinsNotToSample();
	}

	/** 
	 * Load annotated proteins and their occurrence as a HashMap. 
	 * @param inputFile			Text File containing annotated proteins and their occurrence
	 * @return weightedList 	HashMap<String, Integer> mapping Annotated protein: occurrence of annotation
	 */
	private HashMap<String, Integer> loadProteinOccurrenceList(String inputFile){
		HashMap<String, Integer> weightedList = new HashMap<String, Integer>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine();
			while(line != null) {
				String[] col = line.split("\t");
				weightedList.put(col[0], Integer.parseInt(col[1])); // col[0] = protein name, col[1] = occurrence
				line = in.readLine();
			}
			in.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return weightedList;
	}

	/**
	 * Generate the cumulative sum of weights list in the same order as the proteins in network list,
	 * ensuring that the index of the weight corresponds to it's network proteins index for ease of look up. 
	 * 
	 * @return cumulativeWeightList	list of cumulative weights.
	 */
	private long[] computeCumulativeSumOfWeights() {
		long[] cumulativeWeightList = new long[this.proteinsInNetworkList.size()]; // initialize list the size of network list
		long cumulativeWeight = 0; // initialize cumulative weight
		int countMissingProts = 0;
		/* iterate all proteins in the order of the network protein list; 
		 * obtain weight of protein as it's annotation occurrence,
		 * update it's weight as a cumulative weight and set in list */
		System.out.println("Not annotated proteins:");
		for(int i=0; i<proteinsInNetworkList.size(); i++) {
			if(proteinToOccurrenceMap.containsKey(proteinsInNetworkList.get(i).getProteinName())) {
				long currentWeight = this.proteinToOccurrenceMap.get(this.proteinsInNetworkList.get(i).getProteinName());
				cumulativeWeight += currentWeight; 
				cumulativeWeightList[i] = cumulativeWeight;
			} else { 
				long currentWeight = 0;
				cumulativeWeight += currentWeight; 
				cumulativeWeightList[i] = cumulativeWeight;

				System.out.print(proteinsInNetworkList.get(i).getProteinName() + "_" + i + "|");
				countMissingProts++;
			}

		}
		System.out.println("\nNumber of missing proteins: " + countMissingProts + "\n");
		return cumulativeWeightList;
	}

	private HashSet<Integer> listProteinsNotToSample() {
		HashSet<Integer> proteinsNotToSampleList = new HashSet<>();

		for(int i=0; i<this.proteinsInNetworkList.size(); i++) {
			if(!this.proteinToOccurrenceMap.containsKey(this.proteinsInNetworkList.get(i).getProteinName())) {
				proteinsNotToSampleList.add(i);
			}

		}

		return proteinsNotToSampleList;
	}

	/**
	 * Computes distributions for multiple number of proteins in range from go_start to go_stop. 
	 * Every distribution is output in it's own file
	 *
	 * @param nProtToSampleLowerBound		beginning of the range of proteins that will be sampled
	 * @param nProtToSampleUpperBound		end of the range of proteins that will be sampled
	 * @param numOfTimesNetworkIsSampled 	number of times to calculate the distribution for each amount of proteins
	 * @param mcFilePrefix					String - file path prefix for output distribution
	 */
	public void computeMultipleDistributions(int nProtToSampleLowerBound, int nProtToSampleUpperBound, int numOfTimesNetworkIsSampled, String mcFilePrefix) {

		for (int n = nProtToSampleLowerBound; n <= nProtToSampleUpperBound; n++) { // range of proteins to sample

			System.out.println("Computing TPD: " + n);

			String mcFile = mcFilePrefix + "s" + numOfTimesNetworkIsSampled + "_n" + n;
			HashMap<Double, Double> distribution = computeMonteCarloDistribution(n, numOfTimesNetworkIsSampled);

			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(mcFile)));

				out.write("TPD (n = " + n + ")" + "\t" + "Frequency" + "\n");
				for (double dist : distribution.keySet()) {
					out.write(dist + "\t" + distribution.get(dist) + "\n");
					out.flush();

				} 
				out.write("\n");
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

		}

	}
	/**
	 * Computes the distribution of TPD for certain amount of randomly selected proteins.
	 *
	 * @param distanceMatrix 			distance matrix of proteins shortest pairwise distances
	 * @param timesToSampleNetwork   	number of iterations to sample Network
	 * @return the distribution of TPD for certain amount of randomly selected proteins
	 */
	private HashMap<Double, Double> computeMonteCarloDistribution(int numProteinsToSample, int timesToSampleNetwork) {

		HashMap<Double, Double> distribution;

		/* Sample network for x amounts of proteins, y amount of times */
		double[] tpdSampleList = sampleNetwork(numProteinsToSample, timesToSampleNetwork);

		/* Measure the frequencies of sampled total pairwise distances from x amount of proteins */
		distribution = Sampling.computeFrequenciesOfSampledTPDs(tpdSampleList, timesToSampleNetwork);

		/* Sanity check to ensure the sum of TPD frequencies equals 1 */
		Sampling.checkFrequencyTotal(distribution, numProteinsToSample);

		return distribution;
	}

	/**
	 * Sample network X amount of times for Y number of proteins. 
	 * 
	 * @param numProteinsToSample		Proteins to sample in network
	 * @param timesToSampleNetwork		Number of iterations to sample network
	 * @return	list of TPDs measured in network
	 */
	private double[] sampleNetwork(int numProteinsToSample, int timesToSampleNetwork) {
		double[] tpdSampleList = new double[timesToSampleNetwork]; // array to store the results of sampling {list of TPDs}

		/* Randomly select proteins from the network and compute their total pairwise distance
		 * as many times as the network is to be sampled (e.g. 10X7) */
		for (int i = 0; i < timesToSampleNetwork; i++) {

			if(i%10000 == 0) {
				System.out.print(i + ".");
			}

			if(i%100000 == 0) {
				System.out.println();
			}

			/* select proteins from the weighted list (ie. proteins are proportional to their occurrence in annotation list */
			ArrayList<Integer> randomProteins = getRandomWeightedProteinsWithBinarySearch(numProteinsToSample);

			/* compute the total pairwise distance from the proteins selected above */
			switch(clusteringMeasure) {
			case 0: tpdSampleList[i] = TopPercentPairwiseDistance.computeTPD(randomProteins, distanceMatrix);
			break;
			case 1: tpdSampleList[i] = TopPercentPairwiseDistance.getTPPD(randomProteins, distanceMatrix, percentThreshold);
			break;
			case 2: tpdSampleList[i] = TopPercentPairwiseDistance.getCoreTPD(randomProteins, distanceMatrix, percentThreshold);
			break;
			}

		}
		System.out.print("Done\n");
		return tpdSampleList;
	}

	/**
	 * Randomly selects proteins by identifying a random weight in the cumulative distribution
	 *
	 * @param numProteinsToSample number of proteins to select from network
	 * @return array of selected protein indexes
	 */
	@SuppressWarnings("unused")
	private ArrayList<Integer> getRandomWeightedProteins(int numProteinsToSample) {
		ArrayList<Integer> randomProteins = new ArrayList<>();

		/* Selection process occurs until the number of selected proteins equals number of proteins to sample from */ 
		while (randomProteins.size() < numProteinsToSample) {
			/* select a random weight between 0 and sum of cumulative weight */
			long selectedWeight = ThreadLocalRandom.current().nextLong(this.maxCumulativeWeight); // math random return number between 0 and 1 
			/* protein corresponding to selected random weight will be the first protein to be greater or equal to the weight */
			for(int protIndex=0; protIndex<this.cumulativeSumOfWeights.length; protIndex++) {
				if(selectedWeight < this.cumulativeSumOfWeights[protIndex]) {

					if(this.listProteinsNotToSample().contains(protIndex)) {
						System.out.println("ERROR: Sampled Protein " + protIndex);
					}

					if(!randomProteins.contains(protIndex)) {
						randomProteins.add(protIndex);
					}

					break; // when protein index is identified break out of the loop
				} 
			}
		}
		return randomProteins;
	}

	/**
	 * Randomly selects proteins by identifying a random weight in the cumulative distribution
	 * Binary search implemented as described: https://en.wikipedia.org/wiki/Binary_search_algorithm#Procedure_for_finding_the_leftmost_element
	 *
	 *
	 * @param numProteinsToSample number of proteins to select from network
	 * @return array of selected protein indexes
	 */
	private ArrayList<Integer> getRandomWeightedProteinsWithBinarySearch(int numProteinsToSample) {
		ArrayList<Integer> randomProteinsIdxList = new ArrayList<>();

		/* Selection process occurs until the number of selected proteins equals number of proteins to sample from */ 
		while (randomProteinsIdxList.size() < numProteinsToSample) {
			/* select a random weight between 0 and sum of cumulative weight */
			long selectedWeight = ThreadLocalRandom.current().nextLong(this.maxCumulativeWeight); // math random return number between 0 and 1 
			int l = 0;
			int r = this.cumulativeSumOfWeights.length;

			/* protein corresponding to selected random weight will be the first protein to be greater or equal to the weight */
			while(l < r) {
				int m = (int) Math.floor((l+r)/ (double)2);
				if(this.cumulativeSumOfWeights[m] < selectedWeight) {
					l = m + 1;
				} else { 
					r = m;
				}
			}

			if(!randomProteinsIdxList.contains(l)) {
				randomProteinsIdxList.add(l);
			}
		}
		return randomProteinsIdxList;
	}

}

