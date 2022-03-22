package fdr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;

public class FdrCalculator {

	/* p-values associated to regular motifs */
	// private ArrayList<Double> motifsSignificanceScores;
	/* p-values associated to null model */
	//private ArrayList<Double> nullSignificanceScores;

	private String regMotifFile;
	private String nullMotifFile; 
	
	private int regMotifSize;
	private int nullMotifSize;
	
	private double minPval;

	/**
	 * Constructor of the class.
	 *
	 * @param goAnnotations list of goTerms
	 * @param shuffledGoAnnotations list of shuffled goTerms
	 */
	public FdrCalculator(String motifSignificanceScoresFile, String nullSignificanceScoresFile) {
		System.out.println("Loading regular motifs");
		//this.motifsSignificanceScores = loadSignificanceScores(motifSignificanceScoresFile);
		this.regMotifFile = motifSignificanceScoresFile;
		double[] regMotifInfo = identifyLowestPval(regMotifFile);
		System.out.println("Loading null motifs");
		this.nullMotifFile = nullSignificanceScoresFile;
		double[] nullMotifInfo = identifyLowestPval(nullMotifFile);
		//this.nullSignificanceScores = loadSignificanceScores(nullSignificanceScoresFile);

	   
		this.minPval = Math.min(regMotifInfo[0], nullMotifInfo[0]);
		
		this.regMotifSize = (int) regMotifInfo[1];
		this.nullMotifSize = (int) nullMotifInfo[1];
		
        System.out.println("Total tested motifs: " + regMotifSize);
        System.out.println("Total tested null motifs: " + nullMotifSize);
     
	}

	/**
	 * Computes the false discovery rate table.
	 *
	 * @return array of FalseDiscoveryRates for p-values from 1e-180 to 1.
	 */
	public ArrayList<FalseDiscoveryRate> computeFdr() {

		ArrayList<FalseDiscoveryRate> fdrs = new ArrayList<>();

		/* Computes the false discovery rate and #annotations that pass FDR for the certain p-value */
		for (double pVal = this.minPval; pVal < 0.1; pVal = pVal * 2) {
			System.out.println("pval = " + pVal);
			double[] info = computeFDR(pVal); // info [0] = FDR, [1] = annotation that pass threshold\

			fdrs.add(new FalseDiscoveryRate(info[0], pVal, (int) info[1], (int) info[2]));
		}

		for (double pVal = 0.1; pVal <= 1.01; pVal += 0.05) {
			System.out.println("pval = " + pVal);
			double[] info = computeFDR(pVal); // info [0] = FDR, [1] = annotation that pass threshold

			fdrs.add(new FalseDiscoveryRate(info[0], pVal, (int) info[1], (int) info[2]));
		}

		System.out.println("Performing monotonic transformation");
		fdrs = monotonicTransformationForFdr(fdrs); //monotonic transformation

		return fdrs;
	}

	@SuppressWarnings("unused")
	private ArrayList<Double> loadSignificanceScores(String inputFile){

		ArrayList<Double> significantScores = new ArrayList<>();

		try {
			InputStream in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			int count = 1;
			while(line!=null) {
				if(count%10000 == 0) {
					System.out.print(count + ".");
				}
				if(count%100000 ==0) {
					System.out.println();
				}
				significantScores.add(Double.parseDouble(line));
				line = input.readLine();
				count++;
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Sorting motifs");
		Collections.sort(significantScores); // sort list increasing order

		return significantScores;
	}

	/*    *//**
	 * Modifies every goTerm in both shuffled and original lists with TPD.
	 *
	 * @param distanceMatrix distance matrix which is used to calculate TPD.
	 *//*
    public void modifyGoAnnotationsWithTPD(double[][] distanceMatrix) {
    	System.out.println("computing annotation TPDs:");
        Modifier.setClusterTPD(goAnnotations, distanceMatrix); // Annotation Go List
        System.out.println("");
        System.out.println("computing shuffled annotation TPDs:");
        Modifier.setClusterTPD(shuffledGoAnnotations, distanceMatrix); // Shuffled Annotations
        System.out.println("");
    }

	  *//**
	  * Modifies every goTerm in both shuffled and original lists with a p-value. P-value is calculated using normal
	  * approximation (which is built according to previously computed parameters), so it is obligatory to have a file
	  * with normal distribution parameters before calling this function. If you don't have this file, you can call
	  * computeNormalDistributionParameters and then pass the name of the new file to this function.
	  *
	  * @param distributionParametersFilePath path to the file with distribution parameter
	  *//*
    public double modifyGoAnnotationsWithPvalueFromNormalApproximation(String distributionParametersFilePath, int numOfSampling) {

    	double[] minimum_pvals = new double[2];
        minimum_pvals[0] = NormalApproximation.importNormalDistributionParameters(goAnnotations, distributionParametersFilePath, numOfSampling);
        minimum_pvals[1] = NormalApproximation.importNormalDistributionParameters(shuffledGoAnnotations, distributionParametersFilePath, numOfSampling);

        double min_pval = Math.min(minimum_pvals[0], minimum_pvals[1]);

        return min_pval;
    }

	   *//**
	   * Computes all normal distributions parameters and stores them in the file.
	   *
	   * @param distributionFilePath 				path to the file with distributions
	   * @param nProtToSampleUpperBound 			largest amount of protein for which we will compute distribution parameters
	   * @param distributionParametersFilePath 	path to the file with output (computed parameters)
	   *//*
    public void computeNormalDistributionParameters(String distributionFilePath, int nProtToSampleLowerBound, int nProtToSampleUpperBound, String distributionParametersFilePath) {
        Loader.loadMonteCarloDistributions(distributionFilePath, nProtToSampleLowerBound, nProtToSampleUpperBound, distributionParametersFilePath);
    	//Loader.loadDistributions2(distributionFilePath, nProtToSampleLowerBound, nProtToSampleUpperBound, distributionParametersFilePath);
    }
	    */



	private double[] computeFDR(double pThreshold) {
		/***************************************************************************************************
		 * Compute the false discovery rate for goTerms surpassing a certain p-value threshold (pThreshold)
		 *
		 * FDR = nOfShuffledGo_pval that passed / nGo_pval that passed
		 ***************************************************************************************************/

		double[] info = new double[3];

		int motifCount = countMotifsThatPassSignificanceThreshold(regMotifFile, pThreshold);
		int nullCount = countMotifsThatPassSignificanceThreshold(nullMotifFile, pThreshold);

		
		
//		for (int i = 0; i < this.motifsSignificanceScores.size(); i++) {
//
//			/* Count goTerms that have a p-value smaller than threshold */
//			if (this.motifsSignificanceScores.get(i) <= pThreshold) {
//				motifCount++;
//			}
//		}
//
//		for(int i=0; i < this.nullSignificanceScores.size(); i++) {
//			/* Count shuffled goTerms that have a p-value smaller than the threshold*/
//			if (this.nullSignificanceScores.get(i) <= pThreshold) {
//				nullCount++;
//			}
//
//		}
//
//		double fdr = (nullCount/ (double) this.nullSignificanceScores.size()) / ((double) motifCount / (double) this.motifsSignificanceScores.size());
		double fdr = (nullCount / (double) this.nullMotifSize) / ((double) motifCount / (double) this.regMotifSize);
		
		info[0] = fdr;
		info[1] = motifCount;
		info[2] = nullCount;

		return info;
	}

	/*    private int sumOfAnnotationsThatPassFDRthreshold(double p_val) {

        int motifCount = 0;

        for (int i = 0; i < this.motifsSignificanceScores.size(); i++) {
            //System.out.println(goAnnotation.getPvalue());
            if (this.motifsSignificanceScores.get(i) <= p_val) {
                motifCount++;
            }
        }


        return motifCount;
    }*/

	/**
	 * Modifies the FDR list so that FDR increases monotonously.
	 *
	 * @param fdrs list of already computed FDR values.
	 * @return list of FDR values after transformation.
	 */
	private ArrayList<FalseDiscoveryRate> monotonicTransformationForFdr(ArrayList<FalseDiscoveryRate> fdrs){
		for(int i = fdrs.size()-1; i > 0; i--){
			if (fdrs.get(i-1).getFalseDiscoveryRate() > fdrs.get(i).getFalseDiscoveryRate()){
				fdrs.get(i-1).setFalseDiscoveryRate(fdrs.get(i).getFalseDiscoveryRate());
			}
		}
		return fdrs;
	}

	/***
	 * NEW METHODS : modifying FDR calc to limit elements loaded in memory
	 */

	private double[] identifyLowestPval(String pvalFile) {
		double minPval=Double.MAX_VALUE;
		int countLines = 0;
		try {
			InputStream in = new FileInputStream(new File(pvalFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			countLines++;
			while(line!=null) {

				double currentValue = Double.parseDouble(line);
				if(currentValue != 0) {
					if(currentValue < minPval) {
						minPval = currentValue;
					}
				} else {
					System.out.println("zero-value p-val @line: " + countLines);
				}
				
				line = input.readLine();
				countLines++;
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		
		double[] info = new double[] {minPval, countLines};

		return info;
	}

	private int countMotifsThatPassSignificanceThreshold(String inputFile, double threshold) {
		int count = 0;

		try {
			InputStream in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			while(line!=null) {

				double currentValue = Double.parseDouble(line);
				if(currentValue <= threshold) {
					count++;
				}
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}



		return count;
	}

}
