package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.List;

public class PositionConservationInBins {

	//private static String fasta = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/MotifPosition/corrNetTop2_CDS_longestSequences.txt";
	private static String fasta = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/MotifPosition/corrNetTop2_reverse-complement-sequences-CDS.txt";
	//private static String fasta = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/MotifPosition/LocalRandom_CodingSequences.txt";
	
	private static int motifLength = 8;
	private static int bins = 100;

	public static void main (String[] args) {

		String listOfMotifs = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv";
		String outputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/MotifPosition/coreTPD0.4/percentile/corrNet2-400_coreTPD_p0.4_CDS_RevC_motifPositionsByBins_b" + bins + "_";

		getMotifPositionsFromLongestSequences(listOfMotifs, outputFile);

	}

	private static void getMotifPositionsFromLongestSequences(String representativeMotifsFile, String motifOutputPrefixFile) {


		/* Load motifs to test */
		System.out.println("Loading representative motifs");
		List<String> motifs = SequenceUtils.loadRepresentativeMotifs(representativeMotifsFile);

		for(int i=0; i < motifs.size(); i++) {

			System.out.println("Motif : " + (i+1) );
			/* Determine possible instance motifs */
			HashSet<String> possibleInstances = SequenceUtils.getPossibleMotifInstances(motifs.get(i));

			/* initialize motif positions list and considered sequences for normalization */
			int[] motifPositions = new int[bins];
			int[] consideredSequence = new int[bins];
			double[] normalizedPositions = new double[bins];


			/* initialize list to contain nucleotide frequencies */
			//List<List<double[]>> allNucleotideFrequencies = new ArrayList<>();
			//			List<HashSet<String>> proteinsInBin = new ArrayList<>();
			//			for(int j=0; j<1000; j=j+125) {
			//				proteinsInBin.add(new HashSet<String>());
			//			}
			//String proteinName = "";

			/* Check for motif in all FASTA of the filtered FASTA file */
			InputStream in;
			try {
				in = new FileInputStream(new File(fasta));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));
				int consideredSeqs = 0;
				int totalSeqs = 0;
				String line = input.readLine();
				while(line != null) {


					//String[] header = line.split("\\_");
					//proteinName = header[header.length-2] + "_" + header[header.length-1];					}

					/* load sequence; ignore lines that are identifiers */
					if(!line.startsWith(">")) {
						totalSeqs++;

						if(Math.floor(line.length() / bins) >= motifLength) {
							consideredSeqs++; 

							/* format sequence into substrings with bins of almost equal size */						
							String[] sequence = formatSequenceInBins(line); 

							/* search for motif positions */
							motifPositions = searchForMotifPositions(sequence, motifPositions, possibleInstances);

							/* update considered sequences*/
							consideredSequence = updateConsideredSequencesForNormalization(sequence, consideredSequence);

							/* perform nucleotide frequency calculations */
							/* search sequence in blocks of 125 BP; ignore the filler part */
							normalizedPositions = normalizePositions(motifPositions, consideredSequence, normalizedPositions);

							/* store sequence nucleotide frequency */
							//allNucleotideFrequencies.add(sequenceFreq);

							//determineProteinsInBins(sequence, proteinsInBin, possibleInstances, proteinName);
						}
					}
					line = input.readLine();
				}
				input.close();
				System.out.println("considered sequences = " + consideredSeqs + " / " + totalSeqs);
				/* Print positions */
				//String motifOutputFile = motifOutputPrefixFile+ "allSeq_motif" + (i+1);
				String motifNormalizedOutputFile = motifOutputPrefixFile+ "allSeq_Normalized_motif" + (i+1);
				//printMotifPosition(motifPositions, motifOutputFile);
				printNormalizedMotifPosition(normalizedPositions, motifNormalizedOutputFile);

				//double[] normalizedFreqs = normalizeNucleotideFrequencies(motifPositions, consideredSequences);
				//printProteinsPerBin(proteinsInBin, normalizedFreqs, protPerBinFile + "_motif" + (i+1) + ".tsv");

				/* calculate nucleotide frequency for each bin */
				//List<double[]> overallBinFrequency = calculateOverallBinNucleotideFrequency(allNucleotideFrequencies);

				//double[] overallNucleotideFrequency = calculateOverallNucleotideFrequency(overallBinFrequency);

				//printNucleotideFrequency(overallNucleotideFrequency, overallBinFrequency, nucleotideFreqFile + "_motif" + (i+1) + ".tsv");


			} catch (IOException e) {
				e.printStackTrace();
			}

		}
		System.out.println("Done");

	}

	private static String[] formatSequenceInBins(String sequence) {

		String[] sequenceBins = new String[bins];

		int[] binSizes = new int[bins];
		int[] binSizesReodered = new int[bins];

		int minimumBinSize = (int) Math.floor(sequence.length() / (double) bins); 
		int binRemainder = sequence.length() % bins; 

		/* determine size of bins (adding remainder to minimum bin sizes) */
		for(int i=0; i<bins; i++) {

			// determining indexes for unequal bins : https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length
			int binStart = i*minimumBinSize + Math.min(i, binRemainder);
			int binEnd = (i+1) * minimumBinSize + Math.min(i+1, binRemainder);
			int binRange = binEnd - binStart;

			binSizes[i] = binRange;
		}

		/* re-ordering bin sizes to pad edges */
		int countEven = 0;
		int countOdd = bins-1;
		for(int i=0; i<bins; i++) {
			if(i%2 ==0) {
				binSizesReodered[countEven] = binSizes[i];
				countEven++;
			} else { 
				binSizesReodered[countOdd] = binSizes[i];
				countOdd--;
			}
		}

		/* get sequence substring */
		int currentIdx = 0;
		for(int i=0; i<bins; i++) {
			
			int start = currentIdx;
			int end = start + binSizesReodered[i];
			currentIdx += binSizesReodered[i];
			
			sequenceBins[i] = sequence.substring(start, end);
		}
		return sequenceBins;
	}

	private static int[] searchForMotifPositions(String[] sequence, int[] motifPositions, HashSet<String> motifPossibilities) { 

		/* iterate over each bin */
		for(int i=0; i<sequence.length; i++) {

			String subsequence = sequence[i];

			/* search for motif in subsequence */
			for(int j=0; j<subsequence.length()-motifLength; j++) {

				/* update count when motif is found */
				if(motifPossibilities.contains(subsequence.substring(j, j+8))) {
					motifPositions[i] += 1;
				}
			}
		}
		return motifPositions;
	}

	private static int[] updateConsideredSequencesForNormalization(String[] sequence, int[] consideredSequences){

		/* iterate over each bin */
		for(int i=0; i<sequence.length; i++) {

			/* determine number of motifs tested for each subsequence (i.e. bins)  */
			int currentLength = sequence[i].length() - motifLength;
			consideredSequences[i] += currentLength;
		}
		return consideredSequences;
	}

	private static double[] normalizePositions(int[] positionsPerBin, int[] consideredSeq, double[] normalizedPositions) {

		for(int i=0; i<normalizedPositions.length; i++) {

			normalizedPositions[i] = positionsPerBin[i] / (double) consideredSeq[i];
		}
		return normalizedPositions;
	}

	private static void printNormalizedMotifPosition(double[] normalizedPositions, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=0; i<normalizedPositions.length; i++) {
				out.write((i+1) + "\t" + normalizedPositions[i] + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
