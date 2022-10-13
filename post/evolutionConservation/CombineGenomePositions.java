package evolutionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

public class CombineGenomePositions {

	public static void main(String[] args) {

		String positionFilePrefix = "/home/rnade046/projects/rrg-mlaval/rnade046/genomeBrowser/genomePositions/positions/genomePositions_";
		String genomePositionsFile = "allGenomePositions.tsv";
		String bedFile = "genomePositionsOfTestedMotifs.bed";

		combineGenomePositions(positionFilePrefix, genomePositionsFile, bedFile);
	}

	public static void combineGenomePositions(String positionFilePrefix, String allGenomePositionsFile, String bedFile) {

		/* determine chromosomes to check */
		List<String> chromosomes = listChromosomes();
		System.out.println("searching chromosomes");
		
		/* for each chromosomes; load known positions and combine them */
		for(String chr : chromosomes) {
			System.out.println(chr);
			
			HashSet<Integer> positions = new HashSet<>();

			/* load positions from the 1000 position files */
			positions = loadPositions(chr, positionFilePrefix, positions);

			/* combine positions and print*/ 
			printGenomePositions(chr, positions, allGenomePositionsFile);
		}

		/* convert genome positions to bed file */
		System.out.println("printing bed file");
		printBedFile(allGenomePositionsFile, bedFile);
	}

	private static List<String> listChromosomes(){
		List<String> chromosomes = new ArrayList<>();

		for(int i=1; i<=22; i++) {
			chromosomes.add(Integer.toString(i));
		}

		chromosomes.add("X");
		chromosomes.add("Y");

		return chromosomes;
	}

	private static HashSet<Integer> loadPositions(String chr, String positionFilePrefix, HashSet<Integer> positions){

		for(int i=0; i<=999; i++) {
			if(i%50 == 0) {
				System.out.println();
			}
			
			File f = new File(positionFilePrefix + i);
			if(f.exists() && !f.isDirectory()) { 
				try {
					BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(positionFilePrefix + i))));

					System.out.print(i +".");
					String line =in.readLine();
					while(line!=null) {

						/* assess position if lines contributes to current CHR */ 
						if(line.split("\t")[0].equals("chr" + chr + ":")) {
							String[] pos = line.split("\t")[1].split("\\,");
							for(String p : pos) {

								int startIdx = Integer.parseInt(p.split("-")[0]);
								int endIdx = Integer.parseInt(p.split("-")[1]);

								for(int j=startIdx; j<= endIdx; j++) {
									positions.add(j);
								}

							}
							break;
						}
						line =in.readLine();
					}
					in.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		return positions;
	}

	private static void printGenomePositions(String chr, HashSet<Integer> positions, String outputFile) {

		/* print and sort positions; 1 line per chromosomes*/
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(outputFile), true));
			System.out.println("printing");

			out.write("chr" + chr + "\t");

			List<Integer> order = new ArrayList<Integer>(positions);
			Collections.sort(order);

			System.out.println(chr + "\t" + order.size());

			int startIndex = order.get(0); // start index 
			int idx1 = startIndex;
			int idx2 = 0; 
			for(int i=1; i<order.size(); i++) {

				idx2 = order.get(i);

				if(idx2 != (idx1+1)) {
					out.write(startIndex + "-" + idx1 + ",");
					out.flush();

					startIndex = idx2;
				}

				idx1 = idx2;
			}



			out.write(startIndex + "-" + idx1 + "\n");
			out.flush();

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void printBedFile(String allGenomePositions, String bedFile) { 

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(allGenomePositions))));	
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(bedFile)));

			// bed -header
			out.write("Track - genome positions of tested motifs\n");

			String line = in.readLine();
			while(line != null) {
				String chr = line.split("\t")[0];
				String[] positions = line.split("\t")[1].split("\\,");

				for(String p : positions) {

					int startPosition = Integer.parseInt(p.split("-")[0]);
					int endPosition = Integer.parseInt(p.split("-")[1]);

					out.write(chr + "\t" + startPosition + "\t" + (endPosition+1) + "\n");
					out.flush();
				}
				line = in.readLine();
			}

			out.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
