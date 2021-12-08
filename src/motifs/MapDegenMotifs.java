package motifs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;

public class MapDegenMotifs {

	public static void mapDegenMotifsToProteins(String listDegenMotifsToMotifsTestFile, String motifToProteinFile, String degenMotifAnnotationFile) {

		/* Load motifs to proteins */
		HashMap<String, String[]> motifToProteinMap = loadMotifsToProteinsMap(motifToProteinFile);
		System.out.println("Loaded motif to protein map: " + motifToProteinMap.size() + "\n");

		/* Load map degen motifs = motifs line by line */ 
		System.out.println("Loading and mapping annotations: ");
		InputStream in;
		try {
			in = new FileInputStream(new File(listDegenMotifsToMotifsTestFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); //no header
			int countMotifs = 0; 
			while(line != null) {
				
				if(countMotifs%1000 == 0) {
					System.out.print(countMotifs + ".");
				}
				
				if(countMotifs%10000 == 0) {
					System.out.println();
				}
				
				String degenMotif = line.split("\t")[0];
				String[] motifList = line.split("\t")[1].split("\\|");
				
				mapAnnotationFiles(degenMotifAnnotationFile, degenMotif, motifList, motifToProteinMap);

				line = input.readLine();
				countMotifs++;
			}
			System.out.println("Done");
			System.out.println("Total tested motifs: " + countMotifs);
			
			input.close();
		} catch(IOException e) {
			e.printStackTrace();
		}


		
		
	}

	@SuppressWarnings("unused")
	private static HashMap<String, String[]> loadDegenMotifsToTest(String degenMotifsToMotifsFile) {
		HashMap<String, String[]> motifMap = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(degenMotifsToMotifsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				String degenMotif = line.split("\t")[0];
				String[] motifList = line.split("\t")[1].split("\\|");

				motifMap.put(degenMotif, motifList);

				line = input.readLine();
			}

			input.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return motifMap;
	}

	private static HashMap<String, String[]> loadMotifsToProteinsMap(String motifToProteinFile){
		HashMap<String, String[]> motifToProteinMap = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(motifToProteinFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				String motif = line.split("\t")[0];
				// String[] proteinList = line.split("\t")[4].split("\\|");
				String[] proteinList = line.split("\t")[1].split("\\|");

				motifToProteinMap.put(motif, proteinList);

				line = input.readLine();
			}

			input.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return motifToProteinMap;
	}

	private static void mapAnnotationFiles(String outputFile, String degenMotif, String[] motifSet, HashMap<String, String[]> motifToProteinMap) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile), true));

			/* For every motif, identify list of annotating proteins */
			HashSet<String> proteinSet = new HashSet<>();
			for(String motif: motifSet) {
				
				if(motifToProteinMap.containsKey(motif)) {
					for(String protein: motifToProteinMap.get(motif)) {
						proteinSet.add(protein);
					}
				}
			}
			
			/* if at least one protein is annotated by the motif, the motif will be printed to annotation file */
			if(!proteinSet.isEmpty()) {
				out.write(degenMotif + "\t" + proteinSet.size() + "\t");

				for(String protein: proteinSet) {
					out.write(protein + "|");
				}
				out.write("\n");
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
