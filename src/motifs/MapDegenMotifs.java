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
	
		/* Load map degen motifs = motifs */ 
		HashMap<String, String[]> motifToTestMap = loadDegenMotifsToTest(listDegenMotifsToMotifsTestFile);
		System.out.println("Number of loaded degen motifs: " + motifToTestMap.size());
		
		/* Load motifs to proteins */
		HashMap<String, String[]> motifToProteinMap = loadMotifsToProteinsMap(motifToProteinFile);
		System.out.println("Loaded motif to protein map: " + motifToProteinMap.size());
		
		System.out.println("Mapping annotations: ");
		/* Output mapping of degenMotifs to proteins */
		mapAnnotationFiles(degenMotifAnnotationFile, motifToTestMap, motifToProteinMap);
	}
	
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
				String[] proteinList = line.split("\t")[4].split("\\|");
				
				motifToProteinMap.put(motif, proteinList);
				
				line = input.readLine();
			}
			
			input.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return motifToProteinMap;
	}
	
	private static void mapAnnotationFiles(String outputFile, HashMap<String, String[]> degenMotifMap, HashMap<String, String[]> motifToProteinMap) {
		
		 try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			int count = 1;
			for(String degenMotif: degenMotifMap.keySet()) {
				
				if(count%1000 == 0) {
					System.out.print(count + ".");
				}
				
				if(count%10000 == 0) {
					System.out.println();
				}
				
				HashSet<String> proteinSet = new HashSet<>();
				
				for(String motif: degenMotifMap.get(degenMotif)) {
					
					for(String protein: motifToProteinMap.get(motif)) {
						 proteinSet.add(protein);
					}
					
				}
				if(!proteinSet.isEmpty()) {
					out.write(degenMotif + "\t" + proteinSet.size() + "\t");
					
					for(String protein: proteinSet) {
						out.write(protein + "|");
					}
					out.write("\n");
				}
				count++;
			}
			System.out.print("Done");
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
