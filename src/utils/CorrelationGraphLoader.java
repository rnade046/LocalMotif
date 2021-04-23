package utils;
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


public class CorrelationGraphLoader {

	public static void loadGraphFromCorrelationNetwork(String inputRepository, String proteinsInNetworkFile) {
		
		for(double i=0.4; i<=0.8; i+=0.1) {

		HashMap<String, Double> confidentInteractionsMap = getConfidentInteractions(inputRepository, i);
		HashSet<String> confidentProteinSet = getConfidentProteins(confidentInteractionsMap);
		
		printProteinsInNetwork(confidentProteinSet, proteinsInNetworkFile);
		
		System.out.println("Threshold: " + i);
		System.out.println("#Confident interactions: " + confidentInteractionsMap.size());
		System.out.println("#Confident proteins: " + confidentProteinSet.size() + "\n");
		}
	}
	
	
	/**
	 * G
	 * @param inputRepository
	 * @param threshold
	 * @return
	 */
	private static HashMap<String, Double> getConfidentInteractions(String inputRepository, double threshold) {
		HashMap<String, Double> confidentInteractions = new HashMap<>();
		
		try {
			InputStream in = new FileInputStream(new File(inputRepository));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String line = input.readLine(); // header
			line = input.readLine();
			while(line!= null) {
				
				String[] col = line.split("\t");
				
				if(Double.parseDouble(col[2]) >= threshold) {
					String interactionOpt1 = col[0] + "_" + col[1];
					String interactionOpt2 = col[1] + "_" + col[0];
					
					if(!confidentInteractions.containsKey(interactionOpt1) || !confidentInteractions.containsKey(interactionOpt2)) {
						confidentInteractions.put(interactionOpt1, Double.parseDouble(col[2]));
					}
				}
			line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return confidentInteractions;
	}
	
	private static HashSet<String> getConfidentProteins(HashMap<String, Double> confidentInteractions){
		HashSet<String> confidentProteinSet = new HashSet<>();
		
		for(String interaction: confidentInteractions.keySet()) {
			
			String prot1 = interaction.split("\\_")[0];
			String prot2 = interaction.split("\\_")[1];
			
			confidentProteinSet.add(prot2);
			confidentProteinSet.add(prot1);
		}
		
		return confidentProteinSet;
	}
	
	private static void printProteinsInNetwork(HashSet<String> confidentProteinSet, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(String protein: confidentProteinSet) {
				out.write(protein + "\n");
				out.flush();
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
}
