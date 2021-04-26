package opt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

public class checkAnnotations {

	public static void main(String[] args) {

		String outputFile = "";
		
		int numFiles = 992;
		HashMap<Integer,Integer> mapOfProteinFreq = calculateProteinFreq(numFiles);
		outputFreqOfProteinMap(mapOfProteinFreq, outputFile);

	}

	private static HashMap<Integer, Integer> calculateProteinFreq(int numFiles) {

		HashMap<Integer, Integer> mapOfProteinFreq = new HashMap<>();

		for(int i=0; i < numFiles; i++) {

			String annotationFile = "mapDegenMotifsToProteins_" + i + ".tsv"; 

			InputStream in;
			try {
				in = new FileInputStream(new File(annotationFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));
				
				String line = input.readLine();
				
				while(line!=null) {
					
					Integer numOfProt = Integer.parseInt(line.split("\t")[2]);
					
					if(mapOfProteinFreq.containsKey(numOfProt)) {
						mapOfProteinFreq.put(numOfProt, mapOfProteinFreq.get(numOfProt) + 1);
					} else {
						mapOfProteinFreq.put(numOfProt, 1);
					}
					
					line = input.readLine();
				}
				
				input.close();
			} catch(IOException e) {
				e.printStackTrace();
			}
		}

		return mapOfProteinFreq;
	}
	
	private static void outputFreqOfProteinMap(HashMap<Integer, Integer> mapOfProteinFreq, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(int numProt: mapOfProteinFreq.keySet()) {
				out.write(numProt + "\t" + mapOfProteinFreq.get(numProt) + "\n");
				out.flush();
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

}
