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

public class CompareNullModelSignificance {

	public static void main(String[] args) {
		
		String motifClustering_v1File = args[0];
		String motifClustering_v2File = args[1];
		String outputFile = args[2];
		
		/* Load motif clustering file for V1 -- map {motif = p-value} */
		HashMap<String, String> motifToPvalV1map = loadMap(motifClustering_v1File);
		
		/* Iterate through corresponding motif clustering file of V2*/
		InputStream in;
		try {
			in = new FileInputStream(new File(motifClustering_v2File));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			String line = input.readLine(); // no header
			while(line!=null) {
				
				String[] col = line.split("\t"); // [0] = protein name, [3] = p-value
				
				/* If find matching motifs between V1 and V2; output protein | p-val(V1) | p-val (V2) */
				if(motifToPvalV1map.containsKey(col[0])) {
					out.write(col[0] + "\t" + motifToPvalV1map.get(col[0]) + "\t" + col[3] + "\n");
					out.flush();
				}
				line = input.readLine();
			}

			out.close();
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	
	}
	
	private static HashMap<String, String> loadMap(String inputFile) {

		HashMap<String, String> map = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			while(line!=null) {

				String[] col = line.split("\t"); // [0] = protein Name, [3] = p-value
				map.put(col[0], col[3]);

				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return map;
	}

}
