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
import java.util.Map.Entry;
import java.util.Set;

public class MergeComparedAnnotationsDifferences {

	public static void main(String[] args) {

		String wd = args[0];
		String proteinChangesFile = wd + "annotationChanges/ProteinAnnotationChanges_";
		String mergeFile = wd + "merged_annotationChanges_originalXRegex.tsv";

		/* degree map */
		HashMap<String, Integer> degreeMap = loadMap(proteinChangesFile + "0");

		/* merge all annotation changes */ 
		HashMap<String, int[]> protMap = initializeProteinModsMap(degreeMap.keySet());
		protMap = updateProteinAnnotationCountMap(protMap, proteinChangesFile);

		printMergeFile(protMap, degreeMap, mergeFile);
	}

	private static HashMap<String, Integer> loadMap(String inputFile) {

		HashMap<String, Integer> map = new HashMap<>();
		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();
			while(line!=null) {

				String[] col = line.split("\t"); // [0] = protein Name, [1] = integer (#degrees or #motif)
				map.put(col[0], Integer.parseInt(col[1]));

				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return map;
	}

	private static HashMap<String, int[]> initializeProteinModsMap(Set<String> proteinNames){

		HashMap<String, int[]> proteinModifications = new HashMap<>();

		for(String p: proteinNames) {
			proteinModifications.put(p, new int[2]);
		}
		return proteinModifications;
	}

	private static HashMap<String, int[]> updateProteinAnnotationCountMap(HashMap<String, int[]> protMap, String inputFilePrefix) {

		for(int i=0; i<1000; i++) {
			try {
				BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFilePrefix + i))));

				String line = input.readLine(); // header
				line = input.readLine();
				while(line!=null) {

					// [0] = protein Name, [1] = integer (#degrees or #motif), [2] = lost in REGEX, [3] = added in REGEX
					String[] col = line.split("\t"); 

					/* modifications */
					int[] counts = protMap.get(col[0]);
					counts[0] += Integer.parseInt(col[2]); // count removed
					counts[1] += Integer.parseInt(col[3]); // count added

					protMap.put(col[0], counts);
					line = input.readLine();
				}
				input.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return protMap;
	}

	private static void printMergeFile(HashMap<String, int[]> protMap, HashMap<String, Integer> degreeMap, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("ProteinName\t#Degrees\t#ProteinsRemoved\t#ProteinsAdded\n");
			
			for(Entry<String, int[]> prot: protMap.entrySet()) {
				
				out.write(prot.getKey() + "\t" + degreeMap.get(prot.getKey()) + "\t" + prot.getValue()[0] + "\t" + prot.getValue()[1] +"\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
