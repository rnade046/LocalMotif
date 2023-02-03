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

public class CompareNullModels {

	public static void main(String[] args) {

		String wd = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\compare-nullModels\\";
		String networkDegrees = wd + "corrNetTop2-400_degreesInNetwork.tsv";
		String testedMotifs_v1 = wd + "corrNetTop2-400_nullModel_v1_protFreqAnnotation.tsv";
		String testedMotifs_v2 = wd + "corrNetTop2-400_nullModel_v2_protFreqAnnotation.tsv";
		String output = wd + "corrNet2-400_motifCountDifferenceByDegree.tsv";

		/* Load degrees network -- map{protein = # degrees} */
		HashMap<String, Integer> degreesByProteinMap = loadMap(networkDegrees);

		/* Load tested motifs in each null model per protein -- map {protein = # motifs} */
		HashMap<String, Integer> motifCountByProteinInV1Map = loadMapNoHeader(testedMotifs_v1);
		HashMap<String, Integer> motifCountByProteinInV2Map = loadMapNoHeader(testedMotifs_v2);

		/* Output info as protein | # degrees | [motifs_v1 - motifs_v2] */
		outputMotifByDegreeInfo(degreesByProteinMap, motifCountByProteinInV1Map, motifCountByProteinInV2Map, output);
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

	private static HashMap<String, Integer> loadMapNoHeader(String inputFile) {

		HashMap<String, Integer> map = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
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

	private static void outputMotifByDegreeInfo(HashMap<String, Integer> degreesByProteinInNetworkMap, HashMap<String, Integer> motifsByProteinInV1Map, HashMap<String, Integer> motifsByProteinInV2Map, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("Protein\tDegree\tDiff\n");

			for(Entry<String, Integer> mapping: degreesByProteinInNetworkMap.entrySet()) {

				String prot = mapping.getKey();

				/* check exists because some proteins in the network are not annotated by any motifs */
				if(motifsByProteinInV1Map.containsKey(prot) && motifsByProteinInV2Map.containsKey(prot)) {
					int motifCount1 = motifsByProteinInV1Map.get(prot);
					int motifCount2 = motifsByProteinInV2Map.get(prot);

					int difference = motifCount1 - motifCount2;

					out.write(prot + "\t" + mapping.getValue() + "\t" + difference + "\n");
					out.flush();
				}
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
