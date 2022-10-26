import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map.Entry;

public class CompareProteinsToSequenceLength {

	public static void main(String[] args) {

		String degreeConnectionFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/network-topology/corrNetTop2_degreesInNetwork.tsv";
		String refSeqIdFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/corrNetTop2_proteinsInNetwork_info.tsv";
		String fastaFile = "/Users/rnadeau2/Documents/LESMoNlocal/input_files/human_3UTRsequences.txt";
		
		String outputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/network-topology/corrNetTop2-400_seqLengthvsDegree.tsv";
		
		assessProteinsToSequenceLength(degreeConnectionFile, refSeqIdFile, fastaFile, outputFile);
	}

	private static void assessProteinsToSequenceLength(String numberOfConnectionsFile, String refSeqIdFile, String fastaFile, String outputFile) {

		HashMap<String, Integer> degrees = loadDegreeConnections(numberOfConnectionsFile);

		HashMap<String, String> refSeqIds = loadRefSeqIds(refSeqIdFile); // id = protein

		HashMap<String, Integer> proteinLength = assessProteinLength(refSeqIds, fastaFile);	

		printResults(proteinLength, degrees, outputFile);

	}


	private static HashMap<String, Integer> loadDegreeConnections(String numberOfConnectionFile){

		HashMap<String, Integer> degrees = new HashMap<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(numberOfConnectionFile))));

			String line = in.readLine(); // header
			line = in.readLine();
			
			while(line!=null) {
				String[] col = line.split("\t");
				degrees.put(col[0], Integer.parseInt(col[1]));

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return degrees;
	}

	private static HashMap<String, String> loadRefSeqIds(String refSeqIdFile){

		HashMap<String, String> refSeqIds = new HashMap<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(refSeqIdFile))));

			String line = in.readLine(); // no header 
			while(line!=null) {
				String[] col = line.split("\t");
				
				if(col.length > 1 ) {
					String[] ids = col[1].split("\\|");

					for(String id : ids) {
						refSeqIds.put(id, col[0]);
					}	
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return refSeqIds;

	}
	
	private static HashMap<String, Integer> assessProteinLength(HashMap<String, String> refSeqIds, String fastaFile){
		
		HashMap<String, Integer> proteinLength = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));

			String line = in.readLine();
			boolean storeSeq = false;
			String seq = "";
			String id = "";

			while(line!=null) {

				/* store sequence (add multiple lines) */
				if(storeSeq) {
					seq += line;
				}

				if(line.startsWith(">")) {

					/* store length */
					if(!seq.isEmpty()) {
						String protein = refSeqIds.get(id);
						if(proteinLength.containsKey(protein)){
							int currentLength = proteinLength.get(protein);
							if(seq.length() > currentLength) {
								proteinLength.replace(protein, seq.length());
							}
						} else {
							proteinLength.put(protein, seq.length());
						}

					}

					storeSeq = false;
					seq = "";

					String[] col = line.split("_|\\.");
					id = col[2] + "_" + col[3]; // col[2] = type of ID (eg. NM) ; col[3] = number ID (#####)

					if(refSeqIds.containsKey(id)) {
						storeSeq = true;
					}
				}

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinLength;
	}
	
	private static void printResults(HashMap<String, Integer> proteinLengthMap, HashMap<String, Integer> degree, String outputFile) {
		
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			out.write("protein\tdegrees\tlength\n");
			for(Entry<String, Integer> prot: proteinLengthMap.entrySet()) {
				out.write(prot.getKey() + "\t" + degree.get(prot.getKey()) + "\t" + prot.getValue() + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

}
