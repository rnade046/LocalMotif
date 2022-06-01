package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

public class IdLongestSequences {

	public static void main(String[] args) {

		String proteinInfoFile = "C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\corrNetTop2_proteinsInNetwork_info.tsv";
		String fastaFile = "C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\input_files\\human_3UTRsequences.txt";
		String outputFile = "C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\MotifPosition\\corrNetTop2_longestSequence.tsv";

		/* load protein name and associated reFseq id*/ 
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinInfoFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			int protCount = 1;
			String line = in.readLine();
			while(line!=null) {

				System.out.print(protCount + ".");
				if(protCount%50 == 0) {
					System.out.println();
				}

				if(line.split("\t").length > 1) {
					String protName = line.split("\t")[0];
					String[] refSeqIds = line.split("\t")[1].split("\\|");

					/* for each protein; go through FASTA and get sequence length, determine longest sequence; 
					 * print protein names + refSeq Id of longest sequence */
					if(refSeqIds.length>1) {

						String id = getIdForLongestSequence(refSeqIds, fastaFile);
						out.write(protName + "\t" + id + "\n");
						out.flush();
					} else { 
						out.write(protName + "\t" + refSeqIds[0] + "\n");
						out.flush();
					}
				}
				line = in.readLine();
				protCount++;
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static String getIdForLongestSequence(String[] refSeqIds, String fasta) {

		/* convert array to set */
		HashSet<String> ids = new HashSet<>(Arrays.asList(refSeqIds));

		/* initialize map {id = seqLength} */
		HashMap<String, Integer> mapLengthToIds = new HashMap<>();

		/* search for FASTA */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader (new FileInputStream(new File(fasta))));
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
						mapLengthToIds.put(id, seq.length());

						if(mapLengthToIds.size() == ids.size()) {
							break;
						}
					}

					storeSeq = false;
					seq = "";

					String[] col = line.split("_|\\.");
					id = col[2] + "_" + col[3]; // col[2] = type of ID (eg. NM) ; col[3] = number ID (#####)

					if(ids.contains(id)) {
						storeSeq = true;
					}
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		/* id longest sequence*/
		int maxValue = 0;
		String finalId = "";
		for(Entry<String, Integer> map : mapLengthToIds.entrySet()) {

			if(map.getValue() > maxValue) {
				finalId = map.getKey();
			}

		}
		return finalId;
	}
}
