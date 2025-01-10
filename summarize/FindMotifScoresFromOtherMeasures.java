import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map.Entry;

public class FindMotifScoresFromOtherMeasures {

	public static void main(String[] args) {
		
		/* input: 
		 *   + list of motifs to search for 
		 *   + directory to search (number of files to search 0-998 */
		String wd = args[0];
		String motifListFile = wd + "corrNetTop2-400_coreTPD_p0.4_signficantMotifs_p6.94181641095822E-7.tsv";
		String searchFilePrefix = wd + args[1];
		
		String outputFile = wd + "motif_scores_xCoreTPD0.4_n1105.tsv";
		
		/* load motifs to search */
		HashMap<String, Integer> motifs = loadMotifs(motifListFile);
		
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			/* key = motif (eg. VACCTC*B), value = file # where this value is found */
			for(Entry<String,Integer> e : motifs.entrySet()) {
				
				/* search corresponding file for motif details*/
				String details = findMotifDetails(e.getKey(), searchFilePrefix + e.getValue());
				
				/* check contents */
				if(details.isEmpty()) {
					System.out.println("did not find motif : " + e.getKey() + ", in file: " + e.getValue());
				} else { 
					out.write(details + "\n");
					out.flush();
				}
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		
	}

	
	private static HashMap<String, Integer> loadMotifs(String motifsFile){

		HashMap<String, Integer> motifMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(motifsFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {

				String[] col = line.split("\t");
				motifMap.put(col[0], Integer.parseInt(col[4]));

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifMap;
	}
	
	private static String findMotifDetails(String motif, String file) {
		
		String details = "";
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(file))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {

				String[] col = line.split("\t");
				if(col[0].equals(motif)) {
					details = line;
					break;
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return details;
	}


}
