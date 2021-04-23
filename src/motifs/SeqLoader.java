package motifs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;

public class SeqLoader {
	
	/**
	 * Load the list of refseq Ids to test.  
	 * 
	 * @param inputFile		String - file path to list of refseq Ids of corresponding proteins in the protein localization network
	 * @return refSeqIDs	HashSet<String> - set of RefSeq Ids 
	 */
	public static HashSet<String> loadRefSeqIDsToTest(String inputFile){
		HashSet<String> refSeqIDs = new HashSet<String>();
		
		FileInputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String line = input.readLine(); // header
			line = input.readLine();
			while(line != null) {
				refSeqIDs.add(line.split("\t")[1]); // every line is a refSeq ID
				line = input.readLine();
			}
			
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	
		return refSeqIDs;
	}
	
	
	
}
