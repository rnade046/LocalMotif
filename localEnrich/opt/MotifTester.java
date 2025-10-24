package opt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MotifTester {

	public static void main(String[] args) {
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File("src/testMotifs/GLA_seq.txt")));

			String sequence = in.readLine();
			sequence = in.readLine();

			System.out.println(sequence);

			Set<String> formatedMotifs = new HashSet<String>();

			for(int i = 0; i < 1000; i++) {
				BufferedReader in2 = new BufferedReader(new FileReader(new File("src/CompanionFiles_localMotif_n20_2000/localMotif_n20_2000_motifAnnotationsCompanionFile_"+i)));

				String motifLine = in2.readLine();
				while(motifLine!=null) {
					if(motifLine.length()>1) {

						String motif = motifLine.split("\t")[0];
						//	if(Double.parseDouble(motifLine.split("\t")[1])<2000) {
						HashMap<Character, String> characterMap = new HashMap<>();
						characterMap.put('A', "A");
						characterMap.put('C', "C");
						characterMap.put('G', "G");
						characterMap.put('T', "T");
						characterMap.put('R', "[AG]");
						characterMap.put('Y', "[CT]");
						characterMap.put('D', "[ACG]");
						characterMap.put('B', "[ACT]");
						characterMap.put('H', "[AGT]");
						characterMap.put('V', "[CGT]");
						characterMap.put('*', ".");
						
						String formatedMotif = "";
						for(int j = 0; j < 8; j++){
							formatedMotif = formatedMotif + characterMap.get(motif.charAt(j));
							 
						}
						formatedMotifs.add(formatedMotif);
						//	}
					}
					motifLine = in2.readLine();
				}
				//System.out.println(i);
			}
			//String formatedMotif = formatMotifWithRegularExpression(motif);

			//System.out.println(searchSeqForMotif(formatedMotif, sequence));
			int counter=0;
			for(String motif: formatedMotifs){
				Pattern pattern = Pattern.compile(motif);
				Matcher matcher = pattern.matcher(sequence);
				boolean found = matcher.find();
				if(found){
					//	System.out.println(motif);
					counter++;
				}
			}
			System.out.println(counter);
		}catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private static HashSet<String> loadTestedMotifs(String motifPrefix){
		HashSet<String> motifs = new HashSet<>();
		
		for(int i=0; i<1000; i++) {
			InputStream in;
			try {
				in = new FileInputStream(new File(motifPrefix + i));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // no header
				while(line!=null) {

					motifs.add(line.split("\t")[0]);
					line = input.readLine();
				}
				input.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return motifs;
	}
}
