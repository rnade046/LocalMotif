package opt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashSet;

public class FindMotifsAssociatedToProtein {

	public static void main(String[] args) {

		String protein = args[0];
		String annotationPrefix = args[1];
		String output = args[2];

		HashSet<String> motifs = new HashSet<>();
		//for(int i=0; i<numFiles; i++) {
		for(int i=0; i <= 999; i++) {

			try {
				InputStream in = new FileInputStream(new File(annotationPrefix + i + ".tsv"));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				BufferedWriter out = new BufferedWriter(new FileWriter(new File(output)));

				String line = input.readLine(); // no header

				while(line!= null) {

					String[] col = line.split("\t");
					String motif = col[0];
					HashSet<String> proteins = new HashSet<>(Arrays.asList(col[2].split("\\|")));

					if(proteins.contains(protein)) {
						motifs.add(motif);
					}
					line = input.readLine();
				}
				input.close();
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(output)));
			
			for(String m : motifs) {
				out.write(m + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
