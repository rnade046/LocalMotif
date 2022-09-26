package motifInfo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;

public class Localization {

	public static void main(String[] args) {

	String inputLocalizations = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/subcellular_location.tsv";
	String summarizedLocalizatiosn = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/summarized_locations.tsv";
	
	determinePossibleLocalizations(inputLocalizations, summarizedLocalizatiosn);
		
	}


	public static void determinePossibleLocalizations(String inputLocalization, String summarizedLocalizations) {

		HashSet<String> localizations = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputLocalization))));

			String line = in.readLine(); // no-header

			while(line!=null) {

				String[] local = line.split("\t")[3].split("\\;"); // [3] = localization
				for(String l : local) {
					localizations.add(l);
				}
				line = in.readLine();
			}
			in.close();
			
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(summarizedLocalizations)));
			for(String local : localizations) {
				out.write(local + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	
}
