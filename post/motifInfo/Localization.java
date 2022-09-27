package motifInfo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;

public class Localization {

	public static void main(String[] args) {

		String inputLocalizations = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/subcellular_location.tsv";
		String summarizedLocalizatiosn = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/summarized_locations.tsv";

		HashSet<String> local = determinePossibleLocalizations(inputLocalizations, summarizedLocalizatiosn);
		HashSet<String> motifs = loadMotifsToTest("/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/coreTPD_motifsOfInterest.tsv");

		for(String m :motifs) {
			String motifFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/motifsOfInterest/coreTPD_motifOfInterest_localization_"+ m + ".tsv"; 
			updateLocalizationCount(motifFile, local, motifFile);
		}

	}


	public static HashSet<String> determinePossibleLocalizations(String inputLocalization, String summarizedLocalizations) {

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
		return localizations;
	}


	private static HashSet<String> loadMotifsToTest(String motifsToTestFile){

		HashSet<String> motifsToTest = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(motifsToTestFile))));

			String line = in.readLine(); // no-header

			while(line!=null) {

				motifsToTest.add(line.split("\t")[1]);
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifsToTest;
	}
	
	private static void	updateLocalizationCount(String motifLocalizations, HashSet<String> possibleLocalizations, String updateCount) {
		
		HashMap<String, Double> localCount = new HashMap<>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(motifLocalizations))));

			String line = in.readLine(); 
			line = in.readLine();

			while(line!=null) {

				String[] col = line.split("\t");
				localCount.put(col[0], Double.parseDouble(col[1]));
				line = in.readLine();
			}

			in.close();
			
			
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(updateCount)));
			
			out.write("Localization\tCount\n");
			for(String local: possibleLocalizations) {
				if(localCount.containsKey(local)) {
					out.write(local + "\t" + localCount.get(local) + "\n");
				} else {
					out.write(local + "\t0\n");
				}
				out.flush();
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		
	
	
	}

}
