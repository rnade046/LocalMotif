package motifComparaison;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class CISBPFormatter {

	public static void main(String[] args) {

		String db = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies\\Homo_sapiens_2022_02_28_12 53_pm\\pwms_all_motifs\\";
		String tomtomDB = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies\\CISBP_TomTomformatted.txt";

		printFormattedDBForTomTom(db, tomtomDB);

	}

	private static void printFormattedDBForTomTom(String databaseDir, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("MEME version 4\n\n" + "ALPHABET= ACGT\n\n" + 
					"Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"); // header

			File f = new File(databaseDir);
			String[] pathnames = f.list();
			int count = 1; 
			for(String path : pathnames) {

				/* Load the ppm (position probability matrix) for every family*/
				String fileX = databaseDir + path;
				ArrayList<ArrayList<Double>> ppm = loadPPM(fileX);

				if(!ppm.isEmpty()) {
					out.write("MOTIF"  + " Motif" + count + " " + path +"\n");
					out.write("letter-probability matrix: alength= 4 w= " + ppm.get(0).size() +"\n");

					for(int k=0; k<ppm.size(); k++) {
						for(int j=0; j<ppm.get(k).size(); j++) {
							out.write(ppm.get(k).get(j) + " ");
						}
						out.write("\n");
						out.flush();
					}
					out.write("\n");
					count++;
				}
			}
			out.write("#END\n"); // end signal

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static ArrayList<ArrayList<Double>> loadPPM(String inputFile) {

		ArrayList<ArrayList<Double>> matrix = new ArrayList<>(); 

		try {
			InputStream in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();

			while(line != null) {
				String[] col = line.split("\\s++");
				ArrayList<Double> nuclFreq = new ArrayList<>();

				for(int i=1; i<col.length; i++) {  // ignore first element

					if(!col[i].isEmpty()) {
						double value = Double.parseDouble(col[i]);
						nuclFreq.add(value);
					}
				}

				matrix.add(nuclFreq);

				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return matrix;
	}

}
