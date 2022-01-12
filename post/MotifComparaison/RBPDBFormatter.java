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

public class RBPDBFormatter {

	public static void main(String[] args) {

		String db = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies\\RBPDB_human_PFMDir\\";
		String outputFile = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies\\RBPDB_formatted.txt"; 
		String tomtomDB = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies\\RBPDB_TomTomformatted.txt";
		
		printFormattedDBforMotifComp(db, outputFile);
		printFormattedDBForTomTom(db, tomtomDB);
	}

	private static void printFormattedDBforMotifComp(String databaseDir, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("#INCLUSive Motif Model\n#\n"); // header

			File f = new File(databaseDir);
			String[] pathnames = f.list();
			for(String path : pathnames) {

				/* Load the ppm (position probability matrix) for every family*/
				String fileX = databaseDir + path;
				ArrayList<ArrayList<Double>> ppm = loadPPM(fileX);

				out.write("#ID = " + path + "\n");
				out.write("#W = " + ppm.get(0).size() + "\n");
				

				
				for(int k=0; k<ppm.get(0).size(); k++) {
					for(int j=0; j<ppm.size(); j++) {
						out.write(ppm.get(j).get(k) + "\t");
					}
					out.write("\n");
					out.flush();
				}
				out.write("\n");
			}
			out.write("\n"); // end signal

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	private static void printFormattedDBForTomTom(String databaseDir, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("MEME version 4\n\n" + "ALPHABET= ACGT\n\n" + 
			"Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"); // header

			File f = new File(databaseDir);
			String[] pathnames = f.list();
			for(String path : pathnames) {

				/* Load the ppm (position probability matrix) for every family*/
				String fileX = databaseDir + path;
				ArrayList<ArrayList<Double>> ppm = loadPPM(fileX);

				out.write("MOTIF" + path + "\n");
				out.write("letter-probability matrix: alength= 4 w= " + ppm.get(0).size() +"\n");

				for(int k=0; k<ppm.get(0).size(); k++) {
					for(int j=0; j<ppm.size(); j++) {
						out.write(ppm.get(j).get(k) + "\t");
					}
					out.write("\n");
					out.flush();
				}
				out.write("\n");
			}
			out.write("#END\n"); // end signal

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

	private static ArrayList<ArrayList<Double>> loadPPM(String inputFile) {

		ArrayList<ArrayList<Double>> ppm = new ArrayList<>(); 
		ArrayList<ArrayList<Double>> matrix = new ArrayList<>(); 
		boolean convertFreq = false;

		try {
			InputStream in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				String[] col = line.split("\\s++");
				ArrayList<Double> nuclFreq = new ArrayList<>();

				for(int i=0; i<col.length; i++) {

					if(!col[i].isEmpty()) {
						double value = Double.parseDouble(col[i]);
						nuclFreq.add(value);

						if(!convertFreq && value > 1) {
							convertFreq = true;
						}
					}
				}

				matrix.add(nuclFreq);

				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		if(convertFreq) {

			/* Determine total elements */
			int maxInstances = 0;
			for(int i=0; i<matrix.size(); i++) {
				maxInstances += matrix.get(i).get(0);
			}

			/* normalize matrix*/
			for(int i=0; i<matrix.size(); i++) {
				ArrayList<Double> nuclFreq = new ArrayList<>();
				for(int j=0; j<matrix.get(i).size(); j++) {
					nuclFreq.add(matrix.get(i).get(j) / (double) maxInstances);
				}
				ppm.add(nuclFreq);
			}

		} else {
			ppm.addAll(matrix);
		}


		return ppm;
	}
}
