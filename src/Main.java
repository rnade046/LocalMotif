import java.io.File;
import java.util.ArrayList;

import graph.Interaction;
import graph.Protein;
import utils.Calculator;
import utils.DistanceMatrix;
import utils.GraphLoader;
import utils.Loader;
import utils.NetworkProteins;


public class Main {

	public static void main(String[] args) {

		String wd = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\";

		String cellMapFile = wd + "input_files\\saint-latest-HumanCellMap.txt";
		String preyMappingFile = wd + "input_files\\preysMappingFile.tsv";
		String baitMappingFile = wd + "input_files\\baitsMappingFile.txt";

		String distanceMatrixFile = wd + "IO_files\\cellMap_distanceMatrix.txt";

		ArrayList<Interaction> interactionList =  GraphLoader.loadInteractionRepository(cellMapFile, baitMappingFile, preyMappingFile, 0.01);
		System.out.println("Number of interactions: " + interactionList.size());
		
		ArrayList<Protein> proteinList = NetworkProteins.getProteinsInNetwork(interactionList);
		System.out.println("Number of Proteins: " + proteinList.size());
		
		
		/*File f = new File(distanceMatrixFile);
		if(!f.exists() && !f.isDirectory()) {
			DistanceMatrix.computeDistanceMatrix(interactionList, proteinList, distanceMatrixFile);
		}
		double[][] distanceMatrix = Loader.loadDistanceMatrix(distanceMatrixFile, proteinList); */

		
	}

}
