import java.util.ArrayList;

import graph.Interaction;

import utils.GraphLoader;


public class Main {

	public static void main(String[] args) {

		String wd = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\";

		String cellMapFile = wd + "input_files\\saint-latest-HumanCellMap.txt";
		String preyMappingFile = wd + "input_files\\preysMappingFile.tsv";
		String baitMappingFile = wd + "input_files\\baitsMappingFile.txt";

		//String distanceMatrixFile = wd + "IO_files\\cellMap_distanceMatrix.txt";

		ArrayList<Interaction> interactionList =  GraphLoader.loadInteractionRepository(cellMapFile, baitMappingFile, preyMappingFile, 0.01);
	}

}
