package randomSeq;

public class RandomMain {
	
	public static void main(String[] args) {
		
		String wd = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\";
		
		String fastaFile = wd + "input_files\\human_3UTRsequences.txt";
		String rnaIdListFile = wd + "IO_files\\refSeqRNAids-HumanCellMap.tsv";
		String randomFastaFile = wd + "IO_files\\random_humna_3UTRsequences.txt";
		
		RandomizeSequences.generateRandomizedFasta(rnaIdListFile, fastaFile, randomFastaFile);
	}


}
