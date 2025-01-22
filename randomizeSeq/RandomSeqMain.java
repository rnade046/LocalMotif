public class RandomSeqMain {
	
	public static void main(String[] args) {
		
		/* .fasta containing sequences to randomize */
		String fastaFile = args[0];
		
		/* mapping file {proteinName = refSeqId1|refSeqId2|...|} */
		String rnaIdListFile = args[1];
		
		/* output file */
		String randomFastaFile = args[0].substring(0, args[0].lastIndexOf('.')) + "_winShuffled.fasta";
		
		System.out.println("**Generating randomized fasta sequences**");
		RandomizeSequences.generateRandomizedFasta(rnaIdListFile, fastaFile, randomFastaFile);
	}


}
