import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class Test {

	public static void main(String[] args) throws InterruptedException, ExecutionException {
		Timestamp timestamp = new Timestamp(System.currentTimeMillis());
        System.out.println(timestamp + " [start]");
        System.out.println("Memory usage: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())+ " bytes");
        System.out.println();
		
		ArrayList<String> seqList = new ArrayList<>(Arrays.asList("AUUCGUAC", "UAUACAUA", "CGCCAAUU"));
		HashMap<Character, Character[]> iupac = defineIUPAC();
		
		HashMap<String, HashSet<String>> seqMap = generateDegenMotifs(seqList, iupac);

		timestamp = new Timestamp(System.currentTimeMillis());
        System.out.println(timestamp + " [end]");
        System.out.println("Memory usage: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())+ " bytes");
	}

	
	public static HashMap<Character, Character[]> defineIUPAC(){
		HashMap<Character, Character[]> iupac = new HashMap<>();
		
		Character[] a = {'A', 'R', 'D', 'H', 'V', '*'};
		Character[] u = {'U', 'Y', 'B', 'D', 'H', '*'};
		Character[] g = {'G', 'R', 'B', 'D', 'V', '*'};
		Character[] c = {'C', 'Y', 'B', 'H', 'V', '*'};
		
		iupac.put('A', a);
		iupac.put('U', u);
		iupac.put('G', g);
		iupac.put('C', c);
		
		return iupac;
	}
	
	public static HashMap<String, HashSet<String>> generateDegenMotifs(ArrayList<String> seqList, HashMap<Character, Character[]> iupacMap) throws InterruptedException, ExecutionException {
		
		HashMap<String, HashSet<String>> seqMap = new HashMap<>();
		
		List<Callable<HashMap<String, HashSet<String>>>> tasks = new ArrayList<Callable<HashMap<String, HashSet<String>>>>();
		
		for (final String seq: seqList) {
		    Callable<HashMap<String, HashSet<String>>> c = new Callable<HashMap<String, HashSet<String>>>() {
		        @Override
		        public HashMap<String, HashSet<String>> call() throws Exception {
		            return degenerate(seq, iupacMap);
		        }
		    };
		    tasks.add(c);
		}
		
		ExecutorService EXEC = Executors.newCachedThreadPool();
		try {
            long start = System.currentTimeMillis();
            
            List<Future<HashMap<String, HashSet<String>>>> results = EXEC.invokeAll(tasks);
            for (Future<HashMap<String, HashSet<String>>> fr : results) {
               seqMap.putAll(fr.get());
                
            }
            long elapsed = System.currentTimeMillis() - start;
            System.out.println(String.format("Eslapsed time: %d ms", elapsed));
		
		}  finally {
            EXEC.shutdown();
        }
		return seqMap;
	}
		  
	public static HashMap<String, HashSet<String>> degenerate(String sequence, HashMap<Character, Character[]> iupacMap) {
	    
		
		int solutions = 1;
	    for(int i = 0; i < sequence.length(); i++) {
	    	solutions *= iupacMap.get(sequence.charAt(i)).length;
	    }
	    HashSet<Character> degenCharacters = new HashSet<>(Arrays.asList('R', 'D', 'H', 'V', 'B', '*'));
	    HashMap<String, HashSet<String>> seqMap = new HashMap<String, HashSet<String>>();
	    //String[] degenerateMotifs = new String[solutions];
	    HashSet<String> degenerateMotifs = new HashSet<>();
	    for(int i = 0; i < solutions; i++) {
	        int j = 1;
	        String seq = "";
	        int degenCount = 0;
	        for(int k=0; k < sequence.length() ; k++) {
	        	Character[] set = iupacMap.get(sequence.charAt(k));
	        	
	        	char charToAppend = set[(i/j)%set.length];
	            seq += charToAppend;
	            j *= set.length;
	            
	           
	            if(degenCharacters.contains(charToAppend)) {
	            	degenCount++;
	            }
	        }
	        if(degenCount <= 7) {
	        	degenerateMotifs.add(seq);
	        }
	        
	    }
	    seqMap.put(sequence, degenerateMotifs);
	    return seqMap;
	}

}
