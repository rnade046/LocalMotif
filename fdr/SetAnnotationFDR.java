import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import graph.Annotation;

public class SetAnnotationFDR {
    
	public static void modifyClustersFdr(ArrayList<Annotation> clusters, ArrayList<FalseDiscoveryRate> fdrs){
		Collections.sort(fdrs, Comparator.comparingDouble((a) -> a.getPvalue()));
		for(Annotation cluster:clusters){
			setClusterFdr(cluster, fdrs);
		}
	}

	private static void setClusterFdr(Annotation cluster, ArrayList<FalseDiscoveryRate> fdrs){

		int l = 0;
		int r = fdrs.size()-1;
		while (l+1 < r){
			int mid = l + (r-l)/2;
			if (fdrs.get(mid).getPvalue() == cluster.getPvalue()){
				cluster.setFDR(fdrs.get(mid).getFalseDiscoveryRate());
				return;
			}
			else if (fdrs.get(mid).getPvalue() < cluster.getPvalue()){
				l = mid;
			} else{
				r = mid;
			}
		}
		cluster.setFDR(fdrs.get(r).getFalseDiscoveryRate());
	}
	
}
