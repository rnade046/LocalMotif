package fdr;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;

public class ExportFdrDetails {

    public static void exportFDR(ArrayList<FalseDiscoveryRate> falseDiscoveryRates, String fdrExportFile) {

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(fdrExportFile));
            out.write("FDR" + "\t" + "Pval" + "\t" + "#Motifs" + "\t" + "#NullMotifs" + "\n");

            for (FalseDiscoveryRate fdr : falseDiscoveryRates) {
                out.write(fdr.getFalseDiscoveryRate() + "\t" + fdr.getPvalue() + "\t" + fdr.getPassingAnnotation() + "\t" + fdr.getNullAnnotations() + "\n");
                out.flush();
            }

            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

	
}
