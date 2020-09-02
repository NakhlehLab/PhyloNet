package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.VariationalInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.VariationalModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.Prior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.VariationalVariable;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Test on real data.
 */
public class TestRealData {
    public static void main(String[] args) throws IOException {
        /*
         * load data
         */
        Map<String, String> nameTranslate = new HashMap<>();
        nameTranslate.put("human", "H");
        nameTranslate.put("chimp", "C");
        nameTranslate.put("gorilla", "G");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Data/hobolth/output.4.fasta");
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        String taxon = null;
        while ((line = br.readLine()) != null) {
            if (line.charAt(0) == '>') {
                taxon = line.substring(1);
            } else {
                omap.put(nameTranslate.get(taxon), line);
            }
        }
        Utils.DATA = new Alignment(omap);

        /*
         * initialize inference
         */
        ModelTree initModel = HCGModelBuilder.getHCGModelInit();
        Prior prior = new Prior();
        VariationalInference algo = new VariationalInference(initModel, prior);

        algo.run();

        VariationalModel posterior = algo.getVariationalPosterior();

        for (VariationalVariable var:posterior.getNodeHeightVariableList()) {
            System.out.println(var.getMean());
            System.out.println(var.getStandardDeviation());
        }
        for (VariationalVariable var:posterior.getPopSizeVariableList()) {
            System.out.println(var.getMean());
            System.out.println(var.getStandardDeviation());
        }
    }
}
