package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test.multi;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test.HCGModelBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.Prior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variationalMulti.VariationalInferenceMulti;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variationalMulti.VariationalModelMulti;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class TestVariationalMulti {
    public static void main(String[] args) throws IOException {
        testHCG();
    }

    private static void testHCG() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "H");
        msName2SpeciesName.put("2", "C");
        msName2SpeciesName.put("3", "G");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/2/aligned.fasta");
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        String taxon = null;
        while ((line = br.readLine()) != null) {
            if (line.charAt(0) == '>') {
                taxon = line.substring(1);
            } else {
                omap.put(msName2SpeciesName.get(taxon), line);
            }
        }
        Utils.DATA = new Alignment(omap);

        /*
         * initialize inference
         */
        ModelTree initModel = HCGModelBuilder.getHCGModelInit();
        Prior prior = new Prior();
        VariationalInferenceMulti algo = new VariationalInferenceMulti(initModel, prior);

        // log time used by algo.run()
        long startTime = System.currentTimeMillis();
        algo.run();
        long endTime = System.currentTimeMillis();
        long totalTimeMillis = endTime - startTime;

        VariationalModelMulti posterior = algo.getVariationalPosterior();
        System.out.println("Mu: " + Arrays.toString(posterior.getDistribution().getMuVector()));
        System.out.println("Sigma: " + Arrays.deepToString(posterior.getDistribution().getHalfSigmaMatrix()));

        System.out.println("");
        System.out.println("Total execution time: " + totalTimeMillis / 1000.0 + " s");
        System.out.println("Time used to build HMM: " + Utils.buildingTime / 1000.0 + " s");
        System.out.println("Time used in forward algorithm: " + Utils.likelihoodTime / 1000.0 + " s");
    }
}
