package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test;


import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.VariationalInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.VariationalModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.Prior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.VariationalVariable;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Test the variational inference procedure.
 * Created by Xinhao Liu on 3/26/20.
 */
public class TestABCD {
    public static void main(String[] args) throws IOException {
        testABCD();
    }

    private static void testDaC_ACD() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "A");
        msName2SpeciesName.put("2", "B");
        msName2SpeciesName.put("3", "C");
        msName2SpeciesName.put("4", "D");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/ABCD_200000/aligned.fasta");
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
        omap.remove("B");
        Utils.DATA = new Alignment(omap);
        /*
         * initialize inference
         */
        ModelTree initModel = DaCABCDModelBuilder.getACDModelInit();
        Prior prior = new Prior();
        VariationalInference algo = new VariationalInference(initModel, prior);

        // log time used by algo.run()
        long startTime = System.currentTimeMillis();
        algo.run();
        long endTime = System.currentTimeMillis();
        long totalTimeMillis = endTime - startTime;

        VariationalModel posterior = algo.getVariationalPosterior();
        for (VariationalVariable var:posterior.getNodeHeightVariableList()) {
            System.out.println(var.getMean());
            System.out.println(var.getStandardDeviation());
        }
        for (VariationalVariable var:posterior.getPopSizeVariableList()) {
            System.out.println(var.getMean());
            System.out.println(var.getStandardDeviation());
        }

        System.out.println("");
        System.out.println("Total execution time: " + totalTimeMillis / 1000.0 + " s");
        System.out.println("Time used to build HMM: " + Utils.buildingTime / 1000.0 + " s");
        System.out.println("Time used in forward algorithm: " + Utils.likelihoodTime / 1000.0 + " s");
    }

    private static void testABCD() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "A");
        msName2SpeciesName.put("2", "B");
        msName2SpeciesName.put("3", "C");
        msName2SpeciesName.put("4", "D");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/ABCD_200000/aligned.fasta");
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
        ModelTree initModel = ABCDModelBuilder.getABCDModelInit();
        Prior prior = new Prior();
        VariationalInference algo = new VariationalInference(initModel, prior);

        long startTime = System.currentTimeMillis();
        algo.run();
        long endTime = System.currentTimeMillis();
        long totalTimeMillis = endTime - startTime;

        VariationalModel posterior = algo.getVariationalPosterior();
        for (VariationalVariable var:posterior.getNodeHeightVariableList()) {
            System.out.println(var.getMean());
            System.out.println(var.getStandardDeviation());
        }
        for (VariationalVariable var:posterior.getPopSizeVariableList()) {
            System.out.println(var.getMean());
            System.out.println(var.getStandardDeviation());
        }

        System.out.println("");
        System.out.println("Total execution time: " + totalTimeMillis / 1000.0 + " s");
        System.out.println("Time used to build HMM: " + Utils.buildingTime / 1000.0 + " s");
        System.out.println("Time used in forward algorithm: " + Utils.likelihoodTime / 1000.0 + " s");
    }

    private static void testDaC_ABC() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "A");
        msName2SpeciesName.put("2", "B");
        msName2SpeciesName.put("3", "C");
        msName2SpeciesName.put("4", "D");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/ABCD_200000/aligned.fasta");
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
        omap.remove("D");
        Utils.DATA = new Alignment(omap);
        /*
         * initialize inference
         */
        ModelTree initModel = DaCABCDModelBuilder.getABCModelInit();
        Prior prior = new Prior();
        VariationalInference algo = new VariationalInference(initModel, prior);

        // log time used by algo.run()
        long startTime = System.currentTimeMillis();
        algo.run();
        long endTime = System.currentTimeMillis();
        long totalTimeMillis = endTime - startTime;

        VariationalModel posterior = algo.getVariationalPosterior();
        for (VariationalVariable var:posterior.getNodeHeightVariableList()) {
            System.out.println(var.getMean());
            System.out.println(var.getStandardDeviation());
        }
        for (VariationalVariable var:posterior.getPopSizeVariableList()) {
            System.out.println(var.getMean());
            System.out.println(var.getStandardDeviation());
        }

        System.out.println("");
        System.out.println("Total execution time: " + totalTimeMillis / 1000.0 + " s");
        System.out.println("Time used to build HMM: " + Utils.buildingTime / 1000.0 + " s");
        System.out.println("Time used in forward algorithm: " + Utils.likelihoodTime / 1000.0 + " s");
    }
}
