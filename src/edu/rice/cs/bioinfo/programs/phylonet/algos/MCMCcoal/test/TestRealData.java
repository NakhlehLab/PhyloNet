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
        butterfly();
//        /*
//         * load data
//         */
//        Map<String, String> nameTranslate = new HashMap<>();
//        nameTranslate.put("human", "H");
//        nameTranslate.put("chimp", "C");
//        nameTranslate.put("gorilla", "G");
//        // Read file
//        Map<String, String> omap = new HashMap<>();
//        File file = new File("/Users/xinhaoliu/Desktop/Research/Data/hobolth/output.4.fasta");
//        BufferedReader br = new BufferedReader(new FileReader(file));
//        String line;
//        String taxon = null;
//        while ((line = br.readLine()) != null) {
//            if (line.charAt(0) == '>') {
//                taxon = line.substring(1);
//            } else {
//                omap.put(nameTranslate.get(taxon), line);
//            }
//        }
//        Utils.DATA = new Alignment(omap);
//
//        /*
//         * initialize inference
//         */
//        ModelTree initModel = HCGModelBuilder.getHCGModelInit();
//        Prior prior = new Prior();
//        VariationalInference algo = new VariationalInference(initModel, prior);
//
//        algo.run();
//
//        VariationalModel posterior = algo.getVariationalPosterior();
//
//        for (VariationalVariable var:posterior.getNodeHeightVariableList()) {
//            System.out.println(var.getMean());
//            System.out.println(var.getStandardDeviation());
//        }
//        for (VariationalVariable var:posterior.getPopSizeVariableList()) {
//            System.out.println(var.getMean());
//            System.out.println(var.getStandardDeviation());
//        }
    }

    private static void butterfly() throws IOException {
        /*
         * load data
         */
        Map<String, String> nameTranslate = new HashMap<>();
        nameTranslate.put("Hcyd", "cydno");
        nameTranslate.put("Hnum", "numata");
        nameTranslate.put("Htim", "timareta");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/butterfly/threetaxa/Hmel213020_1_Hmel213020_singleCopy.fa");
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
         * Print hyperparameters for bookkeeping
         */
        System.out.println("213020_1");
        System.out.println("2e-9 mutation rate");
        System.out.println("N0 is " + Utils.N0);
        System.out.println(Utils.NUM_BIN + "-" + Utils.CROSS_OVER_RATE);
        System.out.println("NODE_HEIGHT_INIT_STDDEV = " + Utils.NODE_HEIGHT_INIT_STDDEV);
        System.out.println("POP_SIZE_INIT_STDDEV = " + Utils.POP_SIZE_INIT_STDDEV);
        System.out.println("NODE_HEIGHT_MIN_STDDEV = " + Utils.NODE_HEIGHT_MIN_STDDEV);
        System.out.println("POP_SIZE_MIN_STDDEV = " + Utils.POP_SIZE_MIN_STDDEV);
        System.out.println("NODE_HEIGHT_MEAN_LEARNING_RATE = " + Utils.NODE_HEIGHT_MEAN_LEARNING_RATE);
        System.out.println("NODE_HEIGHT_STDDEV_LEARNING_RATE = " + Utils.NODE_HEIGHT_STDDEV_LEARNING_RATE);
        System.out.println("POP_SIZE_MEAN_LEARNING_RATE = " + Utils.POP_SIZE_MEAN_LEARNING_RATE);
        System.out.println("POP_SIZE_STDDEV_LEARNING_RATE = " + Utils.POP_SIZE_STDDEV_LEARNING_RATE);

        /*
         * initialize inference
         */
        ModelTree initModel = NewButterflyModelBuilder.getButterflyModelRealData();
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
