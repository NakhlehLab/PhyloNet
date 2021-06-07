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
 * Test bin number and simulation length using ABC models.
 * Created by Xinhao Liu on 6/16/20.
 */
public class TestABC {
    public static void main(String[] args) throws IOException {
        testABC4();
    }

    private static void testABC1() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "A");
        msName2SpeciesName.put("2", "B");
        msName2SpeciesName.put("3", "C");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/ABC/ABC1_100000/aligned.fasta");
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        String taxon = null;
        while ((line = br.readLine()) != null) {
            if (line.length() == 2) {
                taxon = line.substring(1);
            } else {
                omap.put(msName2SpeciesName.get(taxon), line);
            }
        }
        Utils.DATA = new Alignment(omap);

        System.out.println("ABC1 5000 5");
        /*
         * initialize inference
         */
        ModelTree initModel = ABCModelBuilder.getABC1Model();
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

    private static void testABC2() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "A");
        msName2SpeciesName.put("2", "B");
        msName2SpeciesName.put("3", "C");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/ABC/ABC2_100000/aligned.fasta");
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        String taxon = null;
        while ((line = br.readLine()) != null) {
            if (line.length() == 2) {
                taxon = line.substring(1);
            } else {
                omap.put(msName2SpeciesName.get(taxon), line);
            }
        }
        Utils.DATA = new Alignment(omap);

        System.out.println("ABC2 3000 4");
        /*
         * initialize inference
         */
        ModelTree initModel = ABCModelBuilder.getABC2Model();
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

    private static void testABC4() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "A");
        msName2SpeciesName.put("2", "B");
        msName2SpeciesName.put("3", "C");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/ABC/ABC4_100000/aligned.fasta");
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        String taxon = null;
        while ((line = br.readLine()) != null) {
            if (line.length() == 2) {
                taxon = line.substring(1);
            } else {
                omap.put(msName2SpeciesName.get(taxon), line);
            }
        }
        Utils.DATA = new Alignment(omap);

        /*
         * Print hyperparameters for bookkeeping
         */
        System.out.println("long region simulation");
        System.out.println("ABC4 5000 5");
        System.out.println("1.25e-6 mutation rate");
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
        ModelTree initModel = ABCModelBuilder.getABC4Model();
        Prior prior = new Prior();
        VariationalInference algo = new VariationalInference(initModel, prior);

//        algo.run();
//
//        VariationalModel posterior = algo.getVariationalPosterior();
//        for (VariationalVariable var:posterior.getNodeHeightVariableList()) {
//            System.out.println(var.getMean());
//            System.out.println(var.getStandardDeviation());
//        }
//        for (VariationalVariable var:posterior.getPopSizeVariableList()) {
//            System.out.println(var.getMean());
//            System.out.println(var.getStandardDeviation());
//        }


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
