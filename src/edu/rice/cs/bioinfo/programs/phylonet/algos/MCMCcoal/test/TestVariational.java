package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test;


import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmCore;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.VariationalInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.VariationalInferenceReparam;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.VariationalModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.VariationalModelReparam;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.Prior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.VariationalVariable;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Test the variational inference procedure.
 * Created by Xinhao Liu on 3/26/20.
 */
public class TestVariational {
    public static void main(String[] args) throws IOException {
        testHCG();
        //convergencePlot(3);
        //convergencePlot(4);
        //testHCGreparam();
    }

    private static void testDaC_HCO() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "H");
        msName2SpeciesName.put("2", "C");
        msName2SpeciesName.put("3", "G");
        msName2SpeciesName.put("4", "O");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/HCGO_200000_10recomb50mut/aligned.fasta");
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
        omap.remove("G");
        Utils.DATA = new Alignment(omap);
        /*
         * initialize inference
         */
        ModelTree initModel = DaCModelBuilder.getHCOModelInit();
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

    private static void testDaC_HGO() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "H");
        msName2SpeciesName.put("2", "C");
        msName2SpeciesName.put("3", "G");
        msName2SpeciesName.put("4", "O");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/HCGO_200000_10recomb50mut/aligned.fasta");
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
        omap.remove("C");
        Utils.DATA = new Alignment(omap);
        /*
         * initialize inference
         */
        ModelTree initModel = DaCModelBuilder.getHGOModel();
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

    private static void testDaC_HCG() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "H");
        msName2SpeciesName.put("2", "C");
        msName2SpeciesName.put("3", "G");
        msName2SpeciesName.put("4", "O");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/HCGO_200000_10recomb50mut/aligned.fasta");
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
        omap.remove("O");
        Utils.DATA = new Alignment(omap);
        /*
         * initialize inference
         */
        ModelTree initModel = DaCModelBuilder.getHCGModelInit();
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

    private static void testHCGO() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "H");
        msName2SpeciesName.put("2", "C");
        msName2SpeciesName.put("3", "G");
        msName2SpeciesName.put("4", "O");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/HCGO_200000_10recomb50mut/aligned.fasta");
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
        ModelTree initModel = HCGOModelBuilder.getHCGOModelInit();
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

    private static void testNewButterfly() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "silvaniform");
        msName2SpeciesName.put("2", "melpomene");
        msName2SpeciesName.put("3", "cyd.gal.");
        msName2SpeciesName.put("4", "pachinus");
        // Read file
        Map<String, String> omap = new HashMap<>();
//        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/butterfly_100000_new/aligned.fasta");
        File file = new File("/scratch/xl59/data/butterfly/new.fasta"); // on NOTS
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
        ModelTree initModel = NewButterflyModelBuilder.getButterflyModelInit();
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

    private static void testButterfly() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "silvaniform");
        msName2SpeciesName.put("2", "tim.ssp_nor.");
        msName2SpeciesName.put("3", "mel.ama.");
        msName2SpeciesName.put("4", "mel.agl.");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/butterfly_100000/aligned.fasta");
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
        ModelTree initModel = ButterflyModelBuilder.getButterflyModelMyValue();
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

    private static void testHCG() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "H");
        msName2SpeciesName.put("2", "C");
        msName2SpeciesName.put("3", "G");
        //msName2SpeciesName.put("4", "O");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/69/aligned.fasta");
        //File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/aligned.fasta");
        //File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/HCG10mut_500000/aligned.fasta");
        //File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/10recomb50mut/HCG_100000/aligned.fasta");
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
        //omap.remove("O");
        Utils.DATA = new Alignment(omap);

        System.out.println("HCG_HCG 69");
        /*
         * initialize inference
         */
        ModelTree initModel = HCGModelBuilder.getHCGModelInit();
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

    private static void testHCGreparam() throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "H");
        msName2SpeciesName.put("2", "C");
        msName2SpeciesName.put("3", "G");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/100/aligned.fasta");
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

        System.out.println("HCG_HCG 100");
        /*
         * initialize inference
         */
        ModelTree initModel = HCGModelBuilder.getHCGModelInit();
        Prior prior = new Prior();
        VariationalInferenceReparam algo = new VariationalInferenceReparam(initModel, prior);

        // log time used by algo.run()
        long startTime = System.currentTimeMillis();
        algo.run();
        long endTime = System.currentTimeMillis();
        long totalTimeMillis = endTime - startTime;

        VariationalModelReparam posterior = algo.getVariationalPosterior();
        for (VariationalVariable var:posterior.getBranchLengthVariableList()) {
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

    /**
     * Plot likelihood convergence from trace file on 100 simulated dataset.
     */
    private static void convergencePlot(int experimentIdx) throws IOException {
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "H");
        msName2SpeciesName.put("2", "C");
        msName2SpeciesName.put("3", "G");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/" + experimentIdx + "/aligned.fasta");
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
         * Read trace file
         */
        List<Double> T_HC_trace = new ArrayList<>();
        List<Double> T_HCG_trace = new ArrayList<>();
        List<Double> N_HC_trace = new ArrayList<>();
        List<Double> N_HCG_trace = new ArrayList<>();
        File outfile = new File("/Users/xinhaoliu/Desktop/Research/Experiments/coalHMM/HCG_sim/HCG20datasets/" + experimentIdx + "/1");
        BufferedReader outBr = new BufferedReader(new FileReader(outfile));
        String outLine;
        int index = -1;
        while ((outLine = outBr.readLine()) != null) {
            if (outLine.startsWith("=")) {
                index = 0;
            }
            if (index >= 0) {
                if (index == 1) {
                    T_HC_trace.add(Double.valueOf(outLine.substring(0, outLine.length() - 1)));
                } else if (index == 7) {
                    T_HCG_trace.add(Double.valueOf(outLine.substring(0, outLine.length() - 1)));
                } else if (index == 14) {
                    N_HC_trace.add(Double.valueOf(outLine.substring(0, outLine.length() - 1)));
                } else if (index == 20) {
                    N_HCG_trace.add(Double.valueOf(outLine.substring(0, outLine.length() - 1)));
                }
                index++;
            }
        }

        List<Double> likelihoodTrace = new ArrayList<>();
        for (int i = 0; i < T_HC_trace.size(); i++) {
            double T_HC = T_HC_trace.get(i);
            double T_HCG = T_HCG_trace.get(i);
            double N_HC = N_HC_trace.get(i);
            double N_HCG = N_HCG_trace.get(i);

            int numIter = 10;
            List<Double> likelihoods = new ArrayList<>();
            for (int k = 0; k < numIter; k++) {
                ModelTree model = HCGModelBuilder.getHCGModelSetValue(T_HC, T_HCG, N_HC, N_HCG);
                HmmBuilder builder = new HmmBuilder(model.getTree(), model.getRecombRate());
                HmmCore hmm = builder.build();
                double likelihood = hmm.logLikelihood();
                likelihoods.add(likelihood);
            }

            double avgLikelihood = likelihoods.stream().mapToDouble(x -> x).sum() / numIter;
            likelihoodTrace.add(avgLikelihood);
            //System.out.println("Done with iteration " + i);
        }
        System.out.println(experimentIdx);
        System.out.println(likelihoodTrace);
    }
}
