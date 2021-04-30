package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;

import java.io.*;
import java.util.*;

public class TestModelParam {
    public static void main(String[] args) throws IOException {
        testABC();
    }

    static void testHGO() throws IOException {
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

        int numIter = 5;

        List<Double> trueLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree trueModel = DaCModelBuilder.getHGOModel();
            HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
            HmmCore hmm = builder.build();
            System.out.println(hmm.getStates().size());
            double likelihood = hmm.logLikelihood();
            System.out.println(likelihood);
            trueLikelihoods.add(likelihood);
        }

        double avgTrueLikelihood = trueLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgTrueLikelihood);
        System.out.println("=================");

        // this model is wrong
        List<Double> wrongLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree slightlyWrongModel = DaCModelBuilder.getHGOModelMyValue();
            HmmBuilder slightlyWrongBuilder = new HmmBuilder(slightlyWrongModel.getTree(), slightlyWrongModel.getRecombRate());
            HmmCore slightlyWrongHmm = slightlyWrongBuilder.build();
            System.out.println(slightlyWrongHmm.getStates().size());
            double slightlyWrongLikelihood = slightlyWrongHmm.logLikelihood();
            System.out.println(slightlyWrongLikelihood);
            wrongLikelihoods.add(slightlyWrongLikelihood);
        }

        double avgWrongLikelihood = wrongLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgWrongLikelihood);
        System.out.println("=================");
    }

    static void testHCG() throws IOException {
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
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/3/aligned.fasta");
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
        //omap.remove("O");
        Utils.DATA = new Alignment(omap);

        int numIter = 5;

        List<Double> trueLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree trueModel = HCGModelBuilder.getHCGModel();
            HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
            HmmCore hmm = builder.build();
            System.out.println(hmm.getStates().size());
            double likelihood = hmm.logLikelihood();
            System.out.println(likelihood);
            trueLikelihoods.add(likelihood);
        }

        double avgTrueLikelihood = trueLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgTrueLikelihood);
        System.out.println("=================");

        // this model is wrong
        List<Double> wrongLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree slightlyWrongModel = HCGModelBuilder.getHCGModelMyValue();
            HmmBuilder slightlyWrongBuilder = new HmmBuilder(slightlyWrongModel.getTree(), slightlyWrongModel.getRecombRate());
            HmmCore slightlyWrongHmm = slightlyWrongBuilder.build();
            System.out.println(slightlyWrongHmm.getStates().size());
            double slightlyWrongLikelihood = slightlyWrongHmm.logLikelihood();
            System.out.println(slightlyWrongLikelihood);
            wrongLikelihoods.add(slightlyWrongLikelihood);
        }

        double avgWrongLikelihood = wrongLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgWrongLikelihood);
        System.out.println("=================");
    }

    static void testHCO() throws IOException {
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

        int numIter = 2;

        List<Double> trueLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree trueModel = DaCModelBuilder.getHCOModel();
            HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
            HmmCore hmm = builder.build();
            System.out.println(hmm.getStates().size());
            double likelihood = hmm.logLikelihood();
            System.out.println(likelihood);
            trueLikelihoods.add(likelihood);
        }

        double avgTrueLikelihood = trueLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgTrueLikelihood);
        System.out.println("=================");

        // this model is wrong
        List<Double> wrongLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree slightlyWrongModel = DaCModelBuilder.getHCOModelMyValue();
            HmmBuilder slightlyWrongBuilder = new HmmBuilder(slightlyWrongModel.getTree(), slightlyWrongModel.getRecombRate());
            HmmCore slightlyWrongHmm = slightlyWrongBuilder.build();
            System.out.println(slightlyWrongHmm.getStates().size());
            double slightlyWrongLikelihood = slightlyWrongHmm.logLikelihood();
            System.out.println(slightlyWrongLikelihood);
            wrongLikelihoods.add(slightlyWrongLikelihood);
        }

        double avgWrongLikelihood = wrongLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgWrongLikelihood);
        System.out.println("=================");
    }

    static void testHCGO() throws IOException {
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

        int numIter = 5;

        List<Double> trueLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree trueModel = HCGOModelBuilder.getHCGOModel();
            HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
            HmmCore hmm = builder.build();
            System.out.println(hmm.getStates().size());
            double likelihood = hmm.logLikelihood();
            System.out.println(likelihood);
            trueLikelihoods.add(likelihood);
        }

        double avgTrueLikelihood = trueLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgTrueLikelihood);
        System.out.println("=================");

        // this model is wrong
        List<Double> wrongLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree slightlyWrongModel = HCGOModelBuilder.getHCGOModelMyValue();
            HmmBuilder slightlyWrongBuilder = new HmmBuilder(slightlyWrongModel.getTree(), slightlyWrongModel.getRecombRate());
            HmmCore slightlyWrongHmm = slightlyWrongBuilder.build();
            System.out.println(slightlyWrongHmm.getStates().size());
            double slightlyWrongLikelihood = slightlyWrongHmm.logLikelihood();
            System.out.println(slightlyWrongLikelihood);
            wrongLikelihoods.add(slightlyWrongLikelihood);
        }

        double avgWrongLikelihood = wrongLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgWrongLikelihood);
        System.out.println("=================");
    }

    static void testNewButterfly() throws IOException {
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
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/butterfly_100000_new_upmut/aligned.fasta");
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

        int numIter = 2;

        List<Double> trueLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree trueModel = NewButterflyModelBuilder.getButterflyModel();
            HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
            HmmCore hmm = builder.build();
            System.out.println(hmm.getStates().size());
            double likelihood = hmm.logLikelihood();
//            System.out.println("calculated one likelihood");
            System.out.println(likelihood);
            trueLikelihoods.add(likelihood);
        }

        double avgTrueLikelihood = trueLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgTrueLikelihood);
        System.out.println("=================");

        // this model is wrong
        List<Double> wrongLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree slightlyWrongModel = NewButterflyModelBuilder.getButterflyModelMyValue2();
            HmmBuilder slightlyWrongBuilder = new HmmBuilder(slightlyWrongModel.getTree(), slightlyWrongModel.getRecombRate());
            HmmCore slightlyWrongHmm = slightlyWrongBuilder.build();
            System.out.println(slightlyWrongHmm.getStates().size());
            double slightlyWrongLikelihood = slightlyWrongHmm.logLikelihood();
            System.out.println(slightlyWrongLikelihood);
            wrongLikelihoods.add(slightlyWrongLikelihood);
        }

        double avgWrongLikelihood = wrongLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgWrongLikelihood);
        System.out.println("=================");
    }

    static void testButterfly() throws IOException {
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

        int numIter = 100;

//        List<Double> trueLikelihoods = new ArrayList<>();
//        for (int i = 0; i < numIter; i++) {
//            ModelTree trueModel = ButterflyModelBuilder.getButterflyModel();
//            HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
//            HmmCore hmm = builder.build();
//            double likelihood = hmm.logLikelihood();
//            System.out.println("calculated one likelihood");
////            System.out.println(likelihood);
//            trueLikelihoods.add(likelihood);
//        }
//
//        double avgTrueLikelihood = trueLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
//        System.out.println("=================");
//        System.out.println(avgTrueLikelihood);
//        System.out.println("=================");

        // this model is wrong
        List<Double> wrongLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree slightlyWrongModel = ButterflyModelBuilder.getButterflyModelMyValue();
            HmmBuilder slightlyWrongBuilder = new HmmBuilder(slightlyWrongModel.getTree(), slightlyWrongModel.getRecombRate());
            HmmCore slightlyWrongHmm = slightlyWrongBuilder.build();
            double slightlyWrongLikelihood = slightlyWrongHmm.logLikelihood();
            System.out.println(slightlyWrongLikelihood);
            wrongLikelihoods.add(slightlyWrongLikelihood);
        }

//        double avgWrongLikelihood = wrongLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
//        System.out.println("=================");
//        System.out.println(avgWrongLikelihood);
//        System.out.println("=================");
    }

    static void testABC() throws IOException {
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

        int numIter = 5;

        List<Long> simulationTimes = new ArrayList<>();
        List<Long> forwardAlgoTimes = new ArrayList<>();
        List<Double> numHiddenStates = new ArrayList<>();
        List<Long> computeCountTimes = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree trueModel = ABCModelBuilder.getABC4Model();
//            HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
            HmmBuilderShortLoci builder = new HmmBuilderShortLoci(trueModel.getTree(), trueModel.getRecombRate());
            long buildingStartTime = System.currentTimeMillis();
            HmmCore hmm = builder.build();
            simulationTimes.add(System.currentTimeMillis() - buildingStartTime);
            numHiddenStates.add((double) hmm.getStates().size());
            long likelihoodStartTime = System.currentTimeMillis();
            double likelihood = hmm.logLikelihood();
            forwardAlgoTimes.add(System.currentTimeMillis() - likelihoodStartTime);

            computeCountTimes.add(builder.computeCountTime);
        }

        long avgSimulationTime = simulationTimes.stream().mapToLong(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println("Time used to build HMM: " + avgSimulationTime / 1000.0 + " s");
        System.out.println("=================");
        long avgForwardAlgoTime = forwardAlgoTimes.stream().mapToLong(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println("Time used in forward algorithm: " + avgForwardAlgoTime / 1000.0 + " s");
        System.out.println("=================");
        double avgNumHiddenStates = numHiddenStates.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println("Number of hidden states: " + avgNumHiddenStates);
        System.out.println("=================");
        long avgComputeCountTime = computeCountTimes.stream().mapToLong(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println("Time used in computeCounts(): " + avgComputeCountTime / 1000.0 + " s");
        System.out.println("=================");
    }

//    static void testMSPTime() throws IOException {
//        /*
//         * load data
//         */
//        Map<String, String> msName2SpeciesName = new HashMap<>();
//        msName2SpeciesName.put("1", "A");
//        msName2SpeciesName.put("2", "B");
//        msName2SpeciesName.put("3", "C");
//        // Read file
//        Map<String, String> omap = new HashMap<>();
//        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/ABC/ABC4_100000/aligned.fasta");
//        BufferedReader br = new BufferedReader(new FileReader(file));
//        String line;
//        String taxon = null;
//        while ((line = br.readLine()) != null) {
//            if (line.length() == 2) {
//                taxon = line.substring(1);
//            } else {
//                omap.put(msName2SpeciesName.get(taxon), line);
//            }
//        }
//        Utils.DATA = new Alignment(omap);
//
//        int numIter = 5;
//
//        List<Long> msprimeTimes = new ArrayList<>();
//        for (int i = 0; i < numIter; i++) {
//            ModelTree trueModel = ABCModelBuilder.getABC4Model();
//            HmmBuilderShortLociTestMS builder = new HmmBuilderShortLociTestMS(trueModel.getTree(), trueModel.getRecombRate());
//            HmmCore hmm = builder.build();
//
//            msprimeTimes.add(builder.msprimeTime);
//        }
//
//        long avgmsprimeTime = msprimeTimes.stream().mapToLong(x -> x).sum() / numIter;
//        System.out.println("=================");
//        System.out.println("Time used by msprime: " + avgmsprimeTime / 1000.0 + " s");
//        System.out.println("=================");
//    }

    static void testNewButterflyRealData() throws IOException {
        /*
         * load data
         */
        Map<String, String> nameTranslate = new HashMap<>();
        nameTranslate.put("Hcyd", "cydno");
        nameTranslate.put("Hnum", "numata");
        nameTranslate.put("Htim", "timareta");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/butterfly/threetaxa/Hmel201007_1_Hmel201007_singleCopy.fa");
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

        int numIter = 5;

        List<Double> trueLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree trueModel = NewButterflyModelBuilder.getButterflyModelRealData();
            HmmBuilderShortLoci builder = new HmmBuilderShortLoci(trueModel.getTree(), trueModel.getRecombRate());
            HmmCore hmm = builder.build();
            System.out.println(hmm.getStates().size());
//            for (HiddenState state:hmm.getStates()) {
//                System.out.println(state.getTree().toNewick());
//            }
            double likelihood = hmm.logLikelihood();
            System.out.println(likelihood);
            trueLikelihoods.add(likelihood);
            System.exit(0);
        }

        double avgTrueLikelihood = trueLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgTrueLikelihood);
        System.out.println("=================");

        // this model is wrong
        List<Double> wrongLikelihoods = new ArrayList<>();
        for (int i = 0; i < numIter; i++) {
            ModelTree slightlyWrongModel = NewButterflyModelBuilder.getButterflyModelRealDataMyValue();
            HmmBuilderShortLoci slightlyWrongBuilder = new HmmBuilderShortLoci(slightlyWrongModel.getTree(), slightlyWrongModel.getRecombRate());
            HmmCore slightlyWrongHmm = slightlyWrongBuilder.build();
            System.out.println(slightlyWrongHmm.getStates().size());
//            for (HiddenState state:slightlyWrongHmm.getStates()) {
//                System.out.println(state.getTree().toNewick());
//            }
            double slightlyWrongLikelihood = slightlyWrongHmm.logLikelihood();
            System.out.println(slightlyWrongLikelihood);
            wrongLikelihoods.add(slightlyWrongLikelihood);
        }

        double avgWrongLikelihood = wrongLikelihoods.stream().mapToDouble(x -> x).sum() / numIter;
        System.out.println("=================");
        System.out.println(avgWrongLikelihood);
        System.out.println("=================");
    }
}
