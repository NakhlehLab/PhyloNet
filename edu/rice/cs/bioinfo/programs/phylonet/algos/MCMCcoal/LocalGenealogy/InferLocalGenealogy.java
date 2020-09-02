package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.LocalGenealogy;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmCore;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test.HCGModelBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Posterior decoding of local genealogies.
 *
 * Created by Xinhao Liu on 7/8/20.
 */
public class InferLocalGenealogy {
    public static void main(String[] args) throws IOException, ParseException {
        runAllDataset();
    }

    public static void dataset2() throws IOException, ParseException {
        String[] stateNames = new String[] {"HC1", "HC2", "HG", "CG"};
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
            if (line.length() == 2) {
                taxon = line.substring(1);
            } else {
                omap.put(msName2SpeciesName.get(taxon), line);
            }
        }
        Utils.DATA = new Alignment(omap);

        // Build model and do posterior decoding
        ModelTree trueModel = HCGModelBuilder.getHCGModel();
        HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
        HmmCore hmm = builder.build();
        double[][] posteriorMatrix = hmm.posteriorDecoding();
        double[] HC1 = new double[Utils.DATA.getSiteCount()];
        double[] HC2 = new double[Utils.DATA.getSiteCount()];
        double[] HG = new double[Utils.DATA.getSiteCount()];
        double[] CG = new double[Utils.DATA.getSiteCount()];

        // Initialize confusion matrix as a map of map
        Map<String, Map<String, Integer>> confusionMatrix = new HashMap<>();
        for (String trueState:stateNames) {
            confusionMatrix.put(trueState, new HashMap<>());
            for (String inferredState:stateNames) {
                confusionMatrix.get(trueState).put(inferredState, 0);
            }
        }
        // Parse true states
        List<String> trueStateSequence = ControlFileParser.parse("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/2/control.txt");
        // Iterate over sites
        for (int siteIdx = 0; siteIdx < Utils.DATA.getSiteCount(); siteIdx++) {
            double maxProb = -1;
            String MAP = "ERROR";
            for (int stateIdx = 0; stateIdx < hmm.getNumStates(); stateIdx++) {
                switch (hmm.getState(stateIdx).getStateName()) {
                    case "HC1":
                        HC1[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        if (HC1[siteIdx] > maxProb) {
                            maxProb = HC1[siteIdx];
                            MAP = "HC1";
                        }
                        break;
                    case "HC2":
                        HC2[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        if (HC2[siteIdx] > maxProb) {
                            maxProb = HC2[siteIdx];
                            MAP = "HC2";
                        }
                        break;
                    case "HG":
                        HG[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        if (HG[siteIdx] > maxProb) {
                            maxProb = HG[siteIdx];
                            MAP = "HG";
                        }
                        break;
                    case "CG":
                        CG[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        if (CG[siteIdx] > maxProb) {
                            maxProb = CG[siteIdx];
                            MAP = "CG";
                        }
                        break;
                    default:
                        System.out.println("ERROR IN InferLocalGenealogy!!!!");
                        break;
                }
            }
            String trueState = trueStateSequence.get(siteIdx);
            confusionMatrix.get(trueState).put(MAP, confusionMatrix.get(trueState).get(MAP) + 1);
        }
        // Print results
        for (String trueState:stateNames) {
            for (String inferredState:stateNames) {
                System.out.println("True: " + trueState + ", Inferred: " + inferredState + ", count = " + confusionMatrix.get(trueState).get(inferredState));
            }
        }
//        BufferedWriter HC1writer = new BufferedWriter(new FileWriter("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/2/HC1"));
//        for (double prob:HC1) {
//            HC1writer.write(String.valueOf(prob));
//            HC1writer.newLine();
//        }
//        HC1writer.close();
//
//        BufferedWriter HC2writer = new BufferedWriter(new FileWriter("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/2/HC2"));
//        for (double prob:HC2) {
//            HC2writer.write(String.valueOf(prob));
//            HC2writer.newLine();
//        }
//        HC2writer.close();
//
//        BufferedWriter HGwriter = new BufferedWriter(new FileWriter("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/2/HG"));
//        for (double prob:HG) {
//            HGwriter.write(String.valueOf(prob));
//            HGwriter.newLine();
//        }
//        HGwriter.close();
//
//        BufferedWriter CGwriter = new BufferedWriter(new FileWriter("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/2/CG"));
//        for (double prob:CG) {
//            CGwriter.write(String.valueOf(prob));
//            CGwriter.newLine();
//        }
//        CGwriter.close();
    }

    public static void runAllDataset() throws IOException, ParseException {
        String[] stateNames = new String[] {"HC1", "HC2", "HG", "CG"};
        // Initialize confusion matrix as a map of map
        Map<String, Map<String, Integer>> confusionMatrix = new HashMap<>();
        for (String trueState:stateNames) {
            confusionMatrix.put(trueState, new HashMap<>());
            for (String inferredState:stateNames) {
                confusionMatrix.get(trueState).put(inferredState, 0);
            }
        }

        for (int datasetIdx = 1; datasetIdx <= 100; datasetIdx++) {
            /*
             * load data
             */
            Map<String, String> msName2SpeciesName = new HashMap<>();
            msName2SpeciesName.put("1", "H");
            msName2SpeciesName.put("2", "C");
            msName2SpeciesName.put("3", "G");
            // Read file
            Map<String, String> omap = new HashMap<>();
            File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/" + datasetIdx + "/aligned.fasta");
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

            // Build model and do posterior decoding
            ModelTree trueModel = HCGModelBuilder.getHCGModel();
            HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
            HmmCore hmm = builder.build();
            double[][] posteriorMatrix = hmm.posteriorDecoding();
            double[] HC1 = new double[Utils.DATA.getSiteCount()];
            double[] HC2 = new double[Utils.DATA.getSiteCount()];
            double[] HG = new double[Utils.DATA.getSiteCount()];
            double[] CG = new double[Utils.DATA.getSiteCount()];

            // Parse true states
            List<String> trueStateSequence = ControlFileParser.parse("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/" + datasetIdx + "/control.txt");
            // Iterate over sites
            for (int siteIdx = 0; siteIdx < Utils.DATA.getSiteCount(); siteIdx++) {
                double maxProb = -1;
                String MAP = "ERROR";
                for (int stateIdx = 0; stateIdx < hmm.getNumStates(); stateIdx++) {
                    switch (hmm.getState(stateIdx).getStateName()) {
                        case "HC1":
                            HC1[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                            if (HC1[siteIdx] > maxProb) {
                                maxProb = HC1[siteIdx];
                                MAP = "HC1";
                            }
                            break;
                        case "HC2":
                            HC2[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                            if (HC2[siteIdx] > maxProb) {
                                maxProb = HC2[siteIdx];
                                MAP = "HC2";
                            }
                            break;
                        case "HG":
                            HG[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                            if (HG[siteIdx] > maxProb) {
                                maxProb = HG[siteIdx];
                                MAP = "HG";
                            }
                            break;
                        case "CG":
                            CG[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                            if (CG[siteIdx] > maxProb) {
                                maxProb = CG[siteIdx];
                                MAP = "CG";
                            }
                            break;
                        default:
                            System.out.println("ERROR IN InferLocalGenealogy!!!!");
                            break;
                    }
                }
                String trueState = trueStateSequence.get(siteIdx);
                confusionMatrix.get(trueState).put(MAP, confusionMatrix.get(trueState).get(MAP) + 1);
            }
            System.out.println("Finished dataset " + datasetIdx);
//            // Print results
//            System.out.println("==========================" + "DATASET " + datasetIdx + "==========================");
//            for (String trueState:stateNames) {
//                for (String inferredState:stateNames) {
//                    System.out.println("True: " + trueState + ", Inferred: " + inferredState + ", count = " + confusionMatrix.get(trueState).get(inferredState));
//                }
//            }
        }
        // Print results
        System.out.println("==========================FINAL==========================");
        for (String trueState:stateNames) {
            for (String inferredState:stateNames) {
                System.out.println("True: " + trueState + ", Inferred: " + inferredState + ", count = " + confusionMatrix.get(trueState).get(inferredState));
            }
        }
    }

    public static void writeOutTrueGenealogy(String path) throws IOException, ParseException {
        // Parse true states
        List<String> trueStateSequence = ControlFileParser.parse(path + "/control.txt");
        // Write
        BufferedWriter writer = new BufferedWriter(new FileWriter(path + "/TrueSequence"));
        for (String trueGenealogy:trueStateSequence) {
            writer.write(trueGenealogy);
            writer.newLine();
        }
        writer.close();
    }

    public static void writeOutPosteriorDecoding(String path) throws IOException, ParseException {
        String[] stateNames = new String[] {"HC1", "HC2", "HG", "CG"};
        /*
         * load data
         */
        Map<String, String> msName2SpeciesName = new HashMap<>();
        msName2SpeciesName.put("1", "H");
        msName2SpeciesName.put("2", "C");
        msName2SpeciesName.put("3", "G");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File(path + "/aligned.fasta");
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

        // Build model and do posterior decoding
        ModelTree trueModel = HCGModelBuilder.getHCGModel();
        HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
        HmmCore hmm = builder.build();
        double[][] posteriorMatrix = hmm.posteriorDecoding();
        double[] HC1 = new double[Utils.DATA.getSiteCount()];
        double[] HC2 = new double[Utils.DATA.getSiteCount()];
        double[] HG = new double[Utils.DATA.getSiteCount()];
        double[] CG = new double[Utils.DATA.getSiteCount()];

        // Initialize confusion matrix as a map of map
        Map<String, Map<String, Integer>> confusionMatrix = new HashMap<>();
        for (String trueState:stateNames) {
            confusionMatrix.put(trueState, new HashMap<>());
            for (String inferredState:stateNames) {
                confusionMatrix.get(trueState).put(inferredState, 0);
            }
        }
        // Parse true states
        List<String> trueStateSequence = ControlFileParser.parse(path + "/control.txt");
        // Iterate over sites
        for (int siteIdx = 0; siteIdx < Utils.DATA.getSiteCount(); siteIdx++) {
            double maxProb = -1;
            String MAP = "ERROR";
            for (int stateIdx = 0; stateIdx < hmm.getNumStates(); stateIdx++) {
                switch (hmm.getState(stateIdx).getStateName()) {
                    case "HC1":
                        HC1[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        if (HC1[siteIdx] > maxProb) {
                            maxProb = HC1[siteIdx];
                            MAP = "HC1";
                        }
                        break;
                    case "HC2":
                        HC2[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        if (HC2[siteIdx] > maxProb) {
                            maxProb = HC2[siteIdx];
                            MAP = "HC2";
                        }
                        break;
                    case "HG":
                        HG[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        if (HG[siteIdx] > maxProb) {
                            maxProb = HG[siteIdx];
                            MAP = "HG";
                        }
                        break;
                    case "CG":
                        CG[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        if (CG[siteIdx] > maxProb) {
                            maxProb = CG[siteIdx];
                            MAP = "CG";
                        }
                        break;
                    default:
                        System.out.println("ERROR IN InferLocalGenealogy!!!!");
                        break;
                }
            }
            String trueState = trueStateSequence.get(siteIdx);
            confusionMatrix.get(trueState).put(MAP, confusionMatrix.get(trueState).get(MAP) + 1);
        }
        BufferedWriter HC1writer = new BufferedWriter(new FileWriter(path + "/HC1"));
        for (double prob:HC1) {
            HC1writer.write(String.valueOf(prob));
            HC1writer.newLine();
        }
        HC1writer.close();

        BufferedWriter HC2writer = new BufferedWriter(new FileWriter(path + "/HC2"));
        for (double prob:HC2) {
            HC2writer.write(String.valueOf(prob));
            HC2writer.newLine();
        }
        HC2writer.close();

        BufferedWriter HGwriter = new BufferedWriter(new FileWriter(path + "/HG"));
        for (double prob:HG) {
            HGwriter.write(String.valueOf(prob));
            HGwriter.newLine();
        }
        HGwriter.close();

        BufferedWriter CGwriter = new BufferedWriter(new FileWriter(path + "/CG"));
        for (double prob:CG) {
            CGwriter.write(String.valueOf(prob));
            CGwriter.newLine();
        }
        CGwriter.close();
    }

    public static void writeOutPDandTGForAllDataset() throws IOException, ParseException {
        for (int i = 1; i <= 75; i++) {
            if (i == 2 || i == 6 || i == 25) {
                continue;
            }
            String path = "/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/" + i;
            writeOutPosteriorDecoding(path);
            writeOutTrueGenealogy(path);
            System.out.println("Done with " + i);
        }
    }
}
