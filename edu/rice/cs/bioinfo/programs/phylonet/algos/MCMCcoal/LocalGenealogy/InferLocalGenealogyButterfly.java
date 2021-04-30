package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.LocalGenealogy;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HiddenState;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmCore;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test.HCGModelBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test.NewButterflyModelBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Posterior decoding of local genealogies on real butterfly dataset.
 * Used for analyzing butterfly dataset only for the VICAR paper.
 *
 * Created by Xinhao Liu on 3/10/21.
 */
public class InferLocalGenealogyButterfly {
    public static void main(String[] args) throws IOException, ParseException {
//        writeOutPosteriorDecoding("/Users/xinhaoliu/Desktop/Research/Scripts/butterfly/vcf/210004_1/");
        writeOutPosteriorDecodingWithOutgroup("/Users/xinhaoliu/Desktop/Research/Scripts/butterfly/fourtaxa_vicar_ml/201007_1");
    }

    public static void writeOutPosteriorDecoding(String path) throws IOException, ParseException {
        String[] stateNames = new String[] {"CT1", "CT2", "CN", "TN"};
        /*
         * load data
         */
        Map<String, String> nameTranslate = new HashMap<>();
        nameTranslate.put("Hcyd", "cydno");
        nameTranslate.put("Hnum", "numata");
        nameTranslate.put("Htim", "timareta");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File(path + "/aligned.fasta");
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

        // Build model and do posterior decoding
        ModelTree trueModel = NewButterflyModelBuilder.getButterflyModelRealDataMyValue();
        HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
        HmmCore hmm = builder.build();
        double[][] posteriorMatrix = hmm.posteriorDecoding();
        double[] CT1 = new double[Utils.DATA.getSiteCount()];
        double[] CT2 = new double[Utils.DATA.getSiteCount()];
        double[] CN = new double[Utils.DATA.getSiteCount()];
        double[] TN = new double[Utils.DATA.getSiteCount()];

        // Iterate over sites
        for (int siteIdx = 0; siteIdx < Utils.DATA.getSiteCount(); siteIdx++) {
            for (int stateIdx = 0; stateIdx < hmm.getNumStates(); stateIdx++) {
                switch (hmm.getState(stateIdx).getStateNameButterfly(1.33690707698)) {
                    case "CT1":
                        CT1[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        break;
                    case "CT2":
                        CT2[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        break;
                    case "CN":
                        CN[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        break;
                    case "TN":
                        TN[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        break;
                    default:
                        System.out.println("ERROR IN InferLocalGenealogy!!!!");
                        System.exit(1);
                        break;
                }
            }
        }
        BufferedWriter CT1writer = new BufferedWriter(new FileWriter(path + "/" + Utils.NUM_BIN + "-" + Utils.CROSS_OVER_RATE + "/CT1"));
        for (double prob:CT1) {
            CT1writer.write(String.valueOf(prob));
            CT1writer.newLine();
        }
        CT1writer.close();

        BufferedWriter CT2writer = new BufferedWriter(new FileWriter(path + "/" + Utils.NUM_BIN + "-" + Utils.CROSS_OVER_RATE + "/CT2"));
        for (double prob:CT2) {
            CT2writer.write(String.valueOf(prob));
            CT2writer.newLine();
        }
        CT2writer.close();

        BufferedWriter CNwriter = new BufferedWriter(new FileWriter(path + "/" + Utils.NUM_BIN + "-" + Utils.CROSS_OVER_RATE + "/CN"));
        for (double prob:CN) {
            CNwriter.write(String.valueOf(prob));
            CNwriter.newLine();
        }
        CNwriter.close();

        BufferedWriter TNwriter = new BufferedWriter(new FileWriter(path + "/" + Utils.NUM_BIN + "-" + Utils.CROSS_OVER_RATE + "/TN"));
        for (double prob:TN) {
            TNwriter.write(String.valueOf(prob));
            TNwriter.newLine();
        }
        TNwriter.close();
    }

    public static void writeOutPosteriorDecodingWithOutgroup(String path) throws IOException, ParseException {
        String[] stateNames = new String[] {"CT1", "CT2", "CN", "TN"};
        /*
         * load data
         */
        Map<String, String> nameTranslate = new HashMap<>();
        nameTranslate.put("Hcyd", "cydno");
        nameTranslate.put("Hnum", "numata");
        nameTranslate.put("Htim", "timareta");
        nameTranslate.put("HeraRef", "erato");
        // Read file
        Map<String, String> omap = new HashMap<>();
        File file = new File(path + "/aligned.fasta");
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
        omap.remove("erato");
        Utils.DATA = new Alignment(omap);

        // Build model and do posterior decoding
        ModelTree trueModel = NewButterflyModelBuilder.getButterflyModelRealDataMyValue();
        HmmBuilder builder = new HmmBuilder(trueModel.getTree(), trueModel.getRecombRate());
        HmmCore hmm = builder.build();
        double[][] posteriorMatrix = hmm.posteriorDecoding();
        double[] CT1 = new double[Utils.DATA.getSiteCount()];
        double[] CT2 = new double[Utils.DATA.getSiteCount()];
        double[] CN = new double[Utils.DATA.getSiteCount()];
        double[] TN = new double[Utils.DATA.getSiteCount()];

        // Iterate over sites
        for (int siteIdx = 0; siteIdx < Utils.DATA.getSiteCount(); siteIdx++) {
            for (int stateIdx = 0; stateIdx < hmm.getNumStates(); stateIdx++) {
                switch (hmm.getState(stateIdx).getStateNameButterfly(1.26614376094)) {
                    case "CT1":
                        CT1[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        break;
                    case "CT2":
                        CT2[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        break;
                    case "CN":
                        CN[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        break;
                    case "TN":
                        TN[siteIdx] += posteriorMatrix[stateIdx][siteIdx];
                        break;
                    default:
                        System.out.println("ERROR IN InferLocalGenealogy!!!!");
                        System.exit(1);
                        break;
                }
            }
        }
        BufferedWriter CT1writer = new BufferedWriter(new FileWriter(path + "/" + Utils.NUM_BIN + "-" + Utils.CROSS_OVER_RATE + "/CT1"));
        for (double prob:CT1) {
            CT1writer.write(String.valueOf(prob));
            CT1writer.newLine();
        }
        CT1writer.close();

        BufferedWriter CT2writer = new BufferedWriter(new FileWriter(path + "/" + Utils.NUM_BIN + "-" + Utils.CROSS_OVER_RATE + "/CT2"));
        for (double prob:CT2) {
            CT2writer.write(String.valueOf(prob));
            CT2writer.newLine();
        }
        CT2writer.close();

        BufferedWriter CNwriter = new BufferedWriter(new FileWriter(path + "/" + Utils.NUM_BIN + "-" + Utils.CROSS_OVER_RATE + "/CN"));
        for (double prob:CN) {
            CNwriter.write(String.valueOf(prob));
            CNwriter.newLine();
        }
        CNwriter.close();

        BufferedWriter TNwriter = new BufferedWriter(new FileWriter(path + "/" + Utils.NUM_BIN + "-" + Utils.CROSS_OVER_RATE + "/TN"));
        for (double prob:TN) {
            TNwriter.write(String.valueOf(prob));
            TNwriter.newLine();
        }
        TNwriter.close();
    }
}
