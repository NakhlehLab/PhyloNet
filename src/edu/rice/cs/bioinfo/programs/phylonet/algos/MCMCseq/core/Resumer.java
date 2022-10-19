package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core;
/*
 * @ClassName:   Resumer
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        11/22/21 10:55 AM
 */

import com.google.common.base.Strings;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Tools;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class Resumer {
    /* Constructor */
    public Resumer() {
    }

    public static void main(String[] args) {

        multiRuns(args[0], args[1],
                Integer.parseInt(args[2]),
                Integer.parseInt(args[3]), Integer.parseInt(args[4]),
                Integer.parseInt(args[5]), Integer.parseInt(args[6]),
                Integer.parseInt(args[7]), Integer.parseInt(args[8]), args[9],
                Boolean.parseBoolean(args[10]),
                Boolean.parseBoolean(args[11]),
                Boolean.parseBoolean(args[12]),
                Integer.parseInt(args[13]));

//        multiRuns("/Users/zhen/Desktop/Zhen/research/phylogenetics/butterfly/data/trinets_output/mcmc1/",
//                "/Users/zhen/Desktop/Zhen/research/phylogenetics/butterfly/data/trinets_output/mcmc1/mcmc_1.nex",
//                40000000, 10000000, 2, 5000, 0, 20, 0, null, true, false, true, 50);
    }

    /**
     * Divide one long run of MCMC into multiple subruns.
     * Subrun i+1 will resume from the last state of subrun i.
     * @param basePath
     * @param nexFilePath
     * @param chainLen
     * @param burnInLen
     * @param subrunIndex
     * @param sampleFreq
     * @param lociSetId
     * @param subrunNum
     * @param seedID
     * @param gtrRates
     */
    public static void multiRuns(String basePath, String nexFilePath, int chainLen, int burnInLen, int subrunIndex, int sampleFreq, int lociSetId, int subrunNum, int seedID, String gtrRates, boolean sample_murate, boolean const_popsize, boolean diameter_prior, int locisubset) {
        Utils._NUM_THREADS = Runtime.getRuntime().availableProcessors();
        Utils._CHAIN_LEN = chainLen;
        Utils._BURNIN_LEN = burnInLen;
        Utils._SubrunNum = subrunNum;
        Utils._SubrunIndex = subrunIndex;
        Utils._SAMPLE_FREQUENCY = sampleFreq;
        Utils._SAMPLE_NUM = chainLen / sampleFreq;
        Utils.SAMPLE_MUTATION_RATE = sample_murate;
        Utils._CONST_POP_SIZE = const_popsize;
        Utils._DIAMETER_PRIOR = diameter_prior;

        // set substitution rates
        if ((!Strings.isNullOrEmpty(gtrRates)) && (!"null".equals(gtrRates))) {
            String[] tenRates = gtrRates.split(",");
            if (tenRates.length != 10) {
                throw new RuntimeException();
            }
            double[] baseFreqs = new double[4];
            double[] transRates = new double[6];

            for (int i = 0; i < tenRates.length; i++) {
                double r = Double.parseDouble(tenRates[i].trim());
                if (r < 0) {
                    throw new RuntimeException();
                }
                if (i >= 4) {
                    transRates[i - 4] = r;
                } else {
                    baseFreqs[i] = r;
                }
            }
            Utils._SUBSTITUTION_MODEL = "GTR";
            Utils._BASE_FREQS = baseFreqs;
            Utils._TRANS_RATES = transRates;
        }

        Utils._TRUE_START = false;
        Utils._PreRun = true;
        Utils.resourceNum = 0;
//        Utils._POP_SIZE_MEAN = 0.4;
        Utils._IN_DIRECTORY = basePath;
        Utils._MC3_CHAINS = new ArrayList<>();
        Utils._OUT_DIRECTORY = basePath + "mcmcseq/" + String.valueOf(lociSetId) + "_" + String.valueOf(seedID) + "/";

        Utils.isExistsDir(Utils._OUT_DIRECTORY);
        System.out.println(Utils._OUT_DIRECTORY);
        Utils._CGT_DIRECTORY = Utils._OUT_DIRECTORY;
        Utils._OUT_DIRECTORY += String.valueOf(chainLen) + "_" + String.valueOf(subrunIndex) + "/";
        Utils.isExistsDir(Utils._OUT_DIRECTORY);
        Map<String, Map<String,String>> multiLociSeq = Tools.parseNexusFile(nexFilePath, locisubset);


        List<Alignment> alns = Tools.readSeq(multiLociSeq);
        Collections.sort(alns);

        // Init chain
        if (subrunIndex == 0) {
            Utils._SEED = Tools.getSeed(seedID);
            Utils._START_STATE = basePath + "output/" ;
            Utils._START_TIME = System.currentTimeMillis();
            MC3Core initRun = new MC3Core(alns);
            initRun.firstRun();
            System.out.println(String.format("Total elapsed time : %2.5f s\n",
                    (double) (System.currentTimeMillis() - Utils._START_TIME) / 1000.0));
            System.out.println(Utils._OUT_DIRECTORY);
        }
        // resume the chain from last subrun
        else if (subrunIndex < subrunNum) {
            Utils._START_STATE = Utils._CGT_DIRECTORY + String.valueOf(chainLen) + "_" ;
            Utils._START_TIME = System.currentTimeMillis();
            MC3Core multiRun = new MC3Core(alns,
                    Utils._START_STATE + String.valueOf(Utils._SubrunIndex - 1));
            multiRun.multipleRun();
            System.out.println(String.format("Total elapsed time : %2.5f s",
                    (double) (System.currentTimeMillis() - Utils._START_TIME) / 1000.0));
            System.out.println("Number of processors used: " + Utils._NUM_THREADS);
            System.out.println(Utils._OUT_DIRECTORY);
        }

    }
}
