package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmCore;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.Prior;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.PriorityQueue;

public class TestIllegalModel {
    public static void main(String[] args) throws IOException {
//        Map<String, String> msName2SpeciesName = new HashMap<>();
//        msName2SpeciesName.put("1", "H");
//        msName2SpeciesName.put("2", "C");
//        msName2SpeciesName.put("3", "G");
//        // Read file
//        Map<String, String> omap = new HashMap<>();
//        File file = new File("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/20recomb20mut/HCG_100000/aligned.fasta");
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
//        // now test negative node height and negative pop size
//        // now test negative recomb rate
//        ModelTree model = HCGModelBuilder.getHCGModelIllegal();
//        HmmBuilder builder = new HmmBuilder(model.getTree(), model.getRecombRate());
//        HmmCore hmm = builder.build();
//        double likelihood = hmm.logLikelihood();
//        System.out.println(likelihood);
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

        ModelTree model = ButterflyModelBuilder.getButterflyModelMyValue();
        HmmBuilder builder = new HmmBuilder(model.getTree(), model.getRecombRate());
        HmmCore hmm = builder.build();
        double likelihood = hmm.logLikelihood();
        System.out.println(likelihood);

        Prior prior = new Prior();
        System.out.println(prior.logPrior(model));
    }
}
