package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.distribution;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkPseudoLikelihoodFromGTT;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 2019-04-25
 * Time: 11:55
 * To change this template use File | Settings | File Templates.
 */
public class GeneTreeBrSpeciesNetPseudoLikelihood {
    NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus likelihood = null;
    Map<String,String> allele2species = null;
    Map<String, List<String>> species2alleles = null;

    public GeneTreeBrSpeciesNetPseudoLikelihood() {
        likelihood = new NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus();

        if(Utils._TAXON_MAP != null) {
            allele2species = new HashMap<>();
            for(String species : Utils._TAXON_MAP.keySet()) {
                for(String allele : Utils._TAXON_MAP.get(species)) {
                    allele2species.put(allele, species);
                }
            }

            species2alleles = Utils._TAXON_MAP;
        }

    }

    public double calculateGTDistribution(UltrametricNetwork network, List<UltrametricTree> geneTrees){
        List allTriplets = new ArrayList();
        List tripletFrequencies = new ArrayList();

        List<List<MutableTuple<Tree, Double>>> gts = new ArrayList<>();
        gts.add(new ArrayList<>());
        for(UltrametricTree ut : geneTrees) {
            gts.get(0).add(new MutableTuple<>(Trees.readTree(ut.getTree().toNewick()), 1.0) );
        }
        likelihood.summarizeData(gts, allele2species, allTriplets, tripletFrequencies);

        Network<Object> speciesNetwork = (Network<Object>) Networks.readNetworkWithRootPop(Networks.getFullString( network.getNetwork()));
        return likelihood.computeProbability(speciesNetwork, allTriplets, tripletFrequencies, species2alleles);
    }
}
