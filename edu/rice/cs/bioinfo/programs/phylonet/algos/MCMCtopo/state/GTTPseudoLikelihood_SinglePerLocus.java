package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by dw20 on 5/11/17.
 */
public class GTTPseudoLikelihood_SinglePerLocus extends NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus
        implements GTTLikelihood{

    public double computeProbability(Network<Object> speciesNetwork, List allTriplets,
                                     Map<String,List<String>> species2alleles, List tripleFrequencies) {
        return super.computeProbability(speciesNetwork, allTriplets, tripleFrequencies, species2alleles);
    }

    public void summarizeData(List originalGTs, Map<String,String> allele2species,
                              List gtsForStartingNetwork, List allTriplets, List tripletFrequencies){
        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        for(Object list: originalGTs) {
            for (MutableTuple<Tree, Double> gtTuple : (List<MutableTuple<Tree, Double>>)list) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
                if (existingTuple == null) {
                    exp2tree.put(exp, new MutableTuple<>(gtTuple.Item1, gtTuple.Item2));

                } else {
                    existingTuple.Item2 += gtTuple.Item2;
                }

            }
        }
        gtsForStartingNetwork.addAll(exp2tree.values());
        NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus.computeTripleFrequenciesInGTs(
                gtsForStartingNetwork, allele2species, allTriplets, tripletFrequencies);
    }

}
