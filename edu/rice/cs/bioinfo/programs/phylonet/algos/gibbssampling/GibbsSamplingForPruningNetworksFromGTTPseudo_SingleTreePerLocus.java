package edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityPseudo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by yunyu on 10/1/15.
 */
public class GibbsSamplingForPruningNetworksFromGTTPseudo_SingleTreePerLocus extends GibbsSamplingForPruningNetworksFromGTTPseudo {

    public void summarizeGTs(List<List<MutableTuple<Tree,Double>>> originalGTs, Map<String,String> allele2species, List allTriplets, List tripletFrequencies){
        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        for(Object list: originalGTs) {
            for (MutableTuple<Tree, Double> gtTuple : (List<MutableTuple<Tree, Double>>)list) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
                if (existingTuple == null) {
                    existingTuple = gtTuple;
                    exp2tree.put(exp, existingTuple);

                } else {
                    existingTuple.Item2 += gtTuple.Item2;
                }

            }
        }
        List<MutableTuple<Tree,Double>> gts = new ArrayList<>();
        gts.addAll(exp2tree.values());
        NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus.computeTripleFrequenciesInGTs(gts, allele2species, allTriplets, tripletFrequencies);
    }



    protected double computeLikelihood(double[][] probs, List tripletFrequencies){
        int index = 0;
        double totalProb = 0;
        for(Object o: tripletFrequencies){
            double[] freqs = (double[])o;
            for(int i=0; i<3; i++){
                totalProb += freqs[i] * Math.log(probs[index][i]);
            }
            index++;
        }
        return totalProb;
    }



}
