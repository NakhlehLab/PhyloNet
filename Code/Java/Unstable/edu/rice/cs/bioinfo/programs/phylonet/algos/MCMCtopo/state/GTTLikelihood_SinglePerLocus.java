package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_Rooted;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkLikelihoodFromGTT_SingleTreePerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by wendingqiao on 6/6/15.
 */
public class GTTLikelihood_SinglePerLocus extends NetworkLikelihoodFromGTT_SingleTreePerLocus implements GTTLikelihood {

    public double computeProbability(Network<Object> speciesNetwork, List distinctTrees,
                                     Map<String,List<String>> species2alleles, List gtCorrespondences) {
        return super.computeProbability(speciesNetwork, distinctTrees, gtCorrespondences, species2alleles);
    }


    public void summarizeData(List originalGTs, Map<String,String> allele2species,
                                 List gtsForStartingNetwork, List gtsForInferNetwork, List treeCorrespondences){
        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        Map<String,Integer> exp2ID = new HashMap<String, Integer>();
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
                    Set<Integer> binaryIDs = new HashSet<Integer>();
                    for (Tree btr : Trees.getAllBinaryResolution(gtTuple.Item1)) {
                        String btrExp = Trees.getLexicographicNewickString(btr, allele2species);
                        Integer index = exp2ID.get(btrExp);
                        if (index == null) {
                            index = gtsForInferNetwork.size();
                            gtsForInferNetwork.add(btr);
                            binaryIDs.add(index);
                            exp2ID.put(btrExp, index);
                        } else {
                            binaryIDs.add(index);
                        }
                    }
                    treeCorrespondences.add(new Tuple<MutableTuple<Tree, Double>, Set<Integer>>(gtTuple, binaryIDs));

                } else {
                    existingTuple.Item2 += gtTuple.Item2;
                }

            }
        }

        gtsForStartingNetwork.addAll(exp2tree.values());

    }

}
