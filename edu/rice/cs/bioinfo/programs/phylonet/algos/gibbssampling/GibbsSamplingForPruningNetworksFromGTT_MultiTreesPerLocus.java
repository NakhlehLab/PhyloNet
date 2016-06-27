package edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by yunyu on 10/1/15.
 */
public class GibbsSamplingForPruningNetworksFromGTT_MultiTreesPerLocus extends GibbsSamplingForPruningNetworksFromGTT {
    protected void summarizeGTs(List<List<MutableTuple<Tree,Double>>> originalGTs, Map<String, String> allele2species, List distinctGTs, List treeCorrespondences){
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
                    existingTuple = new MutableTuple<>(gtTuple.Item1, gtTuple.Item2);
                    exp2tree.put(exp, existingTuple);
                    Set<Integer> binaryIDs = new HashSet<Integer>();
                    for (Tree btr : Trees.getAllBinaryResolution(gtTuple.Item1)) {
                        String btrExp = Trees.getLexicographicNewickString(btr, allele2species);
                        Integer index = exp2ID.get(btrExp);
                        if (index == null) {
                            index = distinctGTs.size();
                            distinctGTs.add(btr);
                            binaryIDs.add(index);
                            exp2ID.put(btrExp, index);
                        } else {
                            binaryIDs.add(index);
                        }
                    }
                    treeCorrespondences.add(new Tuple<MutableTuple<Tree, Double>, Set<Integer>>(existingTuple, binaryIDs));
                } else {
                    existingTuple.Item2 += gtTuple.Item2;
                }
            }
        }
    }


    protected double computeLikelihood(double[] probList, List gtCorrespondences){
        double totalProb = 0;
        for(Object o: gtCorrespondences){
            List<MutableTuple<Integer,Double>> locusInfo = (List<MutableTuple<Integer,Double>>)o;
            double totalProbForLocus = 0;
            double totalWeight = 0;
            for(MutableTuple<Integer,Double> one: locusInfo){
                totalProbForLocus += probList[one.Item1] * one.Item2;
                totalWeight += one.Item2;
            }
            totalProbForLocus = totalProbForLocus/totalWeight;
            totalProb += Math.log(totalProbForLocus);
        }
        return totalProb;
    }




}
