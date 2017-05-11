package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by dw20 on 5/11/17.
 */
public class GTTPseudoLikelihood_MultiPerLocus extends NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus {

    public double computeProbability(Network<Object> speciesNetwork, List distinctTrees,
                                     Map<String,List<String>> species2alleles, List gtCorrespondences) {
        return super.computeProbability(speciesNetwork, distinctTrees, gtCorrespondences, species2alleles);
    }


    protected void summarizeData(List originalGTs, Map<String,String> allele2species,
                                 List gtsForStartingNetwork, List allTriplets, List tripletFrequencies){
        int treeID = 0;
        Map<String, Integer> distinctTree2ID = new HashMap<>();
        List<List<MutableTuple<Integer,Double>>> treeCorrespondences = new ArrayList<>();
        int i=0;
        for(Object o: originalGTs) {
            List<MutableTuple<Tree,Double>> treesForOneLocus = (List<MutableTuple<Tree,Double>>)o;
            Map<String, Integer> tree2infoIndex = new HashMap<String, Integer>();
            List<MutableTuple<Integer,Double>> infoList = new ArrayList<MutableTuple<Integer, Double>>();
            for (MutableTuple<Tree, Double> gtTuple : treesForOneLocus) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                Integer existingTreeID = distinctTree2ID.get(exp);
                if (existingTreeID == null) {
                    existingTreeID = treeID;
                    gtsForStartingNetwork.add(new MutableTuple<Tree, Double>(gtTuple.Item1,gtTuple.Item2));
                    distinctTree2ID.put(exp, existingTreeID);
                    tree2infoIndex.put(exp, infoList.size());
                    infoList.add(new MutableTuple(treeID, gtTuple.Item2));
                    treeID++;
                } else {
                    ((MutableTuple<Tree,Double>)(gtsForStartingNetwork.get(existingTreeID))).Item2 += gtTuple.Item2;
                    Integer infoID = tree2infoIndex.get(exp);
                    if(infoID == null){
                        tree2infoIndex.put(exp, infoList.size());
                        infoList.add(new MutableTuple(existingTreeID, gtTuple.Item2));
                    }
                    else{
                        infoList.get(infoID).Item2 += gtTuple.Item2;
                    }
                }
            }
            treeCorrespondences.add(infoList);
        }
        NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus.computeTripleFrequenciesInGTs(
                gtsForStartingNetwork, allele2species, treeCorrespondences, allTriplets, tripletFrequencies);
    }

}
