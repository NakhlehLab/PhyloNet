package edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by yunyu on 10/1/15.
 */
public class GibbsSamplingForPruningNetworksFromGTTPseudo_MultiTreesPerLocus extends GibbsSamplingForPruningNetworksFromGTTPseudo {

    public void summarizeGTs(List<List<MutableTuple<Tree,Double>>> originalGTs, Map<String,String> allele2species, List allTriplets, List tripletFrequencies){
        int treeID = 0;
        Map<String, Integer> distinctTree2ID = new HashMap<>();
        List<MutableTuple<Tree,Double>> distinctGTs = new ArrayList<>();
        List<List<MutableTuple<Integer,Double>>> treeCorrespondences = new ArrayList<>();
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
                    distinctGTs.add(new MutableTuple<Tree, Double>(gtTuple.Item1,gtTuple.Item2));
                    distinctTree2ID.put(exp, existingTreeID);
                    tree2infoIndex.put(exp, infoList.size());
                    infoList.add(new MutableTuple(treeID, gtTuple.Item2));
                    treeID++;
                } else {
                    distinctGTs.get(existingTreeID).Item2 += gtTuple.Item2;
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
        NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus.computeTripleFrequenciesInGTs(distinctGTs, allele2species, treeCorrespondences, allTriplets, tripletFrequencies);
    }



    protected double computeLikelihood(double[][] probs, List tripletFrequencies){
        double totalProb = 0;
        for(Object o: tripletFrequencies){
            List<Tuple<double[][],Double>> freqForOneLocus = (List<Tuple<double[][],Double>>)o;
            double probForOneLocus = 0;
            double totalWeight = 0;
            for(Tuple<double[][],Double> freqs: freqForOneLocus){
                double probForOneTree = 1;
                totalWeight += freqs.Item2;
                for(int i=0; i<probs.length; i++){
                    for(int j=0; j<3; j++){
                        if(probs[i][j] == 0)return Double.NEGATIVE_INFINITY;
                        probForOneTree *= Math.pow(probs[i][j],freqs.Item1[i][j]);
                    }
                }
                probForOneLocus += probForOneTree * freqs.Item2;
            }
            probForOneLocus /= totalWeight;
            totalProb += Math.log(probForOneLocus);

        }
        return totalProb;
    }



}
