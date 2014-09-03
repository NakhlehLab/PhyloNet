/*
 * Copyright (c) 2013 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.algos.network;

import edu.rice.cs.bioinfo.library.programming.Container;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NetworkLikelihoodFromGTTBL extends NetworkLikelihood {
    protected Map<SpeciesPair, Double> _pair2time = null;


    private void computePairwiseCoalesceTime(List<Tree> trees, Map<String,List<String>> species2alleles){
        _pair2time = new HashMap<SpeciesPair, Double>();

        Map<String,String> allele2species = null;
        if(species2alleles!=null){
            allele2species = new HashMap<String, String>();
            for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                for(String allele: entry.getValue()){
                    allele2species.put(allele, entry.getKey());
                }
            }
        }

        for(Tree tree: trees) {
            Map<TNode, Set<String>> node2leaves = new Hashtable<TNode, Set<String>>();
            Map<TNode, Double> node2height = new Hashtable<TNode, Double>();
            for (TNode node : tree.getNodes()) {
                Set<String> taxaUnder = new HashSet<String>();
                double height = 0;
                if (node.isLeaf()) {
                    if(allele2species==null){
                        taxaUnder.add(node.getName());
                    }
                    else{
                        taxaUnder.add(allele2species.get(node.getName()));
                    }
                } else {
                    Iterator children = node.getChildren().iterator();
                    TNode child1 = (TNode) children.next();
                    taxaUnder.addAll(node2leaves.get(child1));
                    TNode child2 = (TNode) children.next();
                    taxaUnder.addAll(node2leaves.get(child2));
                    height = node2height.get(child1) + child1.getParentDistance();
                    height = Math.max(height,node2height.get(child2) + child2.getParentDistance());

                    for (String taxon1 : node2leaves.get(child1)) {
                        for (String taxon2 : node2leaves.get(child2)) {
                            SpeciesPair sp = new SpeciesPair(taxon1,taxon2);
                            Double minTime = _pair2time.get(sp);
                            if(minTime == null || minTime > height){
                                _pair2time.put(sp, height);
                            }
                        }
                    }
                }
                node2leaves.put(node, taxaUnder);
                node2height.put(node, height);
            }
        }

    }



    private void initializeNetwork(Network<Object> speciesNetwork, Map<NetNode, MutableTuple<List<SpeciesPair>, Integer>> node2constraints, Map<NetNode<Object>, Double> node2height){
        Map<NetNode, Integer> node2depth = new Hashtable<NetNode, Integer>();
        Map<NetNode, Integer> node2ID = new Hashtable<NetNode, Integer>();
        int id = 0;

        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            node2ID.put(node, id++);
            if(node.isLeaf()){
                node2height.put(node, 0.0);
                node2depth.put(node, 0);
                continue;
            }
            double upperBound = -1;
            if(node2constraints.containsKey(node)){
                upperBound = node2constraints.get(node).Item1.get(0)._time;
            }
            node2height.put(node, upperBound);
            int maxDepth = 0;
            for(NetNode child: node.getChildren()){
                maxDepth = Math.max(maxDepth, node2depth.get(child));
            }
            node2depth.put(node, maxDepth+1);
        }
        boolean updated;
        do {
            updated = false;
            for(NetNode<Object> node: speciesNetwork.bfs()){
                double minParentHeight = Double.MAX_VALUE;
                for(NetNode<Object> parent: node.getParents()){
                    double parentHeight = node2height.get(parent);
                    if(parentHeight>0){
                        minParentHeight = Math.min(minParentHeight, parentHeight);
                    }
                }
                if(node2height.get(node)>minParentHeight){
                    node2height.put(node, minParentHeight);
                    updated = true;
                }
            }

        }while (updated);

        boolean[][] M = computeM(speciesNetwork, node2ID);

        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            int nodeID = node2ID.get(node);
            double minParent = Double.MAX_VALUE;
            int maxParentDepth = 0;
            double maxChild = 0;
            for(NetNode<Object> relateNode: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
                int relateNodeID = node2ID.get(relateNode);
                if(M[relateNodeID][nodeID]){
                    double parentHeight = node2height.get(relateNode);
                    if(parentHeight>=0){
                        if(minParent > parentHeight){
                            minParent = parentHeight;
                            maxParentDepth = node2depth.get(relateNode);
                        }
                        else if(minParent == parentHeight){
                            maxParentDepth = Math.max(maxParentDepth, node2depth.get(relateNode));
                        }
                    }
                }
                else if(M[nodeID][relateNodeID]){
                    double childHeight = node2height.get(relateNode);
                    if(childHeight>=0){
                        maxChild = Math.max(maxChild, childHeight);
                    }
                    else{
                        throw new RuntimeException();
                    }
                }
            }
            double currentHeight = node2height.get(node);
            if(currentHeight >= minParent || (currentHeight==-1 && minParent!=Double.MAX_VALUE)){
                int depthDiff = maxParentDepth - node2depth.get(node) + 1;
                currentHeight = maxChild + (minParent - maxChild)/depthDiff;
                //currentHeight = Math.round((maxChild + (minParent - maxChild)/depthDiff)*1000000)/1000000.0;
                node2height.put(node, currentHeight);
            }
            else if(currentHeight==-1 && minParent==Double.MAX_VALUE){
                currentHeight = maxChild + 1;
                node2height.put(node, currentHeight);
            }
        }

        double overallMin = 0;


        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            if(node.isLeaf())continue;
            double updatedHeight = node2height.get(node)-overallMin;
            double maxChild = 0;
            for(NetNode child: node.getChildren()){
                maxChild = Math.max(maxChild, node2height.get(child));
            }
            if(updatedHeight == maxChild){
                updatedHeight = maxChild + overallMin;
            }
            node2height.put(node, updatedHeight);
            for(NetNode child: node.getChildren()){
                child.setParentDistance(node, updatedHeight - node2height.get(child));
                if(child.isNetworkNode()){
                    child.setParentProbability(node,0.5);
                }
            }
        }


        //System.out.println(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.network2string(speciesNetwork));

        for(NetNode<Object> node: speciesNetwork.bfs()){
            double height = node2height.get(node);
            if(height<0){
                throw new RuntimeException();
            }
            for(NetNode child: node.getChildren()){
                if(height < node2height.get(child)){
                    throw new RuntimeException();
                }
            }
        }

    }


    protected double findOptimalBranchLength(final Network<Object> speciesNetwork, final Map<String, List<String>> species2alleles, final List gts, final List gtCorrespondence){
        boolean continueRounds = true; // keep trying to improve network

        if(_pair2time == null){
            computePairwiseCoalesceTime(gts, species2alleles);
        }

        final Map<NetNode, MutableTuple<List<SpeciesPair>, Integer>> node2constraints = new Hashtable<NetNode, MutableTuple<List<SpeciesPair>, Integer>>();
        final Map<SpeciesPair, MutableTuple<List<NetNode>, BitSet>> pairHeight2nodes = new Hashtable<SpeciesPair, MutableTuple<List<NetNode>, BitSet>>();
        computeNodeHeightUpperbound(speciesNetwork, node2constraints, pairHeight2nodes);

        final Map<NetNode<Object>, Double> node2height = new Hashtable<NetNode<Object>, Double>();
        initializeNetwork(speciesNetwork, node2constraints, node2height);
        //System.out.println("\n"+network2String(speciesNetwork));
        double initialProb = computeProbability(speciesNetwork, gts, species2alleles, gtCorrespondence);

        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(initialProb);  // records the GTProb of the network at all times
        final Container<Map<NetNode<Object>, Double>> node2heightContainer = new Container<Map<NetNode<Object>, Double>>(node2height);


        int roundIndex = 0;
        for(; roundIndex <_maxRounds && continueRounds; roundIndex++)
        {
            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();
            List<Proc> assigmentActions = new ArrayList<Proc>(); // store adjustment commands here.  Will execute them one by one later.


            for(final NetNode<Object> child : speciesNetwork.getNetworkNodes()) // find every hybrid node
            {


                Iterator<NetNode<Object>> hybridParents = child.getParents().iterator();
                final NetNode hybridParent1 = hybridParents.next();
                final NetNode hybridParent2 = hybridParents.next();

                assigmentActions.add(new Proc()
                {
                    public void execute()
                    {
                        UnivariateFunction functionToOptimize = new UnivariateFunction() {
                            public double value(double suggestedProb) {
                                double incumbentHybridProbParent1 = child.getParentProbability(hybridParent1);
                                child.setParentProbability(hybridParent1, suggestedProb);
                                child.setParentProbability(hybridParent2, 1.0 - suggestedProb);

                                double lnProb = computeProbability(speciesNetwork, gts, species2alleles, gtCorrespondence);
                                if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // change improved GTProb, keep it
                                {

                                    lnGtProbOfSpeciesNetwork.setContents(lnProb);
                                }
                                else // change did not improve, roll back
                                {
                                    child.setParentProbability(hybridParent1, incumbentHybridProbParent1);
                                    child.setParentProbability(hybridParent2, 1.0 - incumbentHybridProbParent1);
                                }
                                return lnProb;
                            }
                        };
                        BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                        try
                        {
                            optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, 0, 1.0);
                        }
                        catch(TooManyEvaluationsException e)  // _maxAssigmentAttemptsPerBranchParam exceeded
                        {
                        }

                    }
                });
            }


            for(final NetNode<Object> node : edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork))
            {
                if(node.isLeaf()){
                    continue;
                }

                assigmentActions.add(new Proc()
                {
                    public void execute()
                    {
                        final Container<Double> minHeight = new Container<Double>(0.0);
                        final Container<Double> maxHeight = new Container<Double>(Double.MAX_VALUE);


                        for(NetNode<Object> child: node.getChildren()){
                            double childHeight = node2heightContainer.getContents().get(child);
                            minHeight.setContents(Math.max(minHeight.getContents(), childHeight));
                        }

                        //TODO
                        MutableTuple<List<SpeciesPair>, Integer> spTuple = node2constraints.get(node);
                        int currentIndex = -1;
                        int maxIndex = -1;
                        boolean hasUpperBound = false;
                        boolean canLower = false;
                        if(spTuple != null){
                            currentIndex = spTuple.Item2;
                            canLower = currentIndex!=0;
                            maxIndex = spTuple.Item2;
                            SpeciesPair sp;
                            boolean canHigher = true;
                            while(canHigher && spTuple.Item1.size()>maxIndex){
                                sp = spTuple.Item1.get(maxIndex);
                                canHigher = pairHeight2nodes.get(sp).Item2.cardinality()>1;
                                if(canHigher){
                                    maxIndex++;
                                }
                            }
                            if(spTuple.Item1.size()>maxIndex){
                                hasUpperBound = true;
                            }
                        }
                        double minParent = Double.MAX_VALUE;
                        if(!node.isRoot()){
                            for(NetNode<Object> parent: node.getParents()){
                                double parentHeight = node2heightContainer.getContents().get(parent);
                                minParent = Math.min(minParent, parentHeight);
                            }
                        }
                        else{
                            minParent = minHeight.getContents() + _maxBranchLength;
                        }

                        if(hasUpperBound){
                            maxHeight.setContents(Math.min(minParent, spTuple.Item1.get(maxIndex)._time));
                        }
                        else{
                            maxHeight.setContents(minParent);
                        }
                        if(canLower){
                            canLower = spTuple.Item1.get(currentIndex-1)._time>minHeight.getContents();
                        }

                        //System.out.println("\nChanging node " + node.getName() + " from " + minHeight.getContents() + " to " + maxHeight.getContents());
                        UnivariateFunction functionToOptimize = new UnivariateFunction() {

                            public double value(double suggestedHeight) {  // brent suggests a new branch length
                                double incumbentHeight = node2heightContainer.getContents().get(node);


                                for(NetNode<Object> child: node.getChildren()){
                                    child.setParentDistance(node, suggestedHeight - node2heightContainer.getContents().get(child));
                                }

                                if(!node.isRoot()){
                                    for(NetNode<Object> parent: node.getParents()){
                                        node.setParentDistance(parent, node2heightContainer.getContents().get(parent) - suggestedHeight);
                                    }
                                }


                                double lnProb = computeProbability(speciesNetwork, gts, species2alleles, gtCorrespondence);

                                //System.out.print("suggest: "+ suggestedHeight + " " + lnProb + " vs. " + lnGtProbOfSpeciesNetwork.getContents() + ": ");
                                if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                                {
                                    lnGtProbOfSpeciesNetwork.setContents(lnProb);
                                    node2heightContainer.getContents().put(node, suggestedHeight);
                                    //System.out.println( " better ");

                                }
                                else  // didn't improve, roll back change
                                {
                                    for(NetNode<Object> child: node.getChildren()){
                                        child.setParentDistance(node, incumbentHeight - node2heightContainer.getContents().get(child));
                                    }
                                    if(!node.isRoot()){
                                        for(NetNode<Object> parent: node.getParents()){
                                            node.setParentDistance(parent, node2heightContainer.getContents().get(parent) - incumbentHeight);
                                        }
                                    }
                                    //System.out.println( " worse ");
                                }
                                return lnProb;
                            }
                        };
                        BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                        try
                        {
                            optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, minHeight.getContents(), maxHeight.getContents());
                        }
                        catch(TooManyEvaluationsException e) // _maxAssigmentAttemptsPerBranchParam exceeded
                        {
                        }

                        if(currentIndex!=maxIndex || canLower){
                            double updatedHeight = node2height.get(node);
                            int updatedIndex = 0;
                            for(; updatedIndex<spTuple.Item1.size(); updatedIndex++){
                                if(updatedHeight < spTuple.Item1.get(updatedIndex)._time){
                                    break;
                                }
                            }
                            if(updatedIndex!=currentIndex){
                                if(updatedIndex>currentIndex){
                                    for(int i=currentIndex; i<updatedIndex; i++){
                                        SpeciesPair sp = spTuple.Item1.get(i);
                                        MutableTuple<List<NetNode>, BitSet> changedSP = pairHeight2nodes.get(sp);
                                        int offBit = changedSP.Item1.indexOf(node);
                                        changedSP.Item2.set(offBit, false);
                                    }
                                }
                                else if(updatedIndex<currentIndex){
                                    currentIndex = Math.min(currentIndex, spTuple.Item1.size()-1);
                                    for(int i=currentIndex; i>=updatedIndex; i--){
                                        SpeciesPair sp = spTuple.Item1.get(i);
                                        MutableTuple<List<NetNode>, BitSet> changedSP = pairHeight2nodes.get(sp);
                                        int offBit = changedSP.Item1.indexOf(node);
                                        changedSP.Item2.set(offBit, true);
                                    }
                                }
                                spTuple.Item2 = updatedIndex;
                            }
                        }
                        //System.out.println(network2String(speciesNetwork) + " " + lnGtProbOfSpeciesNetwork.getContents());
                    }

                });
            }

            Collections.shuffle(assigmentActions);

            for(Proc assigment : assigmentActions)   // for each change attempt, perform attempt
            {
                assigment.execute();
            }

            if( ((double)lnGtProbOfSpeciesNetwork.getContents()) == lnGtProbLastRound)  // if no improvement was made wrt to last around, stop trying to find a better assignment
            {
                continueRounds = false;
            }
            else if (lnGtProbOfSpeciesNetwork.getContents() > lnGtProbLastRound) // improvement was made, ensure it is large enough wrt to improvement threshold to continue searching
            {

                double improvementPercentage = Math.pow(Math.E, (lnGtProbOfSpeciesNetwork.getContents() - lnGtProbLastRound)) - 1.0;  // how much did we improve over last round
                //System.out.println(improvementPercentage + " vs. " + _improvementThreshold);
                if(improvementPercentage < _improvementThreshold  )  // improved, but not enough to keep searching
                {
                    continueRounds = false;
                }
            }
            else
            {
                throw new IllegalStateException("Should never have decreased prob.");
            }
        }
        return lnGtProbOfSpeciesNetwork.getContents();
    }



    private void computeNodeHeightUpperbound(Network network, Map<NetNode, MutableTuple<List<SpeciesPair>, Integer>> node2constraints, Map<SpeciesPair, MutableTuple<List<NetNode>, BitSet>> pairHeight2nodes){
        Map<NetNode, Set<String>> node2leaves = new Hashtable<NetNode, Set<String>>();
        for(Object nodeObject: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeObject;
            Set<String> leafSet = new HashSet<String>();
            if(node.isLeaf()){
                leafSet.add(node.getName());
            }
            else{
                Iterator<NetNode> childIt = node.getChildren().iterator();
                if(node.isNetworkNode()){
                    leafSet.addAll(node2leaves.get(childIt.next()));
                }
                else{
                    NetNode child1 = childIt.next();
                    Set<String> child1Leaves = new HashSet<String>(node2leaves.get(child1));
                    leafSet.addAll(child1Leaves);
                    NetNode child2 = childIt.next();
                    Set<String> child2Leaves = new HashSet<String>(node2leaves.get(child2));
                    leafSet.addAll(child2Leaves);
                    Set<String> temp = new HashSet<String>(child1Leaves);
                    child1Leaves.removeAll(child2Leaves);
                    child2Leaves.removeAll(temp);
                    if(child1Leaves.size()!=0 && child2Leaves.size()!=0){
                        List<SpeciesPair> spList = new ArrayList<SpeciesPair>();
                        for(String leaf1: child1Leaves){
                            for(String leaf2: child2Leaves){
                                SpeciesPair sp = new SpeciesPair(leaf1, leaf2);
                                double time = _pair2time.get(sp);
                                sp.setTime(time);
                                int index = 0;
                                for(SpeciesPair exsp: spList){
                                    if(time < exsp._time){
                                        break;
                                    }
                                    index++;
                                }
                                spList.add(index, sp);
                                MutableTuple<List<NetNode>, BitSet> MutableTuple = pairHeight2nodes.get(sp);
                                if(MutableTuple == null){
                                    MutableTuple = new MutableTuple<List<NetNode>, BitSet>(new ArrayList<NetNode>(),null);
                                    pairHeight2nodes.put(sp, MutableTuple);
                                }
                                MutableTuple.Item1.add(node);
                            }
                        }
                        MutableTuple<List<SpeciesPair>, Integer> MutableTuple = new MutableTuple<List<SpeciesPair>, Integer>(spList, 0);
                        node2constraints.put(node, MutableTuple);
                    }
                }
            }
            node2leaves.put(node, leafSet);
        }
        for(MutableTuple<List<NetNode>, BitSet> tuple: pairHeight2nodes.values()){
            BitSet bs = new BitSet(tuple.Item1.size());
            bs.set(0, tuple.Item1.size(), true);
            tuple.Item2 = bs;
        }
    }

    private boolean[][] computeM(Network<Object> net, Map<NetNode, Integer> node2ID){
        int numNodes = node2ID.size();
        boolean[][] M = new boolean[numNodes][numNodes];
        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(net)){
            int pID = node2ID.get(node);
            //M[pID][pID] = false;
            for(NetNode child: node.getChildren()){
                int cID = node2ID.get(child);
                M[pID][cID] = true;
                for(int i=0; i<numNodes; i++){
                    if(M[cID][i]){
                        M[pID][i] = true;
                    }
                }
            }
        }
        return M;
    }



    private class MyThreadFromScratch extends Thread{
        GeneTreeWithBranchLengthProbabilityYF _gtp;
        double[] _probs;


        public MyThreadFromScratch(GeneTreeWithBranchLengthProbabilityYF gtp, double[] probs){
            _probs = probs;
            _gtp = gtp;
        }


        public void run() {
            _gtp.calculateGTDistribution(_probs);

        }
    }


    protected double computeProbability(Network<Object> speciesNetwork, List geneTrees, Map<String, List<String>> species2alleles, List gtCorrespondences) {
        double[] probArray = new double[geneTrees.size()];
        Thread[] myThreads = new Thread[_numThreads];

        GeneTreeWithBranchLengthProbabilityYF gtp = new GeneTreeWithBranchLengthProbabilityYF(speciesNetwork, geneTrees, species2alleles);
        gtp.setParallel(true);


        for(int i=0; i<_numThreads; i++){
            myThreads[i] = new MyThreadFromScratch(gtp, probArray);
            myThreads[i].start();
        }

        for(int i=0; i<_numThreads; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }

        return calculateFinalLikelihood(probArray, gtCorrespondences);
    }




    abstract protected double calculateFinalLikelihood(double[] probs, List gtCorrespondences);



    class SpeciesPair{
        String _species1;
        String _species2;
        double _time;

        public SpeciesPair(String s1, String s2){
            _species1 = s1;
            _species2 = s2;
        }

        public SpeciesPair(String s1, String s2, double t){
            _species1 = s1;
            _species2 = s2;
            _time = t;
        }

        public void setTime(double t){
            _time = t;
        }

        public int hashCode(){
            return _species1.hashCode() * _species2.hashCode();
        }

        public boolean equals(Object o) {
            if(!(o instanceof SpeciesPair)){
                return false;
            }
            SpeciesPair sp = (SpeciesPair)o;
            if((this._species1.equals(sp._species1) && this._species2.equals(sp._species2)) || (this._species1.equals(sp._species2) && this._species2.equals(sp._species1))){
                return true;
            }
            else{
                return false;
            }
        }

        public String toString(){
            return _species1+"|"+_species2 + ":" + _time;
        }
    }
}
