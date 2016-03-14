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

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
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
    protected Map<UnorderedPair, Double> _pair2time = null;


    private void computePairwiseCoalesceTime(List<Tree> trees, Map<String,List<String>> species2alleles){
        _pair2time = new HashMap<UnorderedPair, Double>();

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
                            UnorderedPair sp = new UnorderedPair(taxon1,taxon2);
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



    private void initializeNetwork(Network<Object> speciesNetwork, Map<NetNode, Double> node2constraints, Map<NetNode<Object>, Double> node2height){
        Map<NetNode, Integer> node2depth = new Hashtable<NetNode, Integer>();
        Map<NetNode, Integer> node2ID = new Hashtable<NetNode, Integer>();
        int id = 0;

        for(NetNode<Object> node: Networks.postTraversal(speciesNetwork)){
            node2ID.put(node, id++);
            if(node.isLeaf()){
                node2height.put(node, 0.0);
                node2depth.put(node, 0);
                continue;
            }
            double upperBound = -1;
            if(node2constraints.get(node)!=Double.POSITIVE_INFINITY){
                upperBound = node2constraints.get(node);
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

        for(NetNode<Object> node: Networks.postTraversal(speciesNetwork)){
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


        //System.out.println(speciesNetwork);

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


    protected double findOptimalBranchLength(final Network<Object> speciesNetwork, final Map<String, List<String>> species2alleles, final List gts, final List gtCorrespondence, final Set<String> singleAlleleSpecies){
        boolean continueRounds = true; // keep trying to improve network

        if(_pair2time == null){
            computePairwiseCoalesceTime(gts, species2alleles);
        }
        //System.out.println("\n"+speciesNetwork);
        final Map<NetNode, Double> node2constraints = new Hashtable<NetNode, Double>();
        computeNodeHeightUpperbound(speciesNetwork, node2constraints);

        final Map<NetNode<Object>, Double> node2height = new Hashtable<NetNode<Object>, Double>();
        initializeNetwork(speciesNetwork, node2constraints, node2height);

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


            for(final NetNode<Object> node : Networks.postTraversal(speciesNetwork))
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

                        maxHeight.setContents(Math.min(minParent, node2constraints.get(node)));


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

        //System.out.print("\n" + lnGtProbOfSpeciesNetwork.getContents() + ": " + speciesNetwork);
        return lnGtProbOfSpeciesNetwork.getContents();
    }



    private void computeNodeHeightUpperbound(Network network, Map<NetNode, Double> node2constraints){
        Map<NetNode, Set<String>> node2taxa = new HashMap<>();
        for(Object o: Networks.postTraversal(network)){
            NetNode node = (NetNode)o;
            Set<String> taxa = new HashSet<>();
            double upperBound = Double.POSITIVE_INFINITY;
            if(node.isLeaf()){
                taxa.add(node.getName());
            }
            else if(node.isNetworkNode()){
                NetNode childNode = (NetNode)node.getChildren().iterator().next();
                if(!childNode.isLeaf())
                    upperBound = node2constraints.get(childNode);
                taxa.addAll(node2taxa.get(childNode));
            }
            else {
                Set<String> intersection = null;
                List<NetNode> childNodes = null;
                for (Object childO : node.getChildren()) {
                    NetNode childNode = (NetNode) childO;
                    if (childNodes == null) {
                        childNodes = new ArrayList<>();
                    }
                    childNodes.add(childNode);
                    if (intersection == null) {
                        intersection = new HashSet<>();
                        intersection.addAll(node2taxa.get(childNode));
                    } else {
                        intersection.retainAll(node2taxa.get(childNode));
                    }

                    taxa.addAll(node2taxa.get(childNode));
                }

                for (int i = 0; i < childNodes.size(); i++) {
                    Set<String> taxa1 = node2taxa.get(childNodes.get(i));
                    for (int j = i + 1; j < childNodes.size(); j++) {
                        Set<String> taxa2 = node2taxa.get(childNodes.get(j));
                        for (String taxon1 : taxa1) {
                            if (intersection.contains(taxon1))
                                continue;
                            for (String taxon2 : taxa2) {
                                if (intersection.contains(taxon2))
                                    continue;
                                upperBound = Math.min(upperBound, _pair2time.get(new UnorderedPair(taxon1, taxon2)));
                            }
                        }

                    }
                }
            }
            if(!node.isLeaf()){
                node2constraints.put(node, upperBound);
            }
            node2taxa.put(node, taxa);
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

        GeneTreeWithBranchLengthProbabilityYF gtp = new GeneTreeWithBranchLengthProbabilityYF(speciesNetwork, geneTrees, species2alleles);
        //
        if(_numThreads==1){
            gtp.calculateGTDistribution(probArray);
        }
        else {
            Thread[] myThreads = new Thread[_numThreads];
            gtp.setParallel(true);

            for (int i = 0; i < _numThreads; i++) {
                myThreads[i] = new MyThreadFromScratch(gtp, probArray);
                myThreads[i].start();
            }

            for (int i = 0; i < _numThreads; i++) {
                try {
                    myThreads[i].join();
                } catch (InterruptedException ignore) {
                }
            }
        }

        double prob = calculateFinalLikelihood(probArray, gtCorrespondences);
        return prob;
    }

    protected void findSingleAlleleSpeciesSet(Network speciesNetwork, Map<String,List<String>> species2alleles, Set<String> singleAlleleSpecies){}


    abstract protected double calculateFinalLikelihood(double[] probs, List gtCorrespondences);



    
}
