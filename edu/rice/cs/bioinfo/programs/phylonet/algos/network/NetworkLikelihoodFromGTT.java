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
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;

import java.util.*;

/**
 * Created by Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class inherits from NetworkLikelihood.
 * It calculates the likelihood of a species network when the input is a collection of gene trees (only topologies are used)
 *
 * See "Maximum Likelihood Inference of Reticulate Evolutionary Histories", Proceedings of the National Academy of Sciences, 2014
 */
public abstract class NetworkLikelihoodFromGTT extends NetworkLikelihood {


    /**
     * This function is to optimize the branch lengths and inheritance probabilities of a given species network
     *
     * @param speciesNetwork            the species network
     * @param species2alleles           mapping from species to alleles which they is sampled from
     * @param distinctTrees             summarized data
     * @param singleAlleleSpecies       to help identify which branch lengths in the species network can be ignored
     * @param gtCorrespondence          relationships between the original data and the data in dataForInferNetwork
     *
     * @return likelihood of the species network after its branch lengths and inheritance probabilities are optimized
     */
    protected double findOptimalBranchLength(final Network<Object> speciesNetwork, final Map<String, List<String>> species2alleles, final List distinctTrees, final List gtCorrespondence, final Set<String> singleAlleleSpecies){
        boolean continueRounds = true; // keep trying to improve network
        for(NetNode<Object> node: speciesNetwork.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent,1.0);
                if(node.isNetworkNode()){
                    node.setParentProbability(parent, 0.5);
                }
            }
        }

        Set<NetNode> node2ignoreForBL = findEdgeHavingNoBL(speciesNetwork, singleAlleleSpecies);
        double initalProb = computeProbabilityForCached(speciesNetwork, distinctTrees, species2alleles, gtCorrespondence);
        if(_printDetails)
            System.out.println(speciesNetwork.toString() + " : " + initalProb);

        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(initalProb);  // records the GTProb of the network at all times

        int roundIndex = 0;
        for(; roundIndex <_maxRounds && continueRounds; roundIndex++)
        {
             /*
             * Prepare a random ordering of network edge examinations each of which attempts to change a branch length or hybrid prob to improve the GTProb score.
             */
            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();
            List<Proc> assigmentActions = new ArrayList<Proc>(); // store adjustment commands here.  Will execute them one by one later.


            for(final NetNode<Object> parent : edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork))
            {

                for(final NetNode<Object> child : parent.getChildren())
                {
                    if(node2ignoreForBL.contains(child)){
                        continue;
                    }

                    assigmentActions.add(new Proc()
                    {
                        public void execute()
                        {

                            UnivariateFunction functionToOptimize = new UnivariateFunction() {
                                public double value(double suggestedBranchLength) {
                                    double incumbentBranchLength = child.getParentDistance(parent);

                                    child.setParentDistance(parent, suggestedBranchLength);
                                    double lnProb = updateProbabilityForCached(speciesNetwork, distinctTrees, gtCorrespondence, child, parent);
                                    if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                                    {
                                        lnGtProbOfSpeciesNetwork.setContents(lnProb);

                                    }
                                    else  // didn't improve, roll back change
                                    {
                                        child.setParentDistance(parent, incumbentBranchLength);
                                    }
                                    return lnProb;
                                }
                            };
                            BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                            try
                            {
                                optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, Double.MIN_VALUE, _maxBranchLength);
                            }
                            catch(TooManyEvaluationsException e) // _maxAssigmentAttemptsPerBranchParam exceeded
                            {
                            }

                            updateProbabilityForCached(speciesNetwork, distinctTrees, gtCorrespondence, child, parent);
                            if(_printDetails)
                                System.out.println(speciesNetwork.toString() + " : " + lnGtProbOfSpeciesNetwork.getContents());

                        }
                    });
                }
            }


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

                                double lnProb = updateProbabilityForCached(speciesNetwork, distinctTrees, gtCorrespondence, child, null);
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
                            if(child.getName().equals("Y"))optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, 0.6, 0.8);
                            else
                            optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, 0, 1.0);
                        }
                        catch(TooManyEvaluationsException e)  // _maxAssigmentAttemptsPerBranchParam exceeded
                        {
                        }
                        updateProbabilityForCached(speciesNetwork, distinctTrees, gtCorrespondence, child, null);
                        if(_printDetails)
                            System.out.println(speciesNetwork.toString() + " : " + lnGtProbOfSpeciesNetwork.getContents());
                    }
                });
            }

            Collections.shuffle(assigmentActions);

            for(Proc assigment : assigmentActions)   // for each change attempt, perform attempt
            {
                assigment.execute();
            }
            if(_printDetails) {
                System.out.println("Round end ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                System.out.println(speciesNetwork.toString() + "\n" + lnGtProbOfSpeciesNetwork.getContents() + "\n");
            }
            if( ((double)lnGtProbOfSpeciesNetwork.getContents()) == lnGtProbLastRound)  // if no improvement was made wrt to last around, stop trying to find a better assignment
            {
                continueRounds = false;
            }
            else if (lnGtProbOfSpeciesNetwork.getContents() > lnGtProbLastRound) // improvement was made, ensure it is large enough wrt to improvement threshold to continue searching
            {

                double improvementPercentage = Math.pow(Math.E, (lnGtProbOfSpeciesNetwork.getContents() - lnGtProbLastRound)) - 1.0;  // how much did we improve over last round
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



    /**
     * This class is for using parallel to compute the likelihood when the intermediate ancestral configurations are not stored
     * It is slow for updating probabilities afterwards but needs less memory
     * */
    private class MyThreadNonCached extends Thread{
        GeneTreeProbabilityYF _gtp;
        Network _speciesNetwork;
        List<Tree> _geneTrees;
        Map<String, List<String>> _species2alleles;
        double[] _probs;


        public MyThreadNonCached(GeneTreeProbabilityYF gtp, Network speciesNetwork, List<Tree> geneTrees, Map<String, List<String>> species2alleles, double[] probs){
            _speciesNetwork = speciesNetwork;
            _geneTrees = geneTrees;
            _species2alleles = species2alleles;
            _probs = probs;
            _gtp = gtp;
        }


        public void run() {
            _gtp.calculateGTDistribution(_speciesNetwork, _geneTrees, _species2alleles, _probs);

        }
    }


    /**
     * This class is for using parallel to compute the likelihood from scratch when the intermediate ancestral configurations are stored
     * It is fast for updating probabilities afterwards but needs more memory
     * */
    private class MyThreadFromScratchForCached extends Thread{
        GeneTreeProbabilityYF_Cached _gtp;
        Network _speciesNetwork;
        List<Tree> _geneTrees;
        Map<String, List<String>> _species2alleles;
        double[] _probs;


        public MyThreadFromScratchForCached(GeneTreeProbabilityYF_Cached gtp, Network speciesNetwork, List<Tree> geneTrees, Map<String, List<String>> species2alleles, double[] probs){
            _speciesNetwork = speciesNetwork;
            _geneTrees = geneTrees;
            _species2alleles = species2alleles;
            _probs = probs;
            _gtp = gtp;
        }


        public void run() {
            _gtp.calculateGTDistribution(_speciesNetwork, _geneTrees, _species2alleles, _probs);

        }
    }



    /**
     * This class is for using parallel to update the likelihood using stored intermediate ancestral configurations
     * */
    private class MyThreadFromNonScratchForCached extends Thread{
        Network _speciesNetwork;
        List<Tree> _gts;
        double[] _probs;
        Set<NetNode> _childNodes;
        Set<NetNode> _parentNodes;
        GeneTreeProbabilityYF_Cached _gtp;


        public MyThreadFromNonScratchForCached(GeneTreeProbabilityYF_Cached gtp, Network speciesNetwork, List<Tree> gts, Set<NetNode> childNodes, Set<NetNode> parentNodes, double[] probs){
            _speciesNetwork = speciesNetwork;
            _gts = gts;
            _probs = probs;
            _childNodes = childNodes;
            _parentNodes = parentNodes;
            _gtp = gtp;
        }


        public void run() {
            _gtp.calculateGTDistribution(_speciesNetwork, _gts, _childNodes, _parentNodes, _probs);

        }
    }


    /**
     * This function is to use single thread to compute the likelihood without caching intermediate ancestral configurations
     *
     * @param speciesNetwork        the species network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param distinctTrees         summarized gene trees
     * @param gtCorrespondences     the correspondences between the summarized gene trees and the original gene trees
     *
     * @return likelihood
     */
    protected double computeProbability(Network<Object> speciesNetwork, List distinctTrees, List gtCorrespondences, Map<String, List<String>> species2alleles) {
        double[] probs = new double[distinctTrees.size()];
        Thread[] myThreads = new Thread[_numThreads];

        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        gtp.setParallel(true);
        gtp.preProcess(speciesNetwork, distinctTrees, true);


        for(int i=0; i<_numThreads; i++){
            myThreads[i] = new MyThreadNonCached(gtp, speciesNetwork, distinctTrees, species2alleles, probs);
            myThreads[i].start();
        }

        for(int i=0; i<_numThreads; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }

        double prob = calculateFinalLikelihood(probs, gtCorrespondences);
        return prob;
    }



    /**
     * This function is to use single thread to compute the likelihood from scratch when caching intermediate ancestral configurations
     *
     * @param speciesNetwork        the species network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param distinctTrees         summarized gene trees
     * @param gtCorrespondences     the correspondences between the summarized gene trees and the original gene trees
     *
     * @return likelihood
     */
    protected double computeProbabilityForCached(Network<Object> speciesNetwork, List distinctTrees, Map<String, List<String>> species2alleles, List gtCorrespondences) {
        double[] probs = new double[distinctTrees.size()];
        Thread[] myThreads = new Thread[_numThreads];

        GeneTreeProbabilityYF_Cached gtp = new GeneTreeProbabilityYF_Cached();
        gtp.setParallel(true);
        gtp.preProcess(speciesNetwork, distinctTrees, true);


        for(int i=0; i<_numThreads; i++){
            myThreads[i] = new MyThreadFromScratchForCached(gtp, speciesNetwork, distinctTrees, species2alleles, probs);
            myThreads[i].start();
        }

        for(int i=0; i<_numThreads; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }

        return calculateFinalLikelihood(probs, gtCorrespondences);
    }



    /**
     * This function is to use single thread to update the likelihood using cached intermediate ancestral configurations
     *
     * @param speciesNetwork        the species network
     * @param geneTrees             summarized gene trees
     * @param gtCorrespondences     the correspondences between the summarized gene trees and the original gene trees
     * @param child                 the tail of the branch whose length or inheritance probability has been changed
     * @param parent                the head of the branch whose length or inheritance probability has been changed
     *
     * @return likelihood
     */
    private double updateProbabilityForCached(Network speciesNetwork, List<Tree> geneTrees, final List gtCorrespondences, NetNode child, NetNode parent) {
        Set<NetNode> childNodes = new HashSet<NetNode>();
        childNodes.add(child);
        Set<NetNode> parentNodes = new HashSet<NetNode>();
        if(parent==null){
            for(Object parentNode: child.getParents()){
                parentNodes.add((NetNode)parentNode);
            }
        }
        else{
            parentNodes.add(parent);
        }

        double[] probs = new double[geneTrees.size()];
        Thread[] myThreads = new Thread[_numThreads];

        GeneTreeProbabilityYF_Cached gtp = new GeneTreeProbabilityYF_Cached();
        gtp.setParallel(true);
        gtp.preProcess(speciesNetwork, geneTrees, false);

        for(int i=0; i<_numThreads; i++){
            myThreads[i] = new MyThreadFromNonScratchForCached(gtp, speciesNetwork, geneTrees, childNodes, parentNodes, probs);
            myThreads[i].start();
        }

        for(int i=0; i<_numThreads; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }

        double probability = calculateFinalLikelihood(probs, gtCorrespondences);
        return probability;
    }



    /**
     * This function is to help find the set of branches whose lengths cannot be estimated so that they can be ignored during the inference
     *
     * @param speciesNetwork        the species network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param singleAlleleSpecies   species that have only one allele sampled from it
     */
    protected void findSingleAlleleSpeciesSet(Network speciesNetwork, Map<String,List<String>> species2alleles, Set<String> singleAlleleSpecies){
        for(Object node: speciesNetwork.getLeaves()){
            String species = ((NetNode)node).getName();
            if(species2alleles == null || species2alleles.get(species).size() == 1)
                singleAlleleSpecies.add(species);
        }
    }


    /**
     * This function is to find the set of branches whose lengths cannot be estimated so that they can be ignored during the inference
     *
     * @param network               the species network
     * @param singleAlleleSpecies   species that have only one allele sampled from it
     */
    private Set<NetNode> findEdgeHavingNoBL(Network network, Set<String> singleAlleleSpecies) {
        Set<NetNode> node2ignore = new HashSet<>();
        Map<NetNode, Set<String>> node2leaves = new HashMap<>();
        for (Object nodeO : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeO;
            Set<String> leaves = new HashSet<>();
            if (node.isLeaf()) {
                leaves.add(node.getName());
            }
            else {
                for (Object childO : node.getChildren()) {
                    NetNode childNode = (NetNode) childO;
                    Set<String> childLeaves = node2leaves.get(childNode);
                    leaves.addAll(childLeaves);
                }
            }
            if(leaves.size()<=1 && singleAlleleSpecies.containsAll(leaves)){
                node2ignore.add(node);
            }
            node2leaves.put(node, leaves);

        }
        return node2ignore;
    }




    /**
     * This function is to calculate the final log likelihood using the correspondences between the summarized data and the original data
     *
     * @param probs               the probabilities of each summarized data respectively
     * @param dataCorrespondences   the correspondences between the summarized data and the original data
     */
    abstract protected double calculateFinalLikelihood(double[] probs, List dataCorrespondences);


}
