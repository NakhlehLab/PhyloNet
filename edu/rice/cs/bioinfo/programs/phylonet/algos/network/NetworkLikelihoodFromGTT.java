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
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NetworkLikelihoodFromGTT extends NetworkLikelihood {

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
                                    //System.out.println(speciesNetwork + ": " + lnProb);
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
                                //System.out.println(speciesNetwork + ": " + lnProb);
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


            // add hybrid probs to hybrid edges
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
        //System.out.println("\n" + lnGtProbOfSpeciesNetwork.getContents() + ": " + speciesNetwork);
        return lnGtProbOfSpeciesNetwork.getContents();
    }


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


    protected double computeProbability(Network<Object> speciesNetwork, List distinctTrees, Map<String, List<String>> species2alleles, List gtCorrespondences) {
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

        _maxNumACs = gtp._maxNumACs;

        double prob = calculateFinalLikelihood(probs, gtCorrespondences);
        return prob;
    }


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

        _maxNumACs = gtp._maxNumACs;

        return calculateFinalLikelihood(probs, gtCorrespondences);
    }



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
        //System.out.println(speciesNetwork.toString() + ": " + probability);
        return probability;
    }


    protected void findSingleAlleleSpeciesSet(Network speciesNetwork, Map<String,List<String>> species2alleles, Set<String> singleAlleleSpecies){
        for(Object node: speciesNetwork.getLeaves()){
            String species = ((NetNode)node).getName();
            if(species2alleles == null || species2alleles.get(species).size() == 1)
                singleAlleleSpecies.add(species);
        }
    }

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


    abstract protected double calculateFinalLikelihood(double[] probs, List gtCorrespondences);


}
