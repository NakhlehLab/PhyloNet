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

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.RichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.phylogenetics.FindRoot;
import edu.rice.cs.bioinfo.library.phylogenetics.GetDirectSuccessors;
import edu.rice.cs.bioinfo.library.phylogenetics.GetInDegree;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkNeighbourhoodRandomWalkGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkWholeNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours.AllNeighboursHillClimberFirstBetter;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours.AllNeighboursHillClimberSteepestAscent;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;

import java.io.ByteArrayInputStream;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferILSNetworkProbabilistically extends MDCOnNetworkYFFromRichNewickJung {
    protected Network[] _optimalNetworks;
    protected double[] _optimalScores;
    protected int _maxRounds;
    protected int _maxTryPerBranch;
    protected double _improvementThreshold;
    protected double _maxBranchLength;
    protected double _Brent1;
    protected double _Brent2;
    protected Long _maxExaminations;
    protected int _diameterLimit;
    protected Network _startNetwork;


    public InferILSNetworkProbabilistically(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }

    public void setSearchParameter(int maxRounds, int maxTryPerBranch, double improvementThreshold, double maxBranchLength, double Brent1, double Brent2, Long maxExaminations, int diameterLimit, Network startNetwork){
        _maxRounds = maxRounds;
        _maxTryPerBranch = maxTryPerBranch;
        _improvementThreshold = improvementThreshold;
        _maxBranchLength = maxBranchLength;
        //_maxBranchLength = 12;
        _Brent1 = Brent1;
        _Brent2 = Brent2;
        _maxExaminations = maxExaminations;
        _diameterLimit = diameterLimit;
        _startNetwork = startNetwork;
    }

    public List<Tuple<Network,Double>> inferNetwork(List<Tree> gts, Map<String,List<String>> species2alleles, int maxReticulations, int numSol){
        _optimalNetworks = new Network[numSol];
        _optimalScores = new double[numSol];
        Arrays.fill(_optimalScores, Double.NEGATIVE_INFINITY);

        List<Tree> distinctTrees = new ArrayList<Tree>();
        List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList = new ArrayList<Tuple3<Tree, Double, List<Integer>>>();
        summarizeGeneTrees(gts, distinctTrees, nbTreeAndCountAndBinaryIDList);
        DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = getStartNetwork(gts, species2alleles,_startNetwork);
        NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
        AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);

        //double[] operationProb = {0.15,0.15,0.2,0.5};
        //NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(operationProb, makeNode, makeEdge);
        //AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);

        Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> scorer = getScoreFunction(distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);
        Comparator<Double> comparator = getDoubleScoreComparator();
        //DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = getStartNetwork(gts, species2alleles,startNetwork);
        HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, _maxExaminations, maxReticulations, _diameterLimit); // search starts here
        //HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, maxExaminations, maxReticulations, new Long(100)); // search starts here

        List<Tuple<Network, Double>> resultList = new ArrayList<Tuple<Network, Double>>();
        for(int i=0; i<numSol; i++){
            resultList.add(new Tuple<Network, Double>(_optimalNetworks[i], _optimalScores[i]));
        }
        //System.out.println("\n #Networks " + result.ExaminationsCount);
        return resultList;
    }


    protected DirectedGraphToGraphAdapter<String,PhyloEdge<String>> getStartNetwork(List<Tree> gts, Map<String,List<String>> species2alleles, Network<Object> startingNetwork){
        if(startingNetwork == null){
            Map<String,String> allele2species = null;
            if(species2alleles!=null){
                allele2species = new HashMap<String, String>();
                for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                    String species = entry.getKey();
                    for(String allele: entry.getValue()){
                        allele2species.put(allele,species);
                    }
                }
            }
            MDCInference_DP mdc = new MDCInference_DP();
            Solution sol;
            if(allele2species==null){
                sol = mdc.inferSpeciesTree(gts, false, 1, false, 100, true, -1).get(0);
            }
            else{
                sol = mdc.inferSpeciesTree(gts, allele2species, false, 1, false, 100, true, -1).get(0);
            }

            Tree startingTree= Trees.generateRandomBinaryResolution(sol._st);
            startingNetwork = string2Network(startingTree.toString());
        }



        int index = 1;
        for(NetNode<Object> node: startingNetwork.dfs()){
            if(node.getName()==null || node.getName().equals("")){
                String name;
                do{
                    name = "i" + (index++);
                }while(startingNetwork.findNode(name)!=null);
                node.setName(name);
            }
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent, Double.NaN);
                node.setParentSupport(parent, Double.NaN);
                node.setParentProbability(parent, Double.NaN);
            }
        }

        String newNetwork = network2String(startingNetwork);

        return makeNetwork(newNetwork);

    }

    protected void summarizeGeneTrees(List<Tree> originalGTs, List<Tree> distinctGTs, List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        for(Tree tr: originalGTs){
            Double weight = ((STINode<Double>)tr.getRoot()).getData();
            if(weight == null){
                weight = 1.0;
            }
            int index = 0;
            Tuple3<Tree, Double, List<Integer>> newTuple = null;
            for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
                if(Trees.haveSameRootedTopology(tr, triple.Item1)){
                    newTuple = new Tuple3<Tree, Double, List<Integer>>(tr, triple.Item2 + weight, triple.Item3);
                    break;
                }
                index++;
            }
            if(newTuple!=null){
                nbTreeAndCountAndBinaryIDList.set(index, newTuple);
            }
            else{
                List<Integer> binaryIDs = new ArrayList<Integer>();
                for(Tree btr: Trees.getAllBinaryResolution(tr)){
                    index = 0;
                    boolean exist = false;
                    for(Tree exTr: distinctGTs){
                        if(Trees.haveSameRootedTopology(btr, exTr)){
                            binaryIDs.add(index);
                            exist = true;
                            break;
                        }
                        index++;
                    }

                    if(!exist){
                        distinctGTs.add(btr);
                        binaryIDs.add(index);
                    }
                }
                newTuple = new Tuple3<Tree, Double, List<Integer>>(tr, weight, binaryIDs);
                nbTreeAndCountAndBinaryIDList.add(newTuple);
            }
        }
    }

    protected Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }


    protected Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                //System.out.println("Start scoring ...");
                //System.out.println(network2String(network));
                //long start = System.currentTimeMillis();
                Network<Object> speciesNetwork = networkNew2Old(network);

                //double score = findUltrametricOptimalBranchLength(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);
                double score = findNonUltrametricOptimalBranchLength(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);


                if(score > _optimalScores[_optimalNetworks.length-1]){
                    boolean exist = false;
                    for(int i=0; i<_optimalNetworks.length; i++){
                        if(_optimalNetworks[i]==null)break;

                        if(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeClusterDistance(speciesNetwork, _optimalNetworks[i])[2]<0.000001 &&
                                edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeTripartitionDistance(speciesNetwork, _optimalNetworks[i])[2]<0.000001){
                            exist = true;
                            break;
                        }
                    }
                    if(!exist){
                        int index = -1;
                        for(int i=0; i<_optimalScores.length; i++){
                            if(score > _optimalScores[i]){
                                index = i;
                                break;
                            }
                        }
                        for(int i=_optimalScores.length-1; i>index; i--){
                            _optimalNetworks[i] = _optimalNetworks[i-1];
                            _optimalScores[i] = _optimalScores[i-1];
                        }
                        _optimalScores[index] = score;
                        _optimalNetworks[index] = string2Network(network2String(speciesNetwork));
                        //System.out.println(network2String(speciesNetwork) + ": "+score);
                    }
                }
                //System.out.println();
                //System.out.println(network2String(speciesNetwork) + ": "+score);
                //System.out.println();
                //System.out.println("End scoring ..." + (System.currentTimeMillis()-start)/1000.0);
                //System.exit(0);
                return score;
            }
        };
    }


    protected double findUltrametricOptimalBranchLength(final Network<Object> speciesNetwork, final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        boolean continueRounds = true; // keep trying to improve network

        Map<NetNode<Object>, Double> node2height = new Hashtable<NetNode<Object>, Double>();
        for(NetNode<Object> parent: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            if(parent.isLeaf()){
                node2height.put(parent, 0.0);
                continue;
            }
            double height = 0;
            for(NetNode<Object> child: parent.getChildren()){
                height = Math.max(height, node2height.get(child));
            }

            height = height + 1;


            node2height.put(parent, height);
            for(NetNode<Object> child: parent.getChildren()){
                if(child.isNetworkNode()){
                    child.setParentProbability(parent,0.5);
                }

                child.setParentDistance(parent, height - node2height.get(child));
            }
        }


        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList));  // records the GTProb of the network at all times
        final Container<Map<NetNode<Object>, Double>> node2heightContainer = new Container<Map<NetNode<Object>, Double>>(node2height);


        int roundIndex = 0;
        for(; roundIndex <_maxRounds && continueRounds; roundIndex++)
        {
            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();
            for(final NetNode<Object> node : edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork))
            {
                if(node.isLeaf()){
                    continue;
                }

                final Container<Double> minHeight = new Container<Double>(0.0);
                final Container<Double> maxHeight = new Container<Double>(Double.MAX_VALUE);


                for(NetNode<Object> child: node.getChildren()){
                    double childHeight = node2heightContainer.getContents().get(child);
                    minHeight.setContents(Math.max(minHeight.getContents(), childHeight));

                }

                if(!node.isRoot()){
                    for(NetNode<Object> parent: node.getParents()){
                        double parentHeight = node2heightContainer.getContents().get(parent);
                        maxHeight.setContents(Math.min(maxHeight.getContents(), parentHeight));
                    }
                }
                else{
                    maxHeight.setContents(minHeight.getContents() + _maxBranchLength);
                }

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

                        double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, node, true);
                        if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                        {
                            lnGtProbOfSpeciesNetwork.setContents(lnProb);
                            node2heightContainer.getContents().put(node, suggestedHeight);

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

                computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, node, true);

            }

            for(final NetNode<Object> child : speciesNetwork.getNetworkNodes()) // find every hybrid node
            {


                Iterator<NetNode<Object>> hybridParents = child.getParents().iterator();
                final NetNode hybridParent1 = hybridParents.next();
                final NetNode hybridParent2 = hybridParents.next();

                UnivariateFunction functionToOptimize = new UnivariateFunction() {
                    public double value(double suggestedProb) {
                        double incumbentHybridProbParent1 = child.getParentProbability(hybridParent1);
                        child.setParentProbability(hybridParent1, suggestedProb);
                        child.setParentProbability(hybridParent2, 1.0 - suggestedProb);

                        double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, false);

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
                computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, false);



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

        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            if(node.isLeaf()){
                if(species2alleles == null || species2alleles.get(node.getName()).size()<2){
                    node.setParentDistance(node.getParents().iterator().next(), Double.NaN) ;
                }
            }
        }

        //System.out.println(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList) + " vs. " + lnGtProbOfSpeciesNetwork.getContents());
        return lnGtProbOfSpeciesNetwork.getContents();
    }



    protected double findNonUltrametricOptimalBranchLength(final Network<Object> speciesNetwork, final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        boolean continueRounds = true; // keep trying to improve network


        for(NetNode<Object> node: speciesNetwork.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent,1.0);
                if(node.isNetworkNode()){
                    //node.setParentDistance(parent,0.0);
                    node.setParentProbability(parent, 0.5);
                }
            }
        }
        //System.out.println();
        //System.out.println(network2String(speciesNetwork));
        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList));  // records the GTProb of the network at all times
        final Container<Integer> callCount = new Container<Integer>(0);


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
                    if(child.isLeaf()){
                        if(species2alleles==null || species2alleles.get(child.getName()).size()<2){
                            continue;
                        }
                    }

                    assigmentActions.add(new Proc()
                    {
                        public void execute()
                        {

                            UnivariateFunction functionToOptimize = new UnivariateFunction() {
                                public double value(double suggestedBranchLength) {
                                    callCount.setContents(callCount.getContents()+1);
                                    double incumbentBranchLength = child.getParentDistance(parent);

                                    child.setParentDistance(parent, suggestedBranchLength);

                                    double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, parent);
                                    //System.out.println("Changing branch ("+parent.getName()+","+child.getName()+") to " + suggestedBranchLength);
                                    //System.out.println(network2String(speciesNetwork)+": " + lnProb);

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

                            computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, parent);
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
                                callCount.setContents(callCount.getContents()+1);
                                double incumbentHybridProbParent1 = child.getParentProbability(hybridParent1);

                                child.setParentProbability(hybridParent1, suggestedProb);
                                child.setParentProbability(hybridParent2, 1.0 - suggestedProb);

                                double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, false);
                                //System.out.println("Changing node probability to "+ suggestedProb);
                                //System.out.println(network2String(speciesNetwork)+": " + lnProb);
                                //System.out.println(Math.abs(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList) - lnProb));
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
                        computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, false);
                        //System.out.println(network2String(speciesNetwork) + " : " + lnGtProbOfSpeciesNetwork.getContents());
                    }
                });


            }


            // add hybrid probs to hybrid edges


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
        //System.out.println(callCount.getContents());
        //System.out.println(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList) + " vs. " + lnGtProbOfSpeciesNetwork.getContents());
        return lnGtProbOfSpeciesNetwork.getContents();
    }


    protected double computeProbability(Network speciesNetwork, List<Tree> geneTrees, Map<String, List<String>> species2alleles, List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList) {
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        List<Double> probList = gtp.calculateGTDistribution(speciesNetwork, geneTrees, species2alleles);
        double total = 0;
        for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double maxProb = 0;
            for(int id: triple.Item3){
                maxProb = Math.max(maxProb, probList.get(id));
            }
            total += Math.log(maxProb) * triple.Item2;
        }
        return total;
    }




    public double computeProbability(Network speciesNetwork, List<Tree> geneTrees, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList, NetNode child, NetNode parent) {
        Set<NetNode> childNodes = new HashSet<NetNode>();
        childNodes.add(child);
        Set<NetNode> parentNodes = new HashSet<NetNode>();
        parentNodes.add(parent);
        return computeProbability(speciesNetwork, geneTrees, nbTreeAndCountAndBinaryIDList, childNodes, parentNodes);
    }

    public double computeProbability(Network speciesNetwork, List<Tree> geneTrees, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList, NetNode node, boolean changeBranchLength) {
        Set<NetNode> childNodes = new HashSet<NetNode>();
        Set<NetNode> parentNodes = new HashSet<NetNode>();
        if(changeBranchLength){
            //parentNodes = new HashSet<NetNode>();
            parentNodes.add(node);
            for(Object childObject: node.getChildren()){
                childNodes.add((NetNode)childObject);
            }

            childNodes.add(node);
            for(Object parentObject: node.getParents()){
                parentNodes.add((NetNode)parentObject);
            }
        }
        else{
            for(Object parentObject: node.getParents()){
                parentNodes.add((NetNode)parentObject);
            }
            childNodes.add(node);
        }
        return computeProbability(speciesNetwork, geneTrees, nbTreeAndCountAndBinaryIDList, childNodes, parentNodes);
    }

    public double computeProbability(Network speciesNetwork, List<Tree> geneTrees, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList, Set<NetNode> childNodes, Set<NetNode> parentNodes) {
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        List<Double> probList = gtp.calculateGTDistribution(speciesNetwork, geneTrees, childNodes, parentNodes);
        double total = 0;
        for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double maxProb = 0;
            for(int id: triple.Item3){
                maxProb = Math.max(maxProb, probList.get(id));
            }
            total += Math.log(maxProb) * triple.Item2;
        }
        return total;
    }

    protected String network2String(Network speciesNetwork){
        RnNewickPrinter<Double> rnNewickPrinter = new RnNewickPrinter<Double>();
        StringWriter sw = new StringWriter();
        rnNewickPrinter.print(speciesNetwork, sw);
        return sw.toString();
    }

    protected String network2String(final DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork){
        Func1<String,String> _getNetworkNodeLabel  = new Func1<String,String>()
        {
            public String execute(String node) {
                return node;
            }
        };

        Func1<String, Iterable<String>> _getDestinationNodes = new Func1<String, Iterable<String>>() {
            public Iterable<String> execute(String node) {
                return new GetDirectSuccessors<String,PhyloEdge<String>>().execute(speciesNetwork, node);
            }
        };


        Func2<String, String, String> _getNetworkDistanceForPrint = new Func2<String, String, String>() {
            public String execute(String parent, String child) {
                PhyloEdge<String> edge = speciesNetwork.getEdge(parent, child);
                if(edge.getBranchLength()==null){
                    return null;
                }
                return edge.getBranchLength()+"";
            }
        };

        Func2<String, String, String> _getProbabilityForPrint = new Func2<String, String, String>() {
            public String execute(String parent, String child) {
                PhyloEdge<String> edge = speciesNetwork.getEdge(parent, child);
                if(edge.getProbability()==null){
                    return null;
                }
                return edge.getProbability()+"";
            }
        };

        Func2<String, String, String> _getSupportForPrint = new Func2<String, String, String>() {
            public String execute(String parent, String child) {
                PhyloEdge<String> edge = speciesNetwork.getEdge(parent, child);
                if(edge.getSupport()==null){
                    return null;
                }
                return edge.getSupport()+"";
            }
        };

        Func1<String, HybridNodeType> _getHybridTypeForPrint = new Func1<String, HybridNodeType>()
        {
            public HybridNodeType execute(String node)
            {
                int inDegree = new GetInDegree<String,PhyloEdge<String>>().execute(speciesNetwork, node);
                return inDegree == 2 ? HybridNodeType.Hybridization : null;
            }
        };

        Func1<String,String> _getHybridNodeIndexForPrint = new Func1<String, String>() {
            List<String> hybridNodes = new ArrayList<String>();

            public String execute(String node) {
                int inDegree = new GetInDegree<String,PhyloEdge<String>>().execute(speciesNetwork, node);
                if(inDegree == 2){
                    int index = hybridNodes.indexOf(node) + 1;
                    if(index == 0){
                        hybridNodes.add(node);
                        return hybridNodes.size()+"";
                    }
                    else{
                        return index + "";
                    }
                }
                else{
                    return null;
                }
            }
        };

        try{
            StringWriter sw = new StringWriter();
            //   new RichNewickPrinterCompact<String>().print(true, "R", _getNetworkNodeLabel, _getDestinationNodes, _getNetworkDistanceForPrint, _getSupportForPrint, _getProbabilityForPrint, _getHybridNodeIndexForPrint, _getHybridTypeForPrint, sw);
            RichNewickPrinterCompact<String> printer = new RichNewickPrinterCompact<String>();
            printer.setGetBranchLength(_getNetworkDistanceForPrint);
            printer.setGetProbability(_getProbabilityForPrint);
            printer.setGetSupport(_getSupportForPrint);

            printer.print(true, new FindRoot<String>().execute(speciesNetwork), _getNetworkNodeLabel, _getDestinationNodes, _getHybridNodeIndexForPrint, _getHybridTypeForPrint, sw);
            sw.flush();
            sw.close();
            return sw.toString();
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return null;
    }




    protected Network networkNew2Old(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork){
        Network<Object> bniNetwork = null;
        try{

            String networkString = network2String(speciesNetwork);
            bniNetwork = string2Network(networkString);
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return bniNetwork;
    }

    protected Network string2Network(String networkString){
        try{
            RichNewickReaderAST reader = new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
            reader.setHybridSumTolerance(BigDecimal.valueOf(0.00001));
            NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
            RichNewickReadResult<Networks> readResult = reader.read(new ByteArrayInputStream(networkString.getBytes()));
            return transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return null;
    }

}
