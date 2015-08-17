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


import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.SimpleHillClimbing;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomTopologyNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NonUltrametricNetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferNetworkMP{
    protected int _maxFailure = 100;
    protected long _maxExaminations = -1;
    protected int _moveDiameter = -1;
    protected int _reticulationDiameter = -1;
    private Set<String> _fixedHybrid = new HashSet<String>();
    protected Network<Object> _startNetwork = null;
    private int _numProcessors = 1;
    protected double[] _operationWeight = {0.1,0.1,0.15,0.55,0.15,0.15};
    private int _numRuns = 1;
    private Long _seed = null;


    public void setSearchParameter(long maxExaminations, int maxFailure, int moveDiameter, int reticulationDiameter, Network startNetwork, Set<String> fixedHybrid, int numProcessors, double[] operationWeight, int numRuns, Long seed){
        _maxExaminations = maxExaminations;
        _moveDiameter = moveDiameter;
        _reticulationDiameter = reticulationDiameter;
        _startNetwork = startNetwork;
        _maxFailure = maxFailure;
        _fixedHybrid = fixedHybrid;
        _numProcessors = numProcessors;
        _operationWeight = operationWeight;
        _numRuns = numRuns;
        _seed = seed;
    }



    private void checkNetworkWithHybrids(Network<Object> startNetwork){
        if(startNetwork == null || _fixedHybrid.size() == 0){
            return;
        }

        if(startNetwork != null){
            for(NetNode<Object> node: startNetwork.getNetworkNodes()){
                NetNode hybridSpecies = node.getChildren().iterator().next();
                if(!(hybridSpecies.isLeaf() && _fixedHybrid.contains(hybridSpecies.getName()))){
                    throw new IllegalArgumentException("The starting network contains hybrid that is not in the specified hybrid set.");
                }
            }
        }
    }


    public List<Tuple<Network,Double>> inferNetwork(List<MutableTuple<Tree,Double>> gts, Map<String,List<String>> species2alleles, int maxReticulations, int numSol){
        LinkedList<Tuple<Network, Double>> resultList = new LinkedList<>();
        String startingNetwork = getStartNetwork(gts, species2alleles,_fixedHybrid, _startNetwork);

        NetworkRandomNeighbourGenerator allNeighboursStrategy = new NetworkRandomNeighbourGenerator(new NetworkRandomTopologyNeighbourGenerator(_operationWeight, maxReticulations, _moveDiameter, _reticulationDiameter, _fixedHybrid, _seed), 1, new NonUltrametricNetworkRandomParameterNeighbourGenerator(), 0, _seed);
        Comparator<Double> comparator = getDoubleScoreComparator();
        SimpleHillClimbing searcher = new SimpleHillClimbing(comparator, allNeighboursStrategy);
        Func1<Network, Double> scorer = getScoreFunction(gts, species2alleles);
        Network speciesNetwork = Networks.readNetwork(startingNetwork);
        searcher.search(speciesNetwork, scorer, numSol, _numRuns, _maxExaminations, _maxFailure, true, resultList); // search starts here
        //To set inheritance probability
        if(_numProcessors != 1){
            for(Tuple<Network, Double> tuple: resultList){
                MDCOnNetworkYF mdc = new MDCOnNetworkYF();
                mdc.countExtraCoal(tuple.Item1, gts, species2alleles, new int[gts.size()]);
            }
        }
        return resultList;
    }



    protected void createHybrid(Network<Object> network, String hybrid){
        List<Tuple<NetNode,NetNode>> edgeList = new ArrayList<Tuple<NetNode,NetNode>>();
        Tuple<NetNode,NetNode> destinationEdge = null;
        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(network)){
            for(NetNode child: node.getChildren()){
                if(child.isLeaf() && child.getName().equals(hybrid)){
                    if(node.isNetworkNode()){
                        return;
                    }
                    destinationEdge = new Tuple<NetNode, NetNode>(node, child);
                }
                else{
                    edgeList.add(new Tuple<NetNode, NetNode>(node, child));
                }
            }

        }

        int numEdges = edgeList.size();
        Tuple<NetNode,NetNode> sourceEdge = edgeList.get((int)(Math.random() * numEdges));
        NetNode insertedSourceNode = new BniNetNode();
        insertedSourceNode.adoptChild(sourceEdge.Item2, NetNode.NO_DISTANCE);
        sourceEdge.Item1.removeChild(sourceEdge.Item2);
        sourceEdge.Item1.adoptChild(insertedSourceNode, NetNode.NO_DISTANCE);
        NetNode insertedDestinationNode = new BniNetNode();
        insertedDestinationNode.adoptChild(destinationEdge.Item2, NetNode.NO_DISTANCE);
        destinationEdge.Item1.removeChild(destinationEdge.Item2);
        destinationEdge.Item1.adoptChild(insertedDestinationNode, NetNode.NO_DISTANCE);
        insertedSourceNode.adoptChild(insertedDestinationNode, NetNode.NO_DISTANCE);
    }


    private String getStartNetwork(List<MutableTuple<Tree,Double>> gts, Map<String,List<String>> species2alleles, Set<String> hybridSpecies, Network<Object> startingNetwork){
        checkNetworkWithHybrids(startingNetwork);

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
                sol = mdc.inferSpeciesTree(gts, false, 1, false, true, -1).get(0);
            }
            else{
                sol = mdc.inferSpeciesTree(gts, allele2species, false, 1, false, true, -1).get(0);
            }
            Tree startingTree= Trees.generateRandomBinaryResolution(sol._st);
            startingNetwork = Networks.readNetwork(startingTree.toString());
        }

        for(String hybrid: hybridSpecies){
            createHybrid(startingNetwork, hybrid);
        }

        Networks.removeAllParameters(startingNetwork);
        Networks.autoLabelNodes(startingNetwork);

        String newNetwork = startingNetwork.toString();

        return newNetwork;

    }


    private Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o2, o1);
            }
        };
    }



    private int[] computeXLParallel(Network network, List<MutableTuple<Tree,Double>> gts, Map<String, List<String>> species2alleles, int[] xls){
        Thread[] myThreads = new Thread[_numProcessors];

        MDCOnNetworkYF mdc = new MDCOnNetworkYF();
        mdc.setParallel(true);
        mdc.processNetwork(network);

        for(int i=0; i<_numProcessors; i++){
            myThreads[i] = new MyThread(network, gts, species2alleles, xls, mdc);
            myThreads[i].start();
        }

        for(int i=0; i<_numProcessors; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }

        return xls;

    }



    private Func1<Network, Double> getScoreFunction(final List<MutableTuple<Tree,Double>> gts, final Map<String, List<String>> species2alleles){
        return new Func1<Network, Double>() {
            public Double execute(Network network) {
                int[] scores = new int[gts.size()];
                if(_numProcessors == 1){
                    MDCOnNetworkYF scorer = new MDCOnNetworkYF();
                    scorer.countExtraCoal(network, gts, species2alleles, scores);
                }
                else{
                    computeXLParallel(network, gts, species2alleles, scores);
                }

                double total = 0;
                Iterator<MutableTuple<Tree,Double>> weightIt = gts.iterator();
                for(int score: scores){
                    total += score * weightIt.next().Item2;
                }

                return total;

            }
        };
    }



    private class MyThread extends Thread{
        Network _speciesNetwork;
        List<MutableTuple<Tree,Double>> _geneTrees;
        Map<String, List<String>> _species2alleles;
        int[] _xls;
        MDCOnNetworkYF _mdc;


        public MyThread(Network speciesNetwork, List<MutableTuple<Tree,Double>> geneTrees, Map<String, List<String>> species2alleles, int[] xls, MDCOnNetworkYF mdc){
            _speciesNetwork = speciesNetwork;
            _geneTrees = geneTrees;
            _species2alleles = species2alleles;
            _xls = xls;
            _mdc = mdc;
        }


        public void run() {
            _mdc.countExtraCoal(_speciesNetwork, _geneTrees, _species2alleles, _xls);
        }
    }


}
