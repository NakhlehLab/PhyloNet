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

import edu.rice.cs.bioinfo.library.programming.UnorderedPair;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NonUltrametricNetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.UltrametricNetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
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
public abstract class InferNetworkMLFromGTTBL extends InferNetworkMLFromGT {
    //Map<UnorderedPair, Double> _pairwiseTimeLimit;

    protected void findSingleAlleleSpeciesSet(Network speciesNetwork, Map<String,List<String>> species2alleles, Set<String> singleAlleleSpecies){

    }

    protected NetworkRandomParameterNeighbourGenerator getNetworkRandomParameterNeighbourGenerator(List dataForNetworkInference, Map<String,String> allele2species, Set<String> singleAlleleSpecies){
        //computePairwiseCoalesceTime(dataForNetworkInference, allele2species);
        return new NonUltrametricNetworkRandomParameterNeighbourGenerator();
    }

/*
    private void computePairwiseCoalesceTime(List<Tree> trees, Map<String,String> allele2species){
        _pairwiseTimeLimit = new HashMap<>();

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
                            UnorderedPair<String> sp = new UnorderedPair(taxon1,taxon2);
                            Double minTime = _pairwiseTimeLimit.get(sp);
                            if(minTime == null || minTime > height){
                                _pairwiseTimeLimit.put(sp, height);
                            }
                        }
                    }
                }
                node2leaves.put(node, taxaUnder);
                node2height.put(node, height);
            }
        }

    }


    protected String getStartNetwork(List gts, Map<String,List<String>> species2alleles, Set<String> hybridSpecies, Network<Object> startingNetwork){
        Network network = Networks.readNetwork(super.getStartNetwork(gts, species2alleles, hybridSpecies, startingNetwork));
        Map<NetNode, Double> node2constraints = new HashMap<>();
        computeNodeHeightUpperbound(network, node2constraints);
        initializeNetwork(network, node2constraints);
        return network.toString();
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
                upperBound = node2constraints.get(node.getChildren().iterator().next());
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
                                upperBound = Math.min(upperBound, _pairwiseTimeLimit.get(new UnorderedPair(taxon1, taxon2)));
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


    private void initializeNetwork(Network<Object> speciesNetwork, Map<NetNode, Double> node2constraints){
        Map<NetNode, Double> node2height = new HashMap<>();
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
            if(node2constraints.containsKey(node)){
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

    */

}
