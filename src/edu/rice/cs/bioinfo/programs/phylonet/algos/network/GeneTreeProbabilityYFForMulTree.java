/*
 * Copyright (c) 2012 Rice University.
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


import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;

import java.util.*;


/**
 * Created by Yun Yu
 * Date: 3/13/12
 * Time: 11:31 AM
 *
 * This class is to compute the probability of observing a gene tree given a Multree based on ancestral configurations
 * This method uses only topologies of gene trees.
 * See "Fast Algorithms and Heuristics for Phylogenomics under ILS and Hybridization”, BMC Bioinformatics, 2013.
 * The method is an extension of "Coalescent-based Species Tree Inference from Gene Tree Topologies Under Incomplete Lineage Sorting by Maximum Likelihood", Evolution, 2012. Also see this paper for computation details. Variable names are the same as those appeared in the paper.
 */
public class GeneTreeProbabilityYFForMulTree {
    Set<NetNode> _firstArticulationNodes;
    boolean _printDetails = false;
    int _netNodeNum;

    /**
     * Sets printing option
     */
    public void setPrintDetails(boolean p){
        _printDetails = p;
    }


    /**
     * Computes the probability of observing a gene tree given a Multree
     *
     * @param network 	the given species network
     * @param mulTree	the Multree corresponding to the network
     * @param netLeaf2treeLeaves	the mapping from the leaves of the species network to the leaves of the Multree
     * @param alleleMapping     the mapping from the names of the species from the list of sampled alleles, which is used when multiple alleles are sampled per species
     * @param gt the given gene tree
     *
     * @return	a list of probabilities corresponding to the list of gene trees.
     */
    public double calculateProbability(Network<Object> network, Tree mulTree, Map<NetNode,List<TNode>> netLeaf2treeLeaves, Map<String, List<String>> alleleMapping, Tree gt){
        if(_printDetails){
            System.out.println("Network: " + network);
            System.out.println("MulTree: " + mulTree);
            System.out.println("Gene tree: " + gt);
            System.out.println("Allele mapping: " + alleleMapping);
            System.out.println();
        }


        if(_firstArticulationNodes == null) {
            _firstArticulationNodes = Networks.getLowestArticulationNodes((Network) network);
        }
        _netNodeNum = network.getReticulationCount();
        double gtProb = 0;
        String[] gtTaxa = gt.getLeaves();
        Map<Integer,Integer> child2parent = new HashMap<Integer, Integer>();
        List<STITreeCluster> gtClusters = new ArrayList<STITreeCluster>();
        processGT(gt, gtTaxa, child2parent, gtClusters);
        boolean[][] R = computeR(gtClusters);

        HashSet<String> gtTaxaSet = new HashSet<String>();
        for(String taxon: gtTaxa){
            gtTaxaSet.add(taxon);
        }

        Map<NetNode, Map<NetNode,List<Configuration>>> node2ACMinus = new HashMap<>();
        int reticulationID = 0;
        for(NetNode node: Networks.postTraversal(network)) {
            if(_printDetails){
                System.out.println("\nOn node " + node.getName());
            }
            Map<NetNode,List<Configuration>> CACs = new HashMap<>();
            //compute CACs
            if(node.isLeaf()){
                List<TNode> mulTreeNodes = netLeaf2treeLeaves.get(node);
                computeLeafNodeCACs(node, mulTreeNodes, alleleMapping, gtTaxa, gtTaxaSet, gtClusters, CACs);
            }
            else if(node.isTreeNode()){
                computeInternalTreeNodeCACs(node, node2ACMinus, CACs);
            }
            else{
                computeReticulationNodeCACs(node, node2ACMinus, CACs, reticulationID);
                reticulationID++;
            }

            if(_printDetails){
                System.out.println("AC:");
                for(Map.Entry<NetNode,List<Configuration>> entry: CACs.entrySet()){
                    if(!node.isRoot()) {
                        System.out.println("To parent node " + entry.getKey().getName() + ":");
                    }
                    for(Configuration config: entry.getValue()){
                        System.out.print(config.toString(gtClusters)+"  ");
                    }
                    System.out.println();
                }
            }

            //compute AC Minus
            if(node.isRoot()){
                gtProb = computeFinalProbAtRoot(mulTree.getRoot(), child2parent, R, gtClusters, CACs.values().iterator().next());
                break;
            }
            else {
                Map<NetNode, List<Configuration>> ACMinusMap = new HashMap<>();
                for (Map.Entry<NetNode, List<Configuration>> CAC : CACs.entrySet()) {
                    List<Configuration> ACMinus = new ArrayList<>();
                    double inheritanceProb = node.getParentProbability(CAC.getKey());
                    if(inheritanceProb==NetNode.NO_PROBABILITY){
                        inheritanceProb = 1;
                    }
                    computeACMinus(CAC.getValue(), node.getParentDistance(CAC.getKey()), inheritanceProb, child2parent, R, gtClusters, ACMinus);
                    ACMinusMap.put(CAC.getKey(), ACMinus);
                }
                node2ACMinus.put(node, ACMinusMap);

                if (_printDetails) {
                    System.out.println("AC Minus:");
                    for (Map.Entry<NetNode, List<Configuration>> entry : ACMinusMap.entrySet()) {
                        System.out.println("To parent node " + entry.getKey().getName() + ":");
                        for (Configuration config : entry.getValue()) {
                            System.out.print(config.toString(gtClusters) + "  ");
                        }
                        System.out.println();
                    }

                }
            }

        }

        if(_printDetails){
            System.out.println("\nThe probability of this gene tree under this alleleMapping is:" + gtProb);
            System.out.println();
            System.out.println();
        }

        return gtProb;
    }


    /**
     * Computes the CACs for leaves
     *
     * @param node 	        the leaf node in the original network
     * @param mulTreeNodes  the corresponding leaf nodes in the multree
     * @param alleleMapping	the mapping from the names of species to the list of sampled alleles, which is used when multiple alleles are sampled per species
     * @param gtTaxa	    taxa in the gene tree
     * @param gtTaxaSet	    taxa in the gene tree
     * @param gtClusters	all clusters in the gene tree
     * @param CACs	        the resulting CACs
     */
    private void computeLeafNodeCACs(NetNode node, List<TNode> mulTreeNodes, Map<String, List<String>> alleleMapping, String[] gtTaxa, Set<String> gtTaxaSet, List<STITreeCluster> gtClusters, Map<NetNode,List<Configuration>> CACs){
        Configuration config = new Configuration();
        for(TNode mtNode: mulTreeNodes){
            if(alleleMapping.containsKey(mtNode.getName())) {
                HashSet<Integer> lineages = new HashSet<>();
                for (String allele : alleleMapping.get(mtNode.getName())) {
                    if (gtTaxaSet.contains(allele)) {
                        STITreeCluster cl = new STITreeCluster(gtTaxa);
                        cl.addLeaf(allele);
                        lineages.add(gtClusters.indexOf(cl));
                    }
                }
                if (lineages.size() != 0) {
                    config.addLineage(mtNode, lineages);
                }
            }
        }
        CACs.put((NetNode)node.getParents().iterator().next(), Arrays.asList(config));
    }



    /**
     * Computes the CACs for internal tree nodes by merging two AC- sets
     *
     * @param node 	        the internal tree node in the original network
     * @param node2ACMinus  ACMinus on the outgoing edges of the node that are used to compute CACs
     * @param CACs	        the resulting CACs
     */
    private void computeInternalTreeNodeCACs(NetNode node, Map<NetNode, Map<NetNode,List<Configuration>>> node2ACMinus, Map<NetNode,List<Configuration>> CACs){
        Iterator<NetNode> childNode = node.getChildren().iterator();
        List<Configuration> AC1 = node2ACMinus.get(childNode.next()).get(node);
        List<Configuration> AC2 = node2ACMinus.get(childNode.next()).get(node);
        List<Configuration> CACList = new ArrayList<>();
        Map<Set<Integer>, Configuration> CACMap = new HashMap<>();
        boolean articulate = _firstArticulationNodes.contains(node);
        for(Configuration config1: AC1){
            for(Configuration config2: AC2){
                if(config1.isCompatible(config2)) {
                    Configuration mergedConfig = new Configuration(config1, config2);
                    if(articulate){
                        mergedConfig.clearCoalEvents();
                        Configuration existingConfig = CACMap.get(mergedConfig.getAllLineages());
                        if(existingConfig == null){
                            CACMap.put(mergedConfig.getAllLineages(), mergedConfig);
                        }
                        else{
                            existingConfig.addProbability(mergedConfig.getProbability());
                        }
                    }
                    else{
                        CACList.add(mergedConfig);
                    }
                }
            }
        }

        if(articulate){
            CACList.addAll(CACMap.values());
        }

        if(node.isRoot()){
            CACs.put(null, CACList);
        }
        else {
            CACs.put((NetNode) node.getParents().iterator().next(), CACList);
        }
    }


    /**
     * Computes the CACs for reticulation nodes by splitting an AC- set
     *
     * @param node 	        the reticulation node in the original network
     * @param node2ACMinus  ACMinus on the outgoing edges of the node that are used to compute CACs
     * @param CACs	        the resulting CACs
     * @param reticulationID    id of the reticulation node
     */
    private void computeReticulationNodeCACs(NetNode node, Map<NetNode, Map<NetNode,List<Configuration>>> node2ACMinus, Map<NetNode,List<Configuration>> CACs, int reticulationID){
        Iterator<NetNode> parentNode = node.getParents().iterator();
        NetNode[] parents = new NetNode[2];
        for(int i=0; i<2; i++){
            parents[i] = parentNode.next();
            CACs.put(parents[i], new ArrayList<Configuration>());
        }
        int id = 1;
        for(Configuration config: node2ACMinus.get(node.getChildren().iterator().next()).get(node)){
            config.updateAllParents();
            Configuration[] newConfigs = new Configuration[2];
            for(int i=0; i<2; i++){
                newConfigs[i] = new Configuration();
                newConfigs[i].setNetNodeChoice(config.getNetNodeChoice());
                newConfigs[i].addNetNodeChoice(reticulationID, id);
            }
            newConfigs[0].setProbability(config.getProbability());
            id++;
            int[] lineageCount = new int[2];
            for(Map.Entry<TNode, HashSet<Integer>> lineages: config.getLineages().entrySet()){
                String parentName = lineages.getKey().getName().split("~")[1];
                for(int i=0; i<2; i++){
                    if(parents[i].getName().equals(parentName)){
                        newConfigs[i].addLineage(lineages.getKey(), lineages.getValue());
                        lineageCount[i] += lineages.getValue().size();
                        break;
                    }
                }
            }
            for(int i=0; i<2; i++){
                CACs.get(parents[i]).add(newConfigs[i]);
            }
        }

    }


    /**
     * Computes the AC- given the CACs on a branch (all possible lineages along with their probabilities when leaving a branch given all possible lineages entering a branch)
     *
     * @param CACs 	the given CACs, which are the lineages entering a branch
     * @param distance	the length of the branch
     * @param hybridProb	the inheritance probability of the branch
     * @param child2parent  maps from a child cluster/node to its parent cluster/node (induced from the gene tree topology)
     * @param R     a matrix that keeps the ancestral relationships of all the nodes in a gene tree
     * @param gtClusters    all clusters in the gene tree
     * @param ACMinus	the resulting AC-, which are the lineages leaving a branch
     */
    private void computeACMinus(List<Configuration> CACs, double distance, double hybridProb, Map<Integer, Integer> child2parent, boolean[][] R, List<STITreeCluster> gtClusters, List<Configuration> ACMinus){
        Map<Integer, Map<Configuration, Configuration>> shape2ACminus = new HashMap<Integer, Map<Configuration, Configuration>>();
        for (Configuration origConfig: CACs) {
            Stack<Configuration> configStack = new Stack<Configuration>();
            Configuration configCopy = new Configuration(origConfig);
            configCopy.clearCoalEvents();
            configStack.push(configCopy);
            Set<Set<Integer>> visitedACMinus = new HashSet<>();
            while (!configStack.empty()) {
                Configuration cconfig = configStack.pop();
                double weight = calculateW(cconfig._coalEvents, R);
                cconfig.setProbability(origConfig.getProbability());
                double prob = Math.max(0, computeProbability(origConfig, cconfig, weight, distance, hybridProb, gtClusters));
                cconfig.setProbability(Math.max(0, prob * origConfig.getProbability()));

                int code = cconfig.getLineages().size() * cconfig.getAllLineages().size();
                Map<Configuration, Configuration> existingACMinus = shape2ACminus.get(code);
                if(existingACMinus == null){
                    existingACMinus = new HashMap<>();
                    existingACMinus.put(cconfig, cconfig);
                    shape2ACminus.put(code, existingACMinus);
                }
                else{
                    Configuration existing = existingACMinus.get(cconfig);
                    if(existing == null){
                        existingACMinus.put(cconfig, cconfig);
                    }
                    else{
                        existing.addProbability(cconfig.getProbability());
                    }
                }

                if (cconfig.getLineageCount() == 1) continue;
                for (Map.Entry<TNode, HashSet<Integer>> lineages : cconfig.getLineages().entrySet()){
                    Map<Integer, List<Integer>> parent2children = new HashMap<Integer, List<Integer>>();
                    for(int i: lineages.getValue()){
                        int pid = child2parent.get(i);
                        List<Integer> children = parent2children.get(pid);
                        if (children == null) {
                            children = new ArrayList<Integer>();
                            parent2children.put(pid, children);
                        }
                        children.add(i);
                    }

                    for (Map.Entry<Integer, List<Integer>> entry : parent2children.entrySet()) {
                        Integer parent = entry.getKey();
                        List<Integer> childrenList = entry.getValue();
                        if (childrenList.size() <= 1) {
                            continue;
                        }
                        Configuration newConfig = new Configuration(cconfig);
                        newConfig.mergeCluster(lineages.getKey(), childrenList.get(0), childrenList.get(1), parent);
                        newConfig.addCoalEvents(cconfig._coalEvents);
                        newConfig.addCoalEvent(parent);
                        Set<Integer> allLineages = newConfig.getAllLineages();
                        if(!visitedACMinus.contains(allLineages)){
                            visitedACMinus.add(allLineages);
                            configStack.push(newConfig);
                        }
                    }
                }
            }
        }

        for (Map<Configuration, Configuration> cc : shape2ACminus.values()) {
            ACMinus.addAll(cc.values());
        }

    }


    /**
     * Computes final probability at the root of the network/multree
     *
     * @param root          the root of the network/multree
     * @param child2parent  maps from a child cluster/node to its parent cluster/node (induced from the gene tree topology)
     * @param R             a matrix that keeps the ancestral relationships of all the nodes in a gene tree
     * @param gtClusters    all clusters in the gene tree
     * @param CACs 	        the ancestral configurations entering a branch
     */
    private double computeFinalProbAtRoot(TNode root, Map<Integer, Integer> child2parent, boolean[][] R, List<STITreeCluster> gtClusters, List<Configuration> CACs){
        double gtProb = 0;
        Configuration rootConfig = new Configuration();
        HashSet<Integer> rootLineageSet = new HashSet<>();
        rootLineageSet.add(gtClusters.size()-1);
        rootConfig.addLineage(root, rootLineageSet);
        for(Configuration preConfig: CACs){
            if(preConfig.getLineageCount()==1){
                gtProb += preConfig._totalProb;
            }else{
                Set<Integer> events = new HashSet<Integer>();
                events.add(gtClusters.size()-1);
                for(int i : preConfig.getAllLineages()) {
                    int p = child2parent.get(i);
                    while (!events.contains(p)) {
                        events.add(p);
                        p = child2parent.get(p);
                    }
                }
                double weight = calculateW(events, R);
                double prob = Math.max(0, computeProbability(preConfig, rootConfig, weight, -1, 1, gtClusters));
                gtProb += Math.max(0, prob * preConfig.getProbability());
            }
        }
        return gtProb;
    }



    /**
     * This function is to pre-process a gene tree, especially getting the ancestral relationships
     *
     * @param gt 	the gene tree to be processed
     * @param child2parent  the resulting mapping from a child cluster/node to its parent cluster/node (induced from the gene tree topology)
     * @param gtClusters    the resulting list of clusters in the gene tree
     */
    private void processGT(Tree gt, String[] gtTaxa, Map<Integer,Integer> child2parent,List<STITreeCluster> gtClusters){
        Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
        int index = 0;
        for (TNode node : gt.postTraverse()) {
            ((STINode<Integer>)node).setData(index);
            BitSet bs = new BitSet();
            if (node.isLeaf()) {
                for (int i = 0; i < gtTaxa.length; i++) {
                    if (node.getName().equals(gtTaxa[i])) {
                        bs.set(i);
                        break;
                    }
                }
            }
            else {
                for (TNode child : node.getChildren()) {
                    bs.or(map.get(child));
                    child2parent.put(((STINode<Integer>)child).getData(), index);
                }
            }
            map.put(node, bs);
            STITreeCluster cl = new STITreeCluster(gtTaxa);
            cl.setCluster(bs);
            gtClusters.add(cl);
            index++;
        }
    }




    /**
     * The function is to calculate the _R matrix for the given tree.
     * Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
     */
    private boolean[][] computeR(List<STITreeCluster> gtClusters){
        boolean[][] R = new boolean[gtClusters.size()][gtClusters.size()];
        for(int i=0; i<gtClusters.size(); i++){
            STITreeCluster cl1 = gtClusters.get(i);
            for(int j=i+1; j<gtClusters.size(); j++){
                STITreeCluster cl2 = gtClusters.get(j);
                if(cl1.containsCluster(cl2)){
                    R[i][j] = true;
                }
                else if(cl2.containsCluster(cl1)){
                    R[j][i] = true;
                }
            }
        }
        return R;
    }



    /**
     * Computes the probability of a given ancestral configuration coalescing into another given one
     *
     * @param preConfig 	the original ancestral configuration
     * @param coalescedConfig     the ancestral configuration that the original one is coalescing into
     * @param w     weight
     * @param distance 	the branch length
     * @param portion	the inheritance probability
     * @param gtClusters    all clusters in the gene tree
     *
     * @return the probability
     */
    private double computeProbability(Configuration preConfig, Configuration coalescedConfig, double w, double distance, double portion, List<STITreeCluster> gtClusters){
        double prob = 1.0;
        int u = preConfig.getLineageCount();
        int v = coalescedConfig.getLineageCount();
        if(distance == 0 && u != v){
            if(_printDetails){
                System.out.println(preConfig.toString(gtClusters)+"->"+coalescedConfig.toString(gtClusters)+ ": g"+u+v+"(0)=0");
            }
            return 0;
        }
        if(u == v && (u == 1 || u == 0)){
            if(_printDetails){
                if(portion==1 || u==0){
                    System.out.println(preConfig.toString(gtClusters)+"->"+coalescedConfig.toString(gtClusters)+ ": g"+u+v+"("+distance+")=1");
                }
                else{
                    System.out.println(preConfig.toString(gtClusters)+"->"+coalescedConfig.toString(gtClusters)+ ": g"+u+v+"("+distance+")*"+portion+"="+portion);
                }
            }
            return prob;
        }
        double gij = gij(distance, u, v);
        if(gij<0){
            return -1;
        }
        double d = calculateD(u,u-v);

        if(_printDetails){
            if(portion!=1){
                System.out.println(preConfig.toString(gtClusters)+"->"+coalescedConfig.toString(gtClusters)+ ": g"+u+v+"("+distance+")*"+w+"/"+d+"*"+portion+"^"+u+"=" + gij*w/d*Math.pow(portion, u));
            }
            else{
                System.out.println(preConfig.toString(gtClusters)+"->"+coalescedConfig.toString(gtClusters)+ ": g"+u+v+"("+distance+")*"+w+"/"+d+"=" + gij*w/d);
            }
        }
        prob *= gij*w/d;
        return prob;
    }



    /**
     * The function is to calculate the g_{ij} function.
     * @param	length	the branch length
     * @param 	i	the number of lineages in
     * @param	j	the number of lineages out
     * @return	the resulting probability
     */
    private double gij(double length, int i, int j){
        if(length == TNode.NO_DISTANCE || length == -1){
            if(j == 1){
                return 1;
            }else{
                return 0;
            }
        }
        if(length==0){
            if(i==j)
                return 1;
            else
                return 0;
        }
        if(i==0){
            return 1;
        }

        double result = 0;
        for(int k=j; k<=i; k++){

            double temp = Math.exp(0.5 * k * (1.0 - k) * length) * (2.0 * k - 1.0) * Math.pow(-1, k - j)*(fact(j, j + k - 2))*(fact(i - k + 1, i));
            double denom = fact(1, j)*(fact(1, k - j))*(fact(i, i + k - 1));
            result += temp/denom;
        }
        return result;
    }


    /**
     * The function is to calculate the number of all possible ordering coalescent events
     * @param	u	the number of lineages entering the branch
     * @param	c	the number of coalescent events happening on the branch
     */
    private double calculateD(int u, int c){
        double d = 1;
        if(c!=0){
            for(int i=1; i<=c; i++){
                d = d * chooseD(u - i + 1, 2);
            }
        }
        return d;
    }


    /**
     * The function is to calculate "N choose K"
     */
    private double chooseD(int N, int K) {
        double ret = 1.0;
        for (int k = 0; k < K; k++) {
            ret = ret*((N-k+0.0)/(k+1));
        }

        return ret;
    }


    /**
     * The function is to calculate the number of ways (weight) that coalescent events on a branch can occur consistently with the gene tree
     * @param coalEvents	the coalescent events happened on the branch
     * @param R     the ancestral relationships of nodes in a gene tree
     */
    private double calculateW(Set<Integer> coalEvents, boolean[][] R){
        double w = 1.0;
        w = w * fact(1, coalEvents.size());
        for (int i: coalEvents) {
            int a = 0;
            for (int j: coalEvents) {
                if(i!=j && R[i][j])
                    a++;
            }
            w = w * (1.0/(1 + a));
        }
        return w;
    }


    /**
     * The function is to calculate factorial
     * @param	start	the first number
     * @param 	end		the last number
     *
     * @return	the resulting factorial
     */
    private double fact(int start, int end){
        double result = 1;
        for(int i=start; i<=end; i++){
            result = result*i;
        }

        return result;
    }


    /**
     * This class is to represent the concept of ancestral configuration
     */
    private class Configuration{
        private Map<TNode,HashSet<Integer>> _lineages;
        private double _totalProb;
        private int[] _netNodeIndex;
        private Set<Integer> _coalEvents;
        private HashSet<Integer> _allLineages;


        public Configuration(){
            _lineages = new HashMap<>();
            _totalProb = 1;
            _netNodeIndex = new int[_netNodeNum];
            Arrays.fill(_netNodeIndex, 0);
            _allLineages = new HashSet<>();
        }

        public Configuration(Configuration config){
            _lineages = new HashMap<>();
            for(Map.Entry<TNode,HashSet<Integer>> entry: config.getLineages().entrySet()){
                _lineages.put(entry.getKey(), (HashSet)entry.getValue().clone());
            }
            _totalProb = config.getProbability();
            _netNodeIndex = config.getNetNodeChoice().clone();
            _allLineages = (HashSet)config.getAllLineages().clone();
        }

        public Configuration(Configuration config1, Configuration config2){
            _lineages = new HashMap<>();
            for(Map.Entry<TNode,HashSet<Integer>> entry: config1.getLineages().entrySet()){
                _lineages.put(entry.getKey().getParent(), (HashSet)entry.getValue().clone());
            }
            for(Map.Entry<TNode,HashSet<Integer>> entry: config2.getLineages().entrySet()){
                HashSet<Integer> existingLineages = _lineages.get(entry.getKey().getParent());
                if(existingLineages == null){
                    existingLineages = new HashSet<>();
                    _lineages.put(entry.getKey().getParent(), existingLineages);
                }
                existingLineages.addAll(entry.getValue());
            }
            _totalProb = Math.max(0, config1._totalProb * config2._totalProb);
            _netNodeIndex = new int[_netNodeNum];
            for(int i=0; i< _netNodeNum; i++){
                if(config1._netNodeIndex[i] == config2._netNodeIndex[i]){
                    _netNodeIndex[i] = config1._netNodeIndex[i];
                }
                else{
                    _netNodeIndex[i] = Math.max(config1._netNodeIndex[i], config2._netNodeIndex[i]);
                }
            }
            _allLineages = (HashSet)config1.getAllLineages().clone();
            _allLineages.addAll(config2.getAllLineages());
        }


        public boolean isCompatible(Configuration config){
            boolean compatible = true;
            for(int i=0; i< _netNodeNum; i++){
                if(_netNodeIndex[i] != config._netNodeIndex[i] && _netNodeIndex[i]!=0 && config._netNodeIndex[i]!=0){
                    compatible = false;
                    break;
                }
            }
            return compatible;
        }

        public void updateAllParents(){
            Map<TNode,HashSet<Integer>> temp = _lineages;
            _lineages = new HashMap<>();
            for(Map.Entry<TNode,HashSet<Integer>> entry: temp.entrySet()){
                _lineages.put(entry.getKey().getParent(), (HashSet)entry.getValue().clone());
            }
        }

        public Map<TNode,HashSet<Integer>> getLineages(){
            return _lineages;
        }

        public HashSet<Integer> getAllLineages(){
            return _allLineages;
        }

        public void addLineage(TNode node, HashSet<Integer> lineages){
            _lineages.put(node, lineages);
            _allLineages.addAll(lineages);
        }


        public int getLineageCount(){
           return _allLineages.size();
        }


        public void mergeCluster(TNode node, int child1, int child2, int parent){
            Set<Integer> lineageSet = _lineages.get(node);
            lineageSet.remove(child1);
            lineageSet.remove(child2);
            lineageSet.add(parent);
            _allLineages.remove(child1);
            _allLineages.remove(child2);
            _allLineages.add(parent);
        }

        public void addProbability(double prob){
            _totalProb += prob;
        }

        public void setProbability(double prob){
            _totalProb = prob;
        }

        public void multiplyProbability(double value){
            _totalProb *= value;
        }

        public double getProbability(){
           return _totalProb;
        }


        public int[] getNetNodeChoice(){
            return _netNodeIndex;
        }

        public void addNetNodeChoice(int net, int index){
            _netNodeIndex[net] = index;
        }

        public void setNetNodeChoice(int[] choice){
            _netNodeIndex = choice.clone();
        }

        public String toString(List<STITreeCluster> gtClusters){
            String exp = "{";
            for(Map.Entry<TNode,HashSet<Integer>> entry: _lineages.entrySet()){
                exp += entry.getKey().getName() + ":[";
                for(int id: entry.getValue()) {
                    exp += gtClusters.get(id) + ",";
                }
                exp = exp.substring(0, exp.lastIndexOf(","));
                exp += "], ";
            }
            if(_lineages.size() != 0) {
                exp = exp.substring(0, exp.lastIndexOf(","));
            }
            exp = exp + "}/" + Arrays.toString(_netNodeIndex);
            exp = exp + ":" + _totalProb;
            return exp;
        }

        public void clearCoalEvents(){
            if(_coalEvents==null){
                _coalEvents = new HashSet<Integer>();
            }
            else{
                _coalEvents.clear();
            }
        }

        public void addCoalEvent(int index){
            if(_coalEvents==null){
                _coalEvents = new HashSet<Integer>();
            }
            _coalEvents.add(index);
        }

        public void addCoalEvents(Set<Integer> events){
            if(_coalEvents==null){
                _coalEvents = new HashSet<Integer>();
            }
            _coalEvents.addAll(events);
        }

        public boolean equals(Object o) {
            if(!(o instanceof Configuration)){
                return false;
            }

            Configuration config = (Configuration) o;
            return config.getAllLineages().equals(_allLineages) && Arrays.equals(config.getNetNodeChoice(), _netNodeIndex);
        }

        public int hashCode(){
            return _allLineages.hashCode() + Arrays.hashCode(_netNodeIndex);
        }

    }



}
