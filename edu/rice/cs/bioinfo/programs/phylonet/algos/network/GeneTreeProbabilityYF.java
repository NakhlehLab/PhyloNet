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


import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
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
 * Created by yy9
 * Date: 5/10/12
 * Time: 2:31 PM
 *
 * This class is to compute the probability of observing a collection of gene trees given a species network based on ancestral configuration
 * This method uses only topologies of gene trees. It does not cache ancestral configurations.
 * See "Fast Algorithms and Heuristics for Phylogenomics under ILS and Hybridization”, BMC Bioinformatics, 2013.
 * The method is an extension of "Coalescent-based Species Tree Inference from Gene Tree Topologies Under Incomplete Lineage Sorting by Maximum Likelihood", Evolution, 2012. Also see this paper for computation details. Variable names are the same as those appeared in the paper.
 */

public class GeneTreeProbabilityYF {
    Set<NetNode> _totalCoverNodes;
    boolean _printDetails = false;
    int _netNodeNum;
    Map<NetNode, Integer> _node2ID;
    boolean _parallel = false;
    int _currentTreeID = 0;
    int _totalTree;
    boolean _preProcessed = false;


    /**
     * Sets parallel computing
     */
    public void setParallel(boolean parallel){
        _parallel = parallel;
    }


    /**
     * Sets printing option
     */
    public void setPrintDetails(boolean p){
        _printDetails = p;
    }


    /**
     * Gets the next triplet to compute, which is used for parallel computing
     */
    public synchronized int getNextTreeID(){
        return _currentTreeID++;
    }


    /**
     * Do all the pre-processing
     */
    public void preProcess(Network network, List<Tree> gts, boolean fromScratch){
        _totalTree = gts.size();
        processNetwork(network, fromScratch);
        _preProcessed = true;
    }


    /**
     * Computes the probability of observing a collection of gene trees given a species network
     *
     * @param network 	the given species network
     * @param gts	the given collection of gene trees
     * @param species2alleles	the mapping from a species to the list of alleles sampled, which is used when multiple alleles are sampled per species
     * @param resultProbs	the resulting probabilities
     */
    public void calculateGTDistribution(Network<Object> network, List<Tree> gts, Map<String, List<String>> species2alleles, double resultProbs[]){
        if(!_preProcessed){
            preProcess(network, gts, true);
        }

        int treeID = 0;
        if(_parallel){
            treeID = getNextTreeID();
        }

        while(treeID < _totalTree){
            Tree gt = gts.get(treeID);
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

            int netNodeIndex = 0;
            Map<NetNode, Map<NetNode,List<Configuration>>> node2ACMinus = new HashMap<>();


            for(NetNode node: Networks.postTraversal(network)){
                if(_printDetails){
                    System.out.println();
                    System.out.println("On node #" + _node2ID.get(node) + " " + node.getName());
                }
                List<Map<Set<Integer>,List<Configuration>>> CACs = new ArrayList<Map<Set<Integer>,List<Configuration>>>();

                //set AC for a node
                if(node.isLeaf()){
                    Map<Set<Integer>,List<Configuration>> sizeOneConfigs = new HashMap<Set<Integer>, List<Configuration>>();
                    Configuration config = new Configuration();
                    if(species2alleles == null){
                        if(gtTaxaSet.contains(node.getName())){
                            STITreeCluster cl = new STITreeCluster(gtTaxa);
                            cl.addLeaf(node.getName());
                            config.addLineage(gtClusters.indexOf(cl));
                        }
                    }
                    else{
                        for(String allele: species2alleles.get(node.getName())){
                            if(gtTaxaSet.contains(allele)){
                                STITreeCluster cl = new STITreeCluster(gtTaxa);
                                cl.addLeaf(allele);
                                config.addLineage(gtClusters.indexOf(cl));
                            }
                        }
                    }
                    config.setTotalProbability(1);
                    List<Configuration> tempList = new ArrayList<Configuration>();
                    tempList.add(config);
                    sizeOneConfigs.put(config._lineages, tempList);
                    CACs.add(sizeOneConfigs);
                }
                else{
                    if(node.isNetworkNode()){
                        Map<NetNode, List<Configuration>> childACMinus = node2ACMinus.get(node.getChildren().iterator().next());
                        Map<Set<Integer>,List<Configuration>> temp = new HashMap<Set<Integer>, List<Configuration>>();
                        temp.put(null,childACMinus.get(node));
                        CACs.add(temp);
                    }
                    else{
                        Iterator<NetNode> childNode = node.getChildren().iterator();
                        List<Configuration> AC1 = node2ACMinus.get(childNode.next()).get(node);
                        List<Configuration> AC2 = node2ACMinus.get(childNode.next()).get(node);

                        List<Integer> configSizeList = new ArrayList<Integer>();
                        boolean total = _totalCoverNodes.contains(node);
                        List<Configuration> mergedConfigList = new ArrayList<Configuration>();
                        for(Configuration config1: AC1){
                            for(Configuration config2: AC2){
                                if(config1.isCompatible(config2)){
                                    Configuration mergedConfig = new Configuration(config1, config2);
                                    int numLin = mergedConfig._lineages.size();
                                    List<Configuration> sameLineageConfigs;
                                    int sizeIndex = configSizeList.indexOf(numLin);
                                    if(sizeIndex == -1){
                                        int pos = 0;
                                        for(Integer size: configSizeList){
                                            if(size >= numLin){
                                                break;
                                            }
                                            pos++;
                                        }
                                        sameLineageConfigs = new ArrayList<Configuration>();
                                        Map<Set<Integer>, List<Configuration>> sameSizelineages2configs = new HashMap<Set<Integer>, List<Configuration>>();
                                        sameSizelineages2configs.put(mergedConfig._lineages,sameLineageConfigs);
                                        CACs.add(pos, sameSizelineages2configs);
                                        configSizeList.add(pos, numLin);
                                    }
                                    else{
                                        Map<Set<Integer>, List<Configuration>> sameSizelineages2configs = CACs.get(sizeIndex);
                                        sameLineageConfigs = sameSizelineages2configs.get(mergedConfig._lineages);
                                        if(sameLineageConfigs==null){
                                            sameLineageConfigs = new ArrayList<Configuration>();
                                            sameSizelineages2configs.put(mergedConfig._lineages,sameLineageConfigs);
                                        }
                                    }

                                    if(total){
                                        mergedConfig.clearNetNodeChoice();
                                        if(sameLineageConfigs.size()==0){
                                            sameLineageConfigs.add(mergedConfig);
                                            mergedConfigList.add(mergedConfig);
                                        }
                                        else{
                                            sameLineageConfigs.get(0).addTotalProbability(mergedConfig._totalProb);
                                        }
                                    }
                                    else{
                                        sameLineageConfigs.add(mergedConfig);
                                        mergedConfigList.add(mergedConfig);
                                    }

                                }

                            }
                        }

                    }
                }

                if(_printDetails){
                    System.out.print("AC: {");
                    for(Map<Set<Integer>,List<Configuration>> lineages2configs: CACs){
                        for(List<Configuration> configList: lineages2configs.values()){
                            for(Configuration config: configList){
                                System.out.print(config.toString(gtClusters)+"  ");
                            }
                        }
                    }
                    System.out.println("}");
                }


                //set AC- for a node
                if(node.isRoot()){
                    Configuration rootConfig = new Configuration();
                    rootConfig.addLineage(gtClusters.size()-1);
                    for(Map<Set<Integer>,List<Configuration>> lineages2configs: CACs){
                        for(List<Configuration> configList: lineages2configs.values()){
                            for(Configuration preConfig: configList){
                                if(preConfig.getLineageCount()==1){
                                    gtProb += preConfig._totalProb;
                                }else{
                                    Set<Integer> events = new HashSet<Integer>();
                                    events.add(gtClusters.size()-1);
                                    for (int i: preConfig._lineages) {
                                        int p = child2parent.get(i);

                                        while(!events.contains(p)){
                                            events.add(p);
                                            p = child2parent.get(p);
                                        }
                                    }
                                    double weight = calculateW(events, R);
                                    double prob = Math.max(0, computeProbability(preConfig, rootConfig, weight, -1, 1, gtClusters));
                                    gtProb += Math.max(0, prob*preConfig._totalProb);
                                }

                            }
                        }

                    }
                    List<Configuration> temp = new ArrayList<Configuration>();
                    temp.add(rootConfig);

                }
                else if(node.isTreeNode()){
                    NetNode parentNode = (NetNode)node.getParents().iterator().next();
                    double distance = node.getParentDistance(parentNode);
                    List<Configuration> ACminus = new ArrayList<Configuration>();
                    computeACMinus(CACs, distance, 1, child2parent, R, gtClusters, ACminus);
                    Map<NetNode, List<Configuration>> ACminusMap = new HashMap<>();
                    ACminusMap.put(parentNode, ACminus);
                    node2ACMinus.put(node, ACminusMap);
                    if(_printDetails){
                        System.out.print("ACminus: {");
                        for(Configuration config: ACminus){
                            System.out.print(config.toString(gtClusters)+"  ");
                        }
                        System.out.println("}");
                    }
                }
                else {
                    List<Integer> configSizeList = new ArrayList<Integer>();
                    List<Map<Set<Integer>,List<Configuration>>> newCACs1 = new ArrayList<Map<Set<Integer>,List<Configuration>>>();
                    List<Map<Set<Integer>,List<Configuration>>> newCACs2 = new ArrayList<Map<Set<Integer>,List<Configuration>>>();
                    Map<Configuration, List<Configuration>> config2splitedConfigs = new HashMap<Configuration, List<Configuration>>();
                    int configIndex = 1;
                    for(Map<Set<Integer>,List<Configuration>> lineages2configs: CACs){
                        for(List<Configuration> configList: lineages2configs.values()){
                            for(Configuration config: configList){
                                List<Configuration> splitedConfigs = new ArrayList<Configuration>();
                                config2splitedConfigs.put(config, splitedConfigs);
                                int[] lineageArray = new int[config.getLineageCount()];
                                int index = 0;
                                for(int lineage: config._lineages){
                                    lineageArray[index++] = lineage;
                                }
                                boolean fEven = (config.getLineageCount()%2)==0;
                                int upper = Math.max(0, config.getLineageCount()/2);
                                Set<Set<Integer>> addedConfigs = null;
                                if(fEven){
                                    addedConfigs = new HashSet<Set<Integer>>();
                                }

                                for(int i=0; i<=upper; i++){
                                    for(boolean[] selectedLineages: getSelected(config.getLineageCount(),i)){
                                        for(int k=0; k<2; k++){
                                            Configuration newConfig = new Configuration();
                                            index = 0;
                                            for(int lin: lineageArray) {
                                                if(selectedLineages[index] && k==0){
                                                    newConfig.addLineage(lin);
                                                }
                                                else if(!selectedLineages[index] && k==1){
                                                    newConfig.addLineage(lin);
                                                }
                                                index ++;
                                            }

                                            if(fEven && i==upper){
                                                if(addedConfigs.contains(newConfig._lineages)){
                                                    break;
                                                }else{
                                                    addedConfigs.add(newConfig._lineages);
                                                }
                                            }

                                            newConfig.setNetNodeChoice(config._netNodeIndex);
                                            newConfig.addNetNodeChoice(netNodeIndex, configIndex);
                                            newConfig.setTotalProbability(Math.sqrt(config._totalProb));

                                            int numLin = k==0?i : config.getLineageCount()-i;
                                            List<Configuration> sameLineageConfigs1;
                                            List<Configuration> sameLineageConfigs2;
                                            int sizeIndex = configSizeList.indexOf(numLin);
                                            if(sizeIndex == -1){
                                                int pos = 0;
                                                for(Integer size: configSizeList){
                                                    if(size >= numLin){
                                                        break;
                                                    }
                                                    pos++;
                                                }
                                                sameLineageConfigs1 = new ArrayList<Configuration>();
                                                Map<Set<Integer>, List<Configuration>> sameSizelineages2configs1 = new HashMap<Set<Integer>, List<Configuration>>();
                                                sameSizelineages2configs1.put(newConfig._lineages,sameLineageConfigs1);
                                                newCACs1.add(pos, sameSizelineages2configs1);

                                                sameLineageConfigs2 = new ArrayList<Configuration>();
                                                Map<Set<Integer>, List<Configuration>> sameSizelineages2configs2 = new HashMap<Set<Integer>, List<Configuration>>();
                                                sameSizelineages2configs2.put(newConfig._lineages,sameLineageConfigs2);
                                                newCACs2.add(pos, sameSizelineages2configs2);
                                                configSizeList.add(pos, numLin);
                                            }
                                            else{
                                                Map<Set<Integer>, List<Configuration>> sameSizelineages2configs1 = newCACs1.get(sizeIndex);
                                                sameLineageConfigs1 = sameSizelineages2configs1.get(newConfig._lineages);
                                                if(sameLineageConfigs1==null){
                                                    sameLineageConfigs1 = new ArrayList<Configuration>();
                                                    sameSizelineages2configs1.put(newConfig._lineages,sameLineageConfigs1);
                                                }

                                                Map<Set<Integer>, List<Configuration>> sameSizelineages2configs2 = newCACs2.get(sizeIndex);
                                                sameLineageConfigs2 = sameSizelineages2configs2.get(newConfig._lineages);
                                                if(sameLineageConfigs2==null){
                                                    sameLineageConfigs2 = new ArrayList<Configuration>();
                                                    sameSizelineages2configs2.put(newConfig._lineages,sameLineageConfigs2);
                                                }
                                            }

                                            Configuration copy = new Configuration(newConfig);
                                            if(config.getLineageCount()==0){
                                                copy.addNetNodeChoice(netNodeIndex, configIndex);
                                            }else{
                                                copy.addNetNodeChoice(netNodeIndex, configIndex+1);
                                            }

                                            if(k==0){
                                                sameLineageConfigs1.add(newConfig);
                                                sameLineageConfigs2.add(copy);

                                            }
                                            else{
                                                sameLineageConfigs2.add(newConfig);
                                                sameLineageConfigs1.add(copy);

                                            }
                                            splitedConfigs.add(newConfig);
                                            splitedConfigs.add(copy);
                                        }
                                        if(config.getLineageCount()==0){
                                            configIndex++;
                                        }else{
                                            configIndex = configIndex + 2;
                                        }                                    }
                                }
                            }
                        }
                    }

                    if(_printDetails){
                        System.out.print("CAC after: {");
                        for(Map<Set<Integer>,List<Configuration>> lineages2configs: newCACs1){
                            for(List<Configuration> configList: lineages2configs.values()){
                                for(Configuration config: configList){
                                    System.out.print(config.toString(gtClusters)+"  ");
                                }
                            }
                        }
                        System.out.println("}");
                        System.out.print("CAC after: {");
                        for(Map<Set<Integer>,List<Configuration>> lineages2configs: newCACs2){
                            for(List<Configuration> configList: lineages2configs.values()){
                                for(Configuration config: configList){
                                    System.out.print(config.toString(gtClusters)+"  ");
                                }
                            }
                        }
                        System.out.println("}");
                    }
                    Iterator<NetNode> parentIt = node.getParents().iterator();
                    NetNode parentNode1 = parentIt.next();
                    double distance1 = node.getParentDistance(parentNode1);
                    double hybridProb1 = node.getParentProbability(parentNode1);
                    hybridProb1 = hybridProb1==NetNode.NO_PROBABILITY?1:hybridProb1;
                    List<Configuration> temp = new ArrayList<Configuration>();
                    for(Map<Set<Integer>,List<Configuration>> lineages2configs: newCACs1){
                        for(List<Configuration> configs: lineages2configs.values()){
                            temp.addAll(configs);
                        }
                    }

                    NetNode parentNode2 = parentIt.next();
                    double distance2 = node.getParentDistance(parentNode2);
                    double hybridProb2 = node.getParentProbability(parentNode2);
                    hybridProb2 = hybridProb2==NetNode.NO_PROBABILITY?1:hybridProb2;
                    temp = new ArrayList<Configuration>();
                    for(Map<Set<Integer>,List<Configuration>> lineages2configs: newCACs2){
                        for(List<Configuration> configs: lineages2configs.values()){
                            temp.addAll(configs);
                        }
                    }

                    List<Configuration> ACminus1 = new ArrayList<Configuration>();
                    List<Configuration> ACminus2 = new ArrayList<Configuration>();
                    computeTwoACMinus(newCACs1, distance1, hybridProb1, newCACs2, distance2, hybridProb2, child2parent, R, gtClusters, ACminus1, ACminus2);

                    Map<NetNode, List<Configuration>> ACminusMap = new HashMap<>();
                    ACminusMap.put(parentNode1, ACminus1);
                    ACminusMap.put(parentNode2, ACminus2);
                    node2ACMinus.put(node, ACminusMap);

                    if(_printDetails){
                        System.out.print("ACminus to " + parentNode1.getName()+ ":  {");
                        for(Configuration config: ACminus1){
                            System.out.print(config.toString(gtClusters)+"  ");
                        }
                        System.out.println("}");
                    }

                    if(_printDetails){
                        System.out.print("ACminus to " + parentNode2.getName()+ ":  {");
                        for(Configuration config: ACminus2){
                            System.out.print(config.toString(gtClusters)+"  ");
                        }
                        System.out.println("}");
                    }
                    netNodeIndex ++;
                }


            }
            resultProbs[treeID] = gtProb;

            if(_printDetails){
                System.out.println("The probability of this gene tree is:" + gtProb);
            }
            if(_parallel){
                treeID = getNextTreeID();
            }
            else{
                treeID++;
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
     * @param ACminus	the resulting AC-, which are the lineages leaving a branch
     */
    private void computeACMinus(List<Map<Set<Integer>,List<Configuration>>> CACs, double distance, double hybridProb, Map<Integer, Integer> child2parent, boolean[][] R, List<STITreeCluster> gtClusters, List<Configuration> ACminus){
        Map<Integer, Map<Configuration, Configuration>> shape2ACminus = new HashMap<Integer, Map<Configuration, Configuration>>();
        for (Map<Set<Integer>, List<Configuration>> lineages2configs : CACs) {
            for (List<Configuration> sameLineageConfigs : lineages2configs.values()) {
                Configuration origConfig = sameLineageConfigs.get(0);
                Stack<Configuration> configStack = new Stack<Configuration>();
                Configuration configCopy = new Configuration(origConfig);
                configCopy.clearCoalEvents();
                configStack.push(configCopy);
                Map<Integer, Set<Set<Integer>>> visitedACminus = new HashMap<Integer, Set<Set<Integer>>>();

                while (!configStack.empty()) {
                    Configuration cconfig = configStack.pop();
                    double weight = calculateW(cconfig._coalEvents, R);
                    double prob = Math.max(0, computeProbability(origConfig, cconfig, weight, distance, hybridProb, gtClusters));
                    int code = cconfig._lineages.size();
                    for (int lin : cconfig._lineages) {
                        if (lin != 0) {
                            code *= lin;
                        }
                    }

                    boolean ffirst = true;
                    for(Configuration config: sameLineageConfigs){
                        double newProb = Math.max(0,prob*config._totalProb);
                        Configuration cconfigCopy;
                        if(ffirst){
                            cconfigCopy = cconfig;
                            ffirst = false;
                        }
                        else{
                            cconfigCopy = new Configuration(cconfig);
                            cconfigCopy.setNetNodeChoice(config._netNodeIndex);
                            cconfigCopy.addCoalEvents(cconfig._coalEvents);
                        }
                        cconfigCopy.setTotalProbability(newProb);
                        Map<Configuration, Configuration> cc = shape2ACminus.get(code);
                        if (cc != null) {
                            Configuration ccExisting = cc.get(cconfigCopy);
                            if (ccExisting != null) {
                                ccExisting.addTotalProbability(cconfigCopy._totalProb);
                            } else {
                                cc.put(cconfigCopy, cconfigCopy);
                            }
                        } else {
                            cc = new HashMap<Configuration, Configuration>();
                            cc.put(cconfigCopy, cconfigCopy);
                            shape2ACminus.put(code, cc);
                        }
                    }

                    if(cconfig.getLineageCount()==1)continue;
                    Map<Integer, List<Integer>> parent2children = new HashMap<Integer, List<Integer>>();
                    for (int i : cconfig._lineages) {
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
                        newConfig.mergeCluster(childrenList.get(0), childrenList.get(1), parent);
                        code = newConfig._lineages.size();
                        for (int lin : newConfig._lineages) {
                            if (lin != 0) {
                                code *= lin;
                            }
                        }
                        Set<Set<Integer>> cs = visitedACminus.get(code);
                        if (cs != null) {
                            if (cs.contains(newConfig._lineages)) {
                                continue;
                            } else {
                                cs.add(newConfig._lineages);
                            }
                        } else {
                            cs = new HashSet<Set<Integer>>();
                            cs.add(newConfig._lineages);
                            visitedACminus.put(code, cs);
                        }

                        newConfig.addCoalEvents(cconfig._coalEvents);
                        newConfig.addCoalEvent(parent);
                        configStack.push(newConfig);

                    }
                }
            }
        }

        for (Map<Configuration, Configuration> cc : shape2ACminus.values()) {
            ACminus.addAll(cc.values());
        }

    }


    /**
     * Computes the two AC- given the two CACs at a reticulation node (all possible lineages along with their probabilities when leaving a branch given all possible lineages entering a branch)
     * Assumes the two reticulation edges incident with the reticulation node are called edge #1 and edge #2 respectively.
     *
     * @param CACs1 	the given CACs on edge #1, which are the lineages entering a branch
     * @param distance1     the length of edge #1
     * @param hybridProb1	the inheritance probability of edge #1
     * @param CACs2 	the given CACs on edge #2, which are the lineages entering a branch
     * @param distance2     the length of edge #2
     * @param hybridProb2	the inheritance probability of edge #2
     * @param child2parent  maps from a child cluster/node to its parent cluster/node (induced from the gene tree topology)
     * @param R     a matrix that keeps the ancestral relationships of all the nodes in a gene tree
     * @param gtClusters    all clusters in the gene tree
     * @param ACminus1	the resulting AC- on edge #1, which are the lineages leaving a branch
     */
    private void computeTwoACMinus(List<Map<Set<Integer>,List<Configuration>>>CACs1, double distance1, double hybridProb1, List<Map<Set<Integer>,List<Configuration>>>CACs2, double distance2, double hybridProb2,Map<Integer, Integer> child2parent, boolean[][] R, List<STITreeCluster> gtClusters, List<Configuration> ACminus1, List<Configuration> ACminus2){
        Map<Integer, Map<Configuration, Configuration>> shape2ACminus1 = new HashMap<Integer, Map<Configuration, Configuration>>();
        Map<Integer, Map<Configuration, Configuration>> shape2ACminus2 = new HashMap<Integer, Map<Configuration, Configuration>>();
        Iterator<Map<Set<Integer>,List<Configuration>>> cacListIt1 = CACs1.iterator();
        Iterator<Map<Set<Integer>,List<Configuration>>> cacListIt2 = CACs2.iterator();
        while(cacListIt1.hasNext()){
            Iterator<List<Configuration>> configIt1 = cacListIt1.next().values().iterator();
            Iterator<List<Configuration>> configIt2 = cacListIt2.next().values().iterator();
            while(configIt1.hasNext()){
                List<Configuration> sameLineageConfigs1 = configIt1.next();
                List<Configuration> sameLineageConfigs2 = configIt2.next();
                Configuration origConfig1 = sameLineageConfigs1.get(0);

                Stack<Configuration> configStack = new Stack<Configuration>();
                Configuration configCopy = new Configuration(origConfig1);
                configCopy.clearCoalEvents();
                configStack.push(configCopy);
                Map<Integer, Set<Set<Integer>>> visitedACminus = new HashMap<Integer, Set<Set<Integer>>>();

                while(!configStack.empty()){
                    Configuration cconfig = configStack.pop();
                    int code = cconfig._lineages.size();
                    for(int lin: cconfig._lineages){
                        if(lin!=0){
                            code *= lin;
                        }
                    }
                    String[] forPrint = new String[2];
                    double weight = calculateW(cconfig._coalEvents, R);
                    double probCommon = computeProbabilityPart1(origConfig1, cconfig, weight, forPrint, gtClusters);
                    double prob = Math.max(0, computeProbabilityPart2(probCommon, origConfig1.getLineageCount(), cconfig.getLineageCount(), distance1, hybridProb1, forPrint));

                    boolean ffirst = true;
                    for(Configuration config: sameLineageConfigs1){
                        double newProb = Math.max(0,prob*config._totalProb);
                        Configuration cconfigCopy;
                        if(ffirst){
                            cconfigCopy = cconfig;
                            ffirst = false;
                        }
                        else{
                            cconfigCopy = new Configuration(cconfig);
                            cconfigCopy.setNetNodeChoice(config._netNodeIndex);
                        }
                        cconfigCopy.setTotalProbability(newProb);
                        Map<Configuration, Configuration> cc = shape2ACminus1.get(code);
                        if (cc != null) {
                            Configuration existing = cc.get(cconfigCopy);
                            if (existing != null) {
                                existing.addTotalProbability(cconfigCopy._totalProb);
                            } else {
                                cc.put(cconfigCopy, cconfigCopy);
                            }
                        } else {
                            cc = new HashMap<Configuration, Configuration>();
                            cc.put(cconfigCopy, cconfigCopy);
                            shape2ACminus1.put(code, cc);
                        }
                    }
                    prob = Math.max(0, computeProbabilityPart2(probCommon, origConfig1.getLineageCount(), cconfig.getLineageCount(), distance2, hybridProb2, forPrint));

                    for(Configuration config: sameLineageConfigs2){
                        double newProb = Math.max(0,prob*config._totalProb);
                        Configuration cconfigCopy = new Configuration(cconfig);
                        cconfigCopy.setNetNodeChoice(config._netNodeIndex);
                        cconfigCopy.setTotalProbability(newProb);
                        Map<Configuration, Configuration> cc = shape2ACminus2.get(code);

                        if (cc != null) {
                            Configuration existing = cc.get(cconfigCopy);
                            if (existing != null) {
                                cc.get(cconfigCopy).addTotalProbability(cconfigCopy._totalProb);
                            } else {
                                cc.put(cconfigCopy, cconfigCopy);
                            }
                        } else {
                            cc = new HashMap<Configuration, Configuration>();
                            cc.put(cconfigCopy, cconfigCopy);
                            shape2ACminus2.put(code, cc);
                        }
                    }
                    Map<Integer, List<Integer>> parent2children = new HashMap<Integer, List<Integer>>();
                    for(int i: cconfig._lineages)
                    {
                        if(child2parent.get(i)==null){
                            continue;
                        }
                        int pid = child2parent.get(i);
                        List<Integer> children = parent2children.get(pid);
                        if(children == null){
                            children = new ArrayList<Integer>();
                            parent2children.put(pid,children);
                        }
                        children.add(i);
                    }

                    for(Map.Entry<Integer, List<Integer>> entry: parent2children.entrySet())
                    {
                        Integer parent = entry.getKey();
                        List<Integer> childrenList = entry.getValue();
                        if( childrenList.size() <= 1 ){
                            continue;
                        }
                        Configuration newConfig = new Configuration(cconfig);
                        newConfig.mergeCluster(childrenList.get(0),childrenList.get(1),parent);
                        code = newConfig._lineages.size();
                        for(int lin: newConfig._lineages){
                            if(lin!=0){
                                code *= lin;
                            }
                        }
                        Set<Set<Integer>> cs = visitedACminus.get(code);
                        if(cs!=null){
                            if( cs.contains(newConfig._lineages)){
                                continue;
                            }
                            else{
                                cs.add(newConfig._lineages) ;
                            }
                        }
                        else{
                            cs = new HashSet<Set<Integer>>();
                            cs.add(newConfig._lineages);
                            visitedACminus.put(code,cs);
                        }

                        newConfig.addCoalEvents(cconfig._coalEvents);
                        newConfig.addCoalEvent(parent);
                        configStack.push(newConfig);
                    }
                }
            }
        }
        for(Map<Configuration,Configuration> cc: shape2ACminus1.values()){
            ACminus1.addAll(cc.values());
        }
        for(Map<Configuration,Configuration> cc: shape2ACminus2.values()){
            ACminus2.addAll(cc.values());
        }
    }


    /**
     * This function is to pre-process a network, including removing binary nodes, getting information of nodes
     */
    private void processNetwork(Network<Object> net, boolean fromScratch){
        if(fromScratch){
            removeBinaryNodes(net);
        }
        _netNodeNum = 0;
        int totalNode = 0;
        _node2ID = new HashMap<NetNode, Integer>();
        List<String> taxa = new ArrayList<String>();
        for(NetNode node: Networks.postTraversal(net)){
            _node2ID.put(node, totalNode++);
            if(node.isLeaf()){
                taxa.add(node.getName());
            }else if(node.isNetworkNode()){
                _netNodeNum++;
            }
        }
        _totalCoverNodes = Networks.getLowestArticulationNodes((Network)net);

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
     * This function is to remove binary nodes of a network
     */
    public static void removeBinaryNodes(Network net)
    {
        // Find all binary nodes.
        List<NetNode> binaryNodes = new ArrayList<NetNode>();
        List<NetNode> degreeTwoNodes = new ArrayList<NetNode>();
        Network<Object> network = net;
        for (NetNode node : network.bfs()) {
            if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                binaryNodes.add(node);
            }

            if (node.getIndeg() == 2 && node.getOutdeg() == 2) {
                degreeTwoNodes.add(node);
            }
        }

        // Remove them.
        for (NetNode<Object> node : binaryNodes) {
            NetNode child = node.getChildren().iterator().next();	// Node's only child.
            NetNode parent = node.getParents().iterator().next();	// Node's only parent.
            double distance = node.getParentDistance(parent) + child.getParentDistance(node);
            double inheritanceProb1 = node.getParentProbability(parent);
            inheritanceProb1 = inheritanceProb1==NetNode.NO_PROBABILITY?1:inheritanceProb1;
            double inheritanceProb2 = child.getParentProbability(node);
            inheritanceProb2 = inheritanceProb2==NetNode.NO_PROBABILITY?1:inheritanceProb2;
            double gamma = inheritanceProb1 * inheritanceProb2;
            parent.removeChild(node);
            node.removeChild(child);
            parent.adoptChild(child, distance);
            child.setParentProbability(parent, gamma);
        }
        for (NetNode<Object> node : degreeTwoNodes) {
            NetNode newNode = new BniNetNode();
            List<NetNode> removedChildren = new ArrayList<NetNode>();
            for(NetNode child: node.getChildren()){
                newNode.adoptChild(child, child.getParentDistance(node));
                child.setParentProbability(node, 1);
                removedChildren.add(child);
            }
            for(NetNode child: removedChildren){
                node.removeChild(child);
            }
            node.adoptChild(newNode, 0);

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
        double prob = 1;
        int u = preConfig.getLineageCount();
        int v = coalescedConfig.getLineageCount();
        prob = Math.pow(portion, u);
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
     * Computes the probability of a given ancestral configuration coalescing into another given one - part1
     * Splitting this into two parts is for computing ACMinus at a reticulation node
     *
     * @param preConfig 	the original ancestral configuration
     * @param coalescedConfig     the ancestral configuration that the original one is coalescing into
     * @param w     weight
     * @param forPrint 	the information for printing
     * @param gtClusters    all clusters in the gene tree
     *
     * @return the probability
     */
    private double computeProbabilityPart1(Configuration preConfig, Configuration coalescedConfig, double w, String[] forPrint, List<STITreeCluster> gtClusters){
        double prob;
        forPrint[0] = preConfig.toString(gtClusters)+"->"+coalescedConfig.toString(gtClusters)+ ": ";
        int u = preConfig.getLineageCount();
        int v = coalescedConfig.getLineageCount();
        if(u == v && (u == 1 || u == 0)){
            return 1;
        }
        double d = calculateD(u,u-v);
        prob = w/d;
        forPrint[1] = "*"+w+"/"+d;
        return prob;
    }


    /**
     * Computes the probability of a given ancestral configuration coalescing into another given one - part2
     * Splitting this into two parts is for computing ACMinus at a reticulation node
     *
     * @param prob 	the probability returned from part1
     * @param u the number of extra lineages entering a branch
     * @param v the number of extra lineages leaving a branch
     * @param distance 	the branch length
     * @param hybridProb	the inheritance probability
     * @param forPrint 	the information for printing
     *
     * @return the probability
     */
    private double computeProbabilityPart2(double prob, int u, int v, double distance, double hybridProb, String[] forPrint){
        prob *= Math.pow(hybridProb, u);
        if(u == v && (u == 1 || u == 0)){
            return prob;
        }

        if(distance == 0 && u != v){
            if(_printDetails){
                System.out.println(forPrint[0] + "g"+u+v+"(0)=0");
            }
            return 0;
        }

        double gij = gij(distance, u, v);
        if(gij<=0){
            return -1;
        }

        if(_printDetails){
            if(hybridProb==1){
                System.out.println(forPrint[0] + "g"+u+v+"("+distance+")"+forPrint[1] + "="+ gij*prob);
            }else{
                System.out.println(forPrint[0] + "g"+u+v+"("+distance+")"+forPrint[1] + "*"+hybridProb+"^"+u +"="+ gij*prob);
            }
        }
        prob *= gij;
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
     * The function is to help split lineages at a reticulation node
     */
    private List<boolean[]> getSelected(int n, int m){
        List<boolean[]> selectedList = new ArrayList<boolean[]>();
        int[] order = new int[m+1];
        for(int i=0; i<=m; i++){
            order[i] = i-1;
        }
        int k = m;
        boolean flag = true;
        while(order[0] == -1){
            if(flag){
                boolean[] bs = new boolean[n];
                for(int i=1; i<=m; i++){
                    bs[order[i]] = true;
                }
                selectedList.add(bs);
                flag = false;
            }

            order[k]++;
            if(order[k] == n){
                order[k--] = 0;
                continue;
            }

            if(k < m){
                order[++k] = order[k-1];
                continue;
            }

            if(k == m)
                flag = true;
        }

        return selectedList;
    }


    /**
     * This class is to represent the concept of ancestral configuration
     */
    private class Configuration{
        private HashSet<Integer> _lineages;
        private double _totalProb;
        int[] _netNodeIndex;
        private Set<Integer> _coalEvents;


        public Configuration(){
            _lineages = new HashSet<Integer>();
            _netNodeIndex = new int[_netNodeNum];
            Arrays.fill(_netNodeIndex, 0);
        }

        public Configuration(Configuration config){
            _lineages = (HashSet)config._lineages.clone();
            _totalProb = config._totalProb;
            _netNodeIndex = config._netNodeIndex.clone();
        }

        public Configuration(Configuration config1, Configuration config2){
            _lineages = (HashSet)config1._lineages.clone();
            _lineages.addAll(config2._lineages);
            _totalProb = Math.max(0,config1._totalProb*config2._totalProb);
            _netNodeIndex = new int[_netNodeNum];
            for(int i=0; i< _netNodeNum; i++){
                if(config1._netNodeIndex[i] == config2._netNodeIndex[i]){
                    _netNodeIndex[i] = config1._netNodeIndex[i];
                }
                else{
                    _netNodeIndex[i] = Math.max(config1._netNodeIndex[i], config2._netNodeIndex[i]);
                }
            }
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


        public void addLineage(int index){
            _lineages.add(index);
        }


        public void mergeCluster(int child1, int child2, int parent){
            _lineages.remove(child1);
            _lineages.remove(child2);
            _lineages.add(parent);
        }

        public void setTotalProbability(double prob){
            _totalProb = prob;
        }


        public int getLineageCount(){
            return _lineages.size();
        }

        public void addTotalProbability(double adds){
            _totalProb += adds;
        }

        public String toString(List<STITreeCluster> gtClusters){
            String exp = "";
            for(int id: _lineages) {
                exp = exp + gtClusters.get(id);
            }
            exp = exp + "/[";
            for(int i=0; i<_netNodeIndex.length; i++){
                exp = exp + _netNodeIndex[i];
                if(i!=_netNodeIndex.length-1){
                    exp = exp + ",";
                }
            }
            exp = exp + "]:" + _totalProb;
            return exp;
        }

        public void addNetNodeChoice(int net, int index){
            _netNodeIndex[net] = index;
        }

        public void setNetNodeChoice(int[] choice){
            _netNodeIndex = choice.clone();
        }

        public void clearNetNodeChoice(){
            Arrays.fill(_netNodeIndex, 0);
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
            return config._lineages.equals(_lineages) && Arrays.equals(config._netNodeIndex, _netNodeIndex);
        }

        public int hashCode(){
            return _lineages.hashCode()+Arrays.hashCode(_netNodeIndex);
        }

    }



}
