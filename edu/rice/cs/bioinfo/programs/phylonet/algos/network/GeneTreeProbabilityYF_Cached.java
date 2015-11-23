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
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 5/10/12
 * Time: 2:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneTreeProbabilityYF_Cached {
    Set<NetNode> _totalCoverNodes;
    boolean _printDetails = false;
    int _netNodeNum;
    //List<STITreeCluster> _gtClusters;
    boolean[][] _M;
    Map<NetNode, Integer> _node2ID;
    boolean _parallel = false;
    int _currentTreeID = 0;
    int _totalTree;
    boolean _preProcessed = false;

    int _maxNumACs = 0;


    public void setParallel(boolean parallel){
        _parallel = parallel;
    }

    public void setTotalNodes(Set<NetNode> totalCoverNodes){
        _totalCoverNodes = totalCoverNodes;
    }

    public void setPrintDetails(boolean p){
        _printDetails = p;
    }

    public synchronized int getNextTreeID(){
        return _currentTreeID++;
    }

    public void preProcess(Network network, List<Tree> gts, boolean fromScratch){
        _totalTree = gts.size();
        processNetwork(network, _totalTree, fromScratch);
        _preProcessed = true;
    }


    public void calculateGTDistribution(Network<CoalescePattern[]> network, List<Tree> gts, Map<String, List<String>> species2alleles, double resultProbs[]){
        if(!_preProcessed){
            preProcess(network, gts, true);
        }

        int treeID = 0;
        if(_parallel){
            treeID = getNextTreeID();
        }


        //List<Long> runningTime = new ArrayList<>();
        //System.out.println(Arrays.toString(resultProbs));
        while(treeID < _totalTree){
            //long startTime = System.currentTimeMillis();
            //System.out.println("Computing tree" + treeID);
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

            for(NetNode<CoalescePattern[]> node: Networks.postTraversal(network)){
                //long start = System.currentTimeMillis();
                CoalescePattern cp = new CoalescePattern();
                node.getData()[treeID] = cp;
                //t6 += (System.currentTimeMillis()-start)/1000.0;
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
                    //start = System.currentTimeMillis();
                    cp.addACs(node.getParents().iterator().next(), tempList);
                    //t5 += (System.currentTimeMillis()-start)/1000.0;
                }
                else{
                    if(node.isNetworkNode()){
                        CoalescePattern childCP = ((NetNode<CoalescePattern[]>)(node.getChildren().iterator().next())).getData()[treeID];
                        Map<Set<Integer>,List<Configuration>> temp = new HashMap<Set<Integer>, List<Configuration>>();
                        List<Configuration> ACMinus = childCP._ACMinuss.get(node);
                        temp.put(null,ACMinus);
                        CACs.add(temp);
                        //cp.addACs(null, ACMinus);
                        /*
                        List<Integer> configSizeList = new ArrayList<Integer>();
                        for(Configuration config: edge2ACminus.remove(edge)){
                            int numLin = config._lineages.size();
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
                                sameSizelineages2configs.put(config._lineages,sameLineageConfigs);
                                CACs.add(pos, sameSizelineages2configs);
                                configSizeList.add(pos, numLin);
                            }
                            else{
                                Map<Set<Integer>, List<Configuration>> sameSizelineages2configs = CACs.get(sizeIndex);
                                sameLineageConfigs = sameSizelineages2configs.get(config._lineages);
                                if(sameLineageConfigs==null){
                                    sameLineageConfigs = new ArrayList<Configuration>();
                                    sameSizelineages2configs.put(config._lineages,sameLineageConfigs);
                                }
                            }
                            sameLineageConfigs.add(config);

                        }
                        */

                    }
                    else{
                        Iterator<NetNode<CoalescePattern[]>> childNode = node.getChildren().iterator();
                        List<Configuration> AC1 = childNode.next().getData()[treeID].getACMinuss(node);
                        List<Configuration> AC2 = childNode.next().getData()[treeID].getACMinuss(node);

                        List<Integer> configSizeList = new ArrayList<Integer>();
                        boolean total = _totalCoverNodes.contains(node);
                        List<Configuration> mergedConfigList = new ArrayList<Configuration>();
                        for(Configuration config1: AC1){
                            for(Configuration config2: AC2){
                                if(config1.isCompatible(config2)){
                                    Configuration mergedConfig = new Configuration(config1, config2);
                                    //mergedConfigs.add(mergedConfig);

                                    /*
                                    if(mergedConfig._totalProb <=0 ){
                                        continue;
                                    }
                                    */
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
                                            //start = System.currentTimeMillis();
                                            sameLineageConfigs.get(0).addChildPair(mergedConfig._childPairList.get(0));
                                            //t2 += (System.currentTimeMillis()-start)/1000.0;
                                        }
                                    }
                                    else{
                                        sameLineageConfigs.add(mergedConfig);
                                        mergedConfigList.add(mergedConfig);
                                    }

                                }

                            }
                        }
                        if(node.isRoot()){
                            //start = System.currentTimeMillis();
                            cp.addACs(null, mergedConfigList);
                            //t5 += (System.currentTimeMillis()-start)/1000.0;
                        }
                        else{
                            //start = System.currentTimeMillis();
                            cp.addACs(node.getParents().iterator().next(), mergedConfigList);
                            //t5 += (System.currentTimeMillis()-start)/1000.0;
                        }
                    }
                }

                int count = 0;
                for(Map<Set<Integer>,List<Configuration>> lineages2configs: CACs){
                    for(List<Configuration> configList: lineages2configs.values()){
                        count += configList.size();
                    }
                }
                _maxNumACs = Math.max(count, _maxNumACs);


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
                                    //start = System.currentTimeMillis();
                                    rootConfig.addUncoalescedConfiguration(preConfig, 1, 1);
                                    //t1 += (System.currentTimeMillis()-start)/1000.0;
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
                                    //start = System.currentTimeMillis();
                                    rootConfig.addUncoalescedConfiguration(preConfig, weight, prob);
                                    //t1 += (System.currentTimeMillis()-start)/1000.0;
                                }

                            }
                        }

                    }
                    List<Configuration> temp = new ArrayList<Configuration>();
                    temp.add(rootConfig);
                    //start = System.currentTimeMillis();
                    cp.addACMinuss(null, temp);
                    //t4 += (System.currentTimeMillis()-start)/1000.0;
                }
                else if(node.isTreeNode()){
                    double distance = node.getParentDistance(node.getParents().iterator().next());
                    List<Configuration> ACminus = new ArrayList<Configuration>();
                    computeACMinus(CACs, distance, 1, child2parent, R, gtClusters, ACminus);
                    //start = System.currentTimeMillis();
                    cp.addACMinuss(node.getParents().iterator().next(),ACminus);

                    _maxNumACs = Math.max(_maxNumACs, ACminus.size());
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
                    //start = System.currentTimeMillis();
                    cp.setConfig2splitedConfigs(config2splitedConfigs);
                    //t3 += (System.currentTimeMillis()-start)/1000.0;
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
                    Iterator<NetNode<CoalescePattern[]>> it = node.getParents().iterator();
                    NetNode<CoalescePattern[]> parentNode1 = it.next();
                    double distance1 = node.getParentDistance(parentNode1);
                    double hybridProb1 = node.getParentProbability(parentNode1);
                    hybridProb1 = hybridProb1==NetNode.NO_PROBABILITY?1:hybridProb1;
                    List<Configuration> temp = new ArrayList<Configuration>();
                    for(Map<Set<Integer>,List<Configuration>> lineages2configs: newCACs1){
                        for(List<Configuration> configs: lineages2configs.values()){
                            temp.addAll(configs);
                        }
                    }
                    //start = System.currentTimeMillis();
                    cp.addACs(parentNode1, temp);
                    //t5 += (System.currentTimeMillis()-start)/1000.0;
                    NetNode<CoalescePattern[]> parentNode2 = it.next();
                    double distance2 = node.getParentDistance(parentNode2);
                    double hybridProb2 = node.getParentProbability(parentNode2);
                    hybridProb2 = hybridProb2==NetNode.NO_PROBABILITY?1:hybridProb2;
                    temp = new ArrayList<Configuration>();
                    for(Map<Set<Integer>,List<Configuration>> lineages2configs: newCACs2){
                        for(List<Configuration> configs: lineages2configs.values()){
                            temp.addAll(configs);
                        }
                    }
                    //start = System.currentTimeMillis();
                    cp.addACs(parentNode2, temp);
                    //t5 += (System.currentTimeMillis()-start)/1000.0;
                    List<Configuration> ACminus1 = new ArrayList<Configuration>();
                    List<Configuration> ACminus2 = new ArrayList<Configuration>();
                    computeTwoACMinus(newCACs1, distance1, hybridProb1, newCACs2, distance2, hybridProb2,child2parent, R, gtClusters, ACminus1, ACminus2);
                    //start = System.currentTimeMillis();
                    cp.addACMinuss(parentNode1, ACminus1);
                    //t4 += (System.currentTimeMillis()-start)/1000.0;
                    if(_printDetails){
                        System.out.print("ACminus to " + parentNode1.getName()+ ":  {");
                        for(Configuration config: ACminus1){
                            System.out.print(config.toString(gtClusters)+"  ");
                        }
                        System.out.println("}");
                    }


                    _maxNumACs = Math.max(_maxNumACs, ACminus1.size());
                    //start = System.currentTimeMillis();
                    cp.addACMinuss(parentNode2, ACminus2);
                    //t4 += (System.currentTimeMillis()-start)/1000.0;
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
            //System.out.println(Thread.currentThread().getId() + " ends #" + treeID);
            if(_parallel){
                treeID = getNextTreeID();
            }
            else{
                treeID++;
            }

            //runningTime.add(System.currentTimeMillis()-startTime);

        }

        //System.out.println(runningTime);
        //System.out.println(Arrays.toString(resultProbs));

    }


    public void calculateGTDistribution(Network<CoalescePattern[]> network, List<Tree> gts, Set<NetNode> editedChilds, Set<NetNode> editedParents, double[] probs){
        if(!_preProcessed){
            preProcess(network, gts, false);
        }

        //int childID = _node2ID.get(editedChild);
        //int parentID = editedParent==null ? childID : _node2ID.get(editedParent);
        Set<Integer> childIDSet = new HashSet<Integer>();
        for(NetNode editedChild: editedChilds){
            childIDSet.add(_node2ID.get(editedChild));
        }
        Set<Integer> parentIDSet = new HashSet<Integer>();
        for(NetNode editedParent: editedParents){
            parentIDSet.add(_node2ID.get(editedParent));
        }


        int treeID = 0;
        if(_parallel){
            treeID = getNextTreeID();
        }

        while(treeID < _totalTree){
            //System.out.println("Updating tree" + treeID);
            double gtProb = 0;
            /*

            String[] gtTaxa = gt.getLeaves();
            Map<Integer,Integer> child2parent = new HashMap<Integer, Integer>();
            processGT(gt, gtTaxa, child2parent);
            computeR();
            */
            List<STITreeCluster> gtClusters = null;
            if(_printDetails){
                Tree gt = gts.get(treeID);
                String[] gtTaxa = gt.getLeaves();
                gt.getClusters(gtTaxa, false);
            }
            for(NetNode<CoalescePattern[]> node: Networks.postTraversal(network)){
                CoalescePattern cp = node.getData()[treeID];
                if(_printDetails){
                    System.out.println();
                    System.out.println("On node #" + _node2ID.get(node) + " " + node.getName());
                }
                int nodeID = _node2ID.get(node);
                boolean unchange = !childIDSet.contains(nodeID);
                if(unchange){
                    //System.out.println(parentIDSet);
                    for(int editedParentID: parentIDSet){
                        if(_M[nodeID][editedParentID]){
                            unchange = false;
                            break;
                        }
                    }
                }

                if(unchange){
                    if(_printDetails){
                        System.out.println("Unchange ");

                        System.out.println("AC:");
                        for(Map.Entry<NetNode, List<Configuration>> node2configList: cp._ACs.entrySet()){
                            System.out.println("To " + node2configList.getKey().getName());
                            for(Configuration config: node2configList.getValue()){
                                System.out.print(config.toString(gtClusters)+"  ");
                            }
                            System.out.println();
                        }
                        System.out.println();

                        System.out.println("AC Minus:");
                        for(Map.Entry<NetNode, List<Configuration>> node2configList: cp._ACMinuss.entrySet()){
                            System.out.println("To " + node2configList.getKey().getName());
                            for(Configuration config: node2configList.getValue()){
                                System.out.print(config.toString(gtClusters)+"  ");
                            }
                            System.out.println();
                        }

                    }
                    continue;
                }

                boolean ACUpdate = false;
                for(NetNode child: node.getChildren()){
                    if(childIDSet.contains(_node2ID.get(child))){
                        ACUpdate = true;
                        break;
                    }
                }

                if(ACUpdate || !childIDSet.contains(nodeID)){
                    if(node.isNetworkNode()){
                        for(Map.Entry<Configuration, List<Configuration>> entry: cp._config2splitedConfigs.entrySet()){
                            double prob = Math.sqrt(entry.getKey()._totalProb);
                            for(Configuration config: entry.getValue()){
                                config.setTotalProbability(prob);
                            }
                        }
                    }else if(node.isTreeNode()){
                        NetNode parent;
                        if(node.isRoot()){
                            parent = null;
                        }
                        else{
                            parent = node.getParents().iterator().next();
                        }
                        for(Configuration config: cp.getACs(parent)){
                            config.setTotalProbability(config.getChildrenProbProduction());
                        }
                    }

                    if(_printDetails){
                        System.out.println("AC:");
                        for(Map.Entry<NetNode, List<Configuration>> node2configList: cp._ACs.entrySet()){
                            if(node2configList.getKey()!=null)
                                System.out.println("To " + node2configList.getKey().getName());
                            for(Configuration config: node2configList.getValue()){
                                System.out.print(config.toString(gtClusters)+"  ");
                            }
                            System.out.println();
                        }
                        System.out.println();
                    }

                    if(!node.isRoot()){
                        for(NetNode<CoalescePattern[]> parent: node.getParents()){
                            for(Configuration coalescedConfig: cp.getACMinuss(parent)){
                                coalescedConfig._totalProb = 0;
                                for(CoalescingInfo info: coalescedConfig._coalesingInfo){
                                    coalescedConfig.addTotalProbability(Math.max(0,info._uncoalescedConfig._totalProb*info._coalProb));
                                }
                            }
                        }
                    }

                    if(_printDetails){
                        System.out.println("AC Minus:");
                        for(Map.Entry<NetNode, List<Configuration>> node2configList: cp._ACMinuss.entrySet()){
                            if(node2configList.getKey()!=null)
                                System.out.println("To " + node2configList.getKey().getName());
                            for(Configuration config: node2configList.getValue()){
                                System.out.print(config.toString(gtClusters)+"  ");
                            }
                            System.out.println();
                        }
                    }
                }

                if(childIDSet.contains(nodeID)){
                    if(_printDetails){
                        System.out.println("AC:");
                        for(Map.Entry<NetNode, List<Configuration>> node2configList: cp._ACs.entrySet()){
                            System.out.println("To " + node2configList.getKey().getName());
                            for(Configuration config: node2configList.getValue()){
                                System.out.print(config.toString(gtClusters)+"  ");
                            }
                            System.out.println();
                        }
                        System.out.println();
                    }
                    for(NetNode parent: node.getParents()){
                        if(!editedParents.contains(parent))continue;
                        for(Configuration coalescedConfig: cp.getACMinuss(parent)){
                            coalescedConfig._totalProb = 0;
                            Iterator<CoalescingInfo> coalescedInfos = coalescedConfig._coalesingInfo.iterator();
                            while(coalescedInfos.hasNext()){
                                CoalescingInfo info = coalescedInfos.next();
                                double inheritanceProb = node.getParentProbability(parent);
                                inheritanceProb = inheritanceProb==NetNode.NO_PROBABILITY?1:inheritanceProb;
                                double prob = Math.max(0,computeProbability(info._uncoalescedConfig, coalescedConfig, info._coalWeight, node.getParentDistance(parent), inheritanceProb, gtClusters));
                                info._coalProb = prob;
                                coalescedConfig.addTotalProbability(Math.max(0, prob*info._uncoalescedConfig._totalProb));
                            }
                        }
                    }
                    if(_printDetails){
                        System.out.println("AC Minus:");
                        for(Map.Entry<NetNode, List<Configuration>> node2configList: cp._ACMinuss.entrySet()){
                            System.out.println("To " + node2configList.getKey().getName());
                            for(Configuration config: node2configList.getValue()){
                                System.out.print(config.toString(gtClusters)+"  ");
                            }
                            System.out.println();
                        }
                    }
                }


                if(node.isRoot()){
                    Configuration rootConfig = cp.getACMinuss(null).get(0);
                    Iterator<CoalescingInfo> coalescedInfos = rootConfig._coalesingInfo.iterator();
                    while(coalescedInfos.hasNext()){
                        CoalescingInfo info = coalescedInfos.next();
                        gtProb += Math.max(0, info._uncoalescedConfig._totalProb*info._coalProb);
                    }
                }
            }
            probs[treeID] = gtProb;
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
                    //boolean allNegative = true;
                    double weight = calculateW(cconfig._coalEvents, R);
                    double prob = Math.max(0, computeProbability(origConfig, cconfig, weight, distance, hybridProb, gtClusters));
                    /*
                    cconfig._prob = prob * origConfig._prob;
                    if (cconfig._prob > 0) {
                        allNegative = false;
                    }
                    */
                    int code = cconfig._lineages.size();
                    for (int lin : cconfig._lineages) {
                        if (lin != 0) {
                            code *= lin;
                        }
                    }

                    boolean ffirst = true;
                    for(Configuration config: sameLineageConfigs){
                        double newProb = Math.max(0,prob*config._totalProb);
                        /*
                        if(newProb<=0){
                            continue;
                        }
                        allNegative = false;
                        */
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
                                //long start = System.currentTimeMillis();
                                ccExisting.addUncoalescedConfiguration(config, weight, prob);
                                //t1 += (System.currentTimeMillis()-start)/1000.0;
                            } else {
                                cc.put(cconfigCopy, cconfigCopy);
                                //long start = System.currentTimeMillis();
                                cconfigCopy.addUncoalescedConfiguration(config, weight, prob);
                                //t1 += (System.currentTimeMillis()-start)/1000.0;
                            }
                        } else {
                            cc = new HashMap<Configuration, Configuration>();
                            cc.put(cconfigCopy, cconfigCopy);
                            shape2ACminus.put(code, cc);
                            //long start = System.currentTimeMillis();
                            cconfigCopy.addUncoalescedConfiguration(config, weight, prob);
                            //t1 += (System.currentTimeMillis()-start)/1000.0;
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
                //Configuration origConfig2 = sameLineageConfigs2.get(0);

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
                    //boolean allNegative = true;
                    //double prob = computeProbability(origConfig1, cconfig, distance1, hybridProb1);
                    String[] forPrint = new String[2];
                    double weight = calculateW(cconfig._coalEvents, R);
                    double probCommon = computeProbabilityPart1(origConfig1, cconfig, weight, forPrint, gtClusters);
                    double prob = Math.max(0, computeProbabilityPart2(probCommon, origConfig1.getLineageCount(), cconfig.getLineageCount(), distance1, hybridProb1, forPrint));

                    boolean ffirst = true;
                    for(Configuration config: sameLineageConfigs1){
                        double newProb = Math.max(0,prob*config._totalProb);
                        /*
                        if(newProb<=0){
                            continue;
                        }
                        allNegative = false;
                        */
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
                                //long start = System.currentTimeMillis();
                                existing.addUncoalescedConfiguration(config, weight, prob);
                                //t1 += (System.currentTimeMillis()-start)/1000.0;
                            } else {
                                cc.put(cconfigCopy, cconfigCopy);
                                //long start = System.currentTimeMillis();
                                cconfigCopy.addUncoalescedConfiguration(config, weight, prob);
                                //t1 += (System.currentTimeMillis()-start)/1000.0;
                            }
                        } else {
                            cc = new HashMap<Configuration, Configuration>();
                            cc.put(cconfigCopy, cconfigCopy);
                            shape2ACminus1.put(code, cc);
                            //long start = System.currentTimeMillis();
                            cconfigCopy.addUncoalescedConfiguration(config, weight, prob);
                            //t1 += (System.currentTimeMillis()-start)/1000.0;
                        }
                    }

                    prob = Math.max(0, computeProbabilityPart2(probCommon, origConfig1.getLineageCount(), cconfig.getLineageCount(), distance2, hybridProb2, forPrint));

                    //prob = computeProbability(origConfig2, cconfig, distance2, hybridProb2);
                    for(Configuration config: sameLineageConfigs2){
                        double newProb = Math.max(0,prob*config._totalProb);
                        /*
                        if(newProb<=0){
                            continue;
                        }
                        allNegative = false;
                        */
                        Configuration cconfigCopy = new Configuration(cconfig);
                        cconfigCopy.setNetNodeChoice(config._netNodeIndex);
                        cconfigCopy.setTotalProbability(newProb);
                        Map<Configuration, Configuration> cc = shape2ACminus2.get(code);

                        if (cc != null) {
                            Configuration existing = cc.get(cconfigCopy);
                            if (existing != null) {
                                cc.get(cconfigCopy).addTotalProbability(cconfigCopy._totalProb);
                                //long start = System.currentTimeMillis();
                                existing.addUncoalescedConfiguration(config, weight, prob);
                                //t1 += (System.currentTimeMillis()-start)/1000.0;
                            } else {
                                cc.put(cconfigCopy, cconfigCopy);
                                //long start = System.currentTimeMillis();
                                cconfigCopy.addUncoalescedConfiguration(config, weight, prob);
                                //t1 += (System.currentTimeMillis()-start)/1000.0;
                            }
                        } else {
                            cc = new HashMap<Configuration, Configuration>();
                            cc.put(cconfigCopy, cconfigCopy);
                            shape2ACminus2.put(code, cc);
                            //long start = System.currentTimeMillis();
                            cconfigCopy.addUncoalescedConfiguration(config, weight, prob);
                            //t1 += (System.currentTimeMillis()-start)/1000.0;
                        }
                    }
                    /*
                    if(allNegative){
                        continue;
                    }
*                   */

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

    private void processNetwork(Network<CoalescePattern[]> net, int totalNumGTs, boolean fromScratch){
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
        if(fromScratch){
            _totalCoverNodes = Networks.getLowestArticulationNodes((Network)net);
            for(NetNode<CoalescePattern[]> node: Networks.postTraversal(net)){
                node.setData(new CoalescePattern[totalNumGTs]);
            }
        }
        else{
            computeM(net, totalNode);
        }
    }


    private void computeM(Network<CoalescePattern[]> net, int numEdge){
        _M = new boolean[numEdge][numEdge];
        for(NetNode<CoalescePattern[]> node: Networks.postTraversal(net)){

            int pID = _node2ID.get(node);
            _M[pID][pID] = true;
            for(NetNode child: node.getChildren()){
                int cID = _node2ID.get(child);
                _M[pID][cID] = true;
                for(int i=0; i<numEdge; i++){
                    if(_M[cID][i]){
                        _M[pID][i] = true;
                    }
                }
            }
        }
    }


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

        for (NetNode<CoalescePattern[]> node : degreeTwoNodes) {
            NetNode<CoalescePattern[]> newNode = new BniNetNode<CoalescePattern[]>();
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


/*
    public static Set<NetNode> computeNodeCoverage(Network net){
        //List<Integer> leaves = new ArrayList<Integer>();
        Network<Object> network = net;
        List<NetNode> allTotalNodes = new ArrayList<NetNode>();
        Set<NetNode> totalCoverNodes = new HashSet<NetNode>();
        for(NetNode<Object> node: Networks.postTraversal(network)){
            if(node.isLeaf()){
                //leaves.add(id);
                allTotalNodes.add(node);
            }

            else if(node.isRoot()){
                boolean ftotal = true;
                for(NetNode child: node.getChildren()){
                    if(!allTotalNodes.contains(child)){
                        ftotal = false;
                        break;
                    }
                }
                if(!ftotal){
                    totalCoverNodes.add(node);
                }

            }
            else if(node.isTreeNode()){
                boolean ftotal = true;
                for(NetNode child: node.getChildren()){
                    if(!allTotalNodes.contains(child)){
                        ftotal = false;
                        break;
                    }
                }
                if(ftotal){
                    allTotalNodes.add(node);
                }else{
                    NetNode parent = node.getParents().iterator().next();
                    double distance = node.getParentDistance(parent);
                    parent.removeChild(node);
                    boolean disconnect = isValidNetwork(net, parent);
                    parent.adoptChild(node, distance);
                    if (disconnect) {
                        totalCoverNodes.add(node);
                        allTotalNodes.add(node);
                    }
                }

            }
        }
        return totalCoverNodes;
    }

    private static boolean isValidNetwork(Network<Object> net, NetNode ignoreNode){
        Set<NetNode> visited = new HashSet<NetNode>();
        Set<NetNode> seen = new HashSet<NetNode>();
        for(NetNode<Object> node: net.bfs()){
            if(node.getIndeg()==1 && node.getOutdeg()==1 && node!=ignoreNode){
                return false;
            }
            visited.add(node);
            for(NetNode parent: node.getParents()){
                seen.add(parent);
            }
            for(NetNode child: node.getChildren()){
                seen.add(child);
            }
        }
        return visited.size()==seen.size();
    }
*/


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

    private double calculateD(int u, int c){
        double d = 1;
        if(c!=0){
            for(int i=1; i<=c; i++){
                d = d * chooseD(u - i + 1, 2);
            }
        }
        return d;
    }

    private double chooseD(int N, int K) {
        double ret = 1.0;
        for (int k = 0; k < K; k++) {
            //ret = ret.multiply(BigInteger.valueOf(N-k))
            //.divide(BigInteger.valueOf(k+1));
            ret = ret*((N-k+0.0)/(k+1));
        }

        return ret;
    }


    private double calculateW(Set<Integer> coalEvents, boolean[][] R){
        //System.out.println(coalEvents);
        double w = 1.0;
        w = w*fact(1, coalEvents.size());
        //System.out.print(fact(1, coalEvents.cardinality()));
        for (int i: coalEvents) {
            int a = 0;
            for (int j: coalEvents) {
                if(i!=j && R[i][j])
                    a++;
            }
            //System.out.print(" * " + 1.0/(1 + a));
            w = w*(1.0/(1 + a));
        }
        //System.out.println();
        return w;
    }


    private double calculateD(Set<Integer> coalEvents, boolean[][] R){
        double w = 1.0;
        w = w*fact(1, coalEvents.size());
        //System.out.print(fact(1, coalEvents.cardinality()));
        for (int i: coalEvents) {
            int a = 0;
            for (int j: coalEvents) {
                if(i!=j && R[i][j])
                    a++;
            }
            //System.out.print(" * " + 1.0/(1 + a));
            w = w*(1.0/(1 + a));
        }
        //System.out.println();
        return w;
    }

    private double fact(int start, int end){
        double result = 1;
        for(int i=start; i<=end; i++){
            result = result*i;
        }

        return result;
    }

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


    public class CoalescePattern{
        private Map<NetNode, List<Configuration>> _ACMinuss;
        private Map<NetNode, List<Configuration>> _ACs; //for network node
        private Map<Configuration, List<Configuration>> _config2splitedConfigs;


        public CoalescePattern(){
            _ACMinuss = new HashMap<NetNode, List<Configuration>>();
            _ACs = new HashMap<NetNode, List<Configuration>>();
        }


        public void setConfig2splitedConfigs(Map<Configuration, List<Configuration>> map){
            _config2splitedConfigs = map;
        }

        public void addACMinuss(NetNode parent, List<Configuration> ACMinuss){
            _ACMinuss.put(parent, ACMinuss);
        }

        public void addACs(NetNode parent, List<Configuration> ACs){
            /*
            List<Configuration> existings = _ACs.get(parent);
            if(existings == null){
                List<Configuration> copy = new ArrayList<Configuration>();
                copy.addAll(ACs);
                _ACs.put(parent, copy);
            }
            else{
                existings.addAll(ACs);
            }
            */
            _ACs.put(parent, ACs);
        }

        public List<Configuration> getACMinuss(NetNode parent){
            return _ACMinuss.get(parent);
        }

        public List<Configuration> getACs(NetNode parent){
            return _ACs.get(parent);
        }
    }


    private class CoalescingInfo{
        Configuration _uncoalescedConfig;
        double _coalWeight;
        double _coalProb;

        public CoalescingInfo(Configuration c, double w, double p){
            _uncoalescedConfig = c;
            _coalWeight = w;
            _coalProb = p;
        }
    }




    private class Configuration{
        private HashSet<Integer> _lineages;
        private double _totalProb;
        int[] _netNodeIndex;
        private Set<Integer> _coalEvents;
        private List<CoalescingInfo> _coalesingInfo;  //uncoalescedConfigs, coalWeights, probsToUncoalescedConfig
        private List<Configuration[]> _childPairList;

        public Configuration(){
            _lineages = new HashSet<Integer>();
            _netNodeIndex = new int[_netNodeNum];
            Arrays.fill(_netNodeIndex, 0);
            _coalesingInfo = new ArrayList<CoalescingInfo>();
        }

        public Configuration(Configuration config){
            _lineages = (HashSet)config._lineages.clone();
            _totalProb = config._totalProb;
            _netNodeIndex = config._netNodeIndex.clone();
            _coalesingInfo = new ArrayList<CoalescingInfo>();
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
            _childPairList = new ArrayList<Configuration[]>();
            Configuration[] configPair = new Configuration[2];
            configPair[0] = config1;
            configPair[1] = config2;
            _childPairList.add(configPair);
            _coalesingInfo = new ArrayList<CoalescingInfo>();
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

        public void addUncoalescedConfiguration(Configuration config, double weight, double prob){
            _coalesingInfo.add(new CoalescingInfo(config, weight, prob));
        }


        public void addChildPair(Configuration[] pair){
            _childPairList.add(pair);
        }

        public double getChildrenProbProduction(){
            double sum = 0;
            for(Configuration[] configPair: _childPairList){
                sum += Math.max(0,configPair[0]._totalProb * configPair[1]._totalProb);
            }
            return sum;
        }

        /*
        public void addChildConfigs(Configuration config1, Configuration config2){
            if(_childConfigs==null){
                _childConfigs = new ArrayList<Configuration[]>();
            }
            Configuration[] configPair = new Configuration[2];
            configPair[0] = config1;
            configPair[1] = config2;
            _childConfigs.add(configPair);
        }
        */

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
