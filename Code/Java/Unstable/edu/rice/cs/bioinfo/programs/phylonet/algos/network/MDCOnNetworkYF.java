package edu.rice.cs.bioinfo.programs.phylonet.algos.network;


import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 5/9/12
 * Time: 1:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class MDCOnNetworkYF {
    String[] _netTaxa;
    Set<Integer> _firstIndependentNodes;
    Set<Integer> _allIndependentNodes;
    boolean _printDetail = false;
    int _netNodeNum;
    int _totalNodeNum;
    double[][] _totalNetNodeLinNum;
    boolean _parallel = false;
    int _currentTreeID = -1;
    boolean _preProcessed = false;



    public void setParallel(boolean parallel){
        _parallel = parallel;
    }

    public void setIndependentNodes(Set<Integer> firstIndependentNodes, Set<Integer> allIndependentNodes){
        _firstIndependentNodes = firstIndependentNodes;
        _allIndependentNodes = allIndependentNodes;

    }


    public double[][] getNetNodeLinNum(){
        return _totalNetNodeLinNum;
    }

    public double[] getHybridProbabilities(){
        double[] probabilities = new double[_totalNetNodeLinNum.length];
        int index = 0;
        for(double[] lineageNum: _totalNetNodeLinNum){
            double total = lineageNum[0]+lineageNum[1];
            //System.out.println(lineageNum[0]+"/"+total);
            if(total == 0){
                probabilities[index] = 0;
            }
            else{
                probabilities[index] = lineageNum[0]/total;
            }
            index++;
        }
        return probabilities;
    }

    public synchronized int getNextTreeID(){
        //System.out.println("In calculation:" + Thread.currentThread().getId() + ":");
        _currentTreeID ++;
        //System.out.println("In calculation:" + Thread.currentThread().getId() + ":" + _currentTreeID);
        return _currentTreeID;
    }


    public void countExtraCoal(Network<Integer> network, List<MutableTuple<Tree,Double>> gts, Map<String, List<String>> species2alleles, int[] xls){
        if(!_preProcessed){
            processNetwork(network);
        }

        int treeID = 0;
        if(_parallel){
            treeID = getNextTreeID();
        }

        while(treeID < gts.size()){
            MutableTuple<Tree,Double> tuple = gts.get(treeID);
            Tree gt = tuple.Item1;
            double weight = tuple.Item2;
            String[] gtTaxa = gt.getLeaves();
            Map<Integer,Integer> child2parent = new HashMap<Integer, Integer>();
            Map<Integer,Integer> node2outdegree = new HashMap<Integer, Integer>();
            List<STITreeCluster> gtClusters = new ArrayList<STITreeCluster>();
            processGT(gt, gtTaxa, child2parent, node2outdegree, gtClusters);
            Map<Integer,Integer> addedNode2resolvedDegree = new HashMap<Integer, Integer>();
            int numGTNode = gtClusters.size();

            if(!checkLeafAgreement(species2alleles, gtTaxa)){
                throw new RuntimeException("Gene tree " + gt + " has leaf that the network doesn't have.");
            }

            HashSet<String> gtTaxaSet = new HashSet<String>();
            Collections.addAll(gtTaxaSet, gtTaxa);

            HashMap<BitSet, List<Configuration>> edge2ACminus = new HashMap<BitSet, List<Configuration>>();
            int netNodeIndex = 0;
            int xl = Integer.MAX_VALUE;
            for(NetNode<Integer> node: Networks.postTraversal(network)){

                if(_printDetail){
                    System.out.println();
                    System.out.println("On node #" + node.getData() + " " + node.getName());
                }
                Map<Set<Integer>,List<Configuration>> CACs = new HashMap<Set<Integer>,List<Configuration>>();

                //set AC for a node
                if(node.isLeaf()){
                    //Map<Set<Integer>,List<Configuration>> sizeOneConfigs = new HashMap<Set<Integer>, List<Configuration>>();
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
                    config.setExtraLineage(0);
                    List<Configuration> tempList = new ArrayList<Configuration>();
                    tempList.add(config);
                    //sizeOneConfigs.put(config._lineages, tempList);
                    CACs.put(config._lineages, tempList);
                }
                else{
                    if(node.getOutdeg() == 1){
                        Iterator<NetNode<Integer>> childNode = node.getChildren().iterator();
                        BitSet edge = new BitSet();
                        edge.set(node.getData());
                        edge.set(childNode.next().getData());
                        CACs.put(null,edge2ACminus.remove(edge));
                    }
                    else{
                        Iterator<NetNode<Integer>> childNode = node.getChildren().iterator();
                        BitSet edge1 = new BitSet();
                        edge1.set(node.getData());
                        edge1.set(childNode.next().getData());
                        List<Configuration> AC1 = edge2ACminus.remove(edge1);
                        BitSet edge2 = new BitSet();
                        edge2.set(node.getData());
                        edge2.set(childNode.next().getData());
                        List<Configuration> AC2 = edge2ACminus.remove(edge2);
                        edge2ACminus.remove(edge1);
                        edge2ACminus.remove(edge2);
                        boolean firstIndependent = _firstIndependentNodes.contains(node.getData());
                        boolean independent = _allIndependentNodes.contains(node.getData());

                        BitSet[][] netNodeLineages=null;
                        if(firstIndependent){
                            //netNodeLineages = new boolean[_netNodeNum][2][_gtClusters.size()];
                            netNodeLineages = new BitSet[_netNodeNum][2];
                            for(int i=0; i<_netNodeNum; i++){
                                for(int j=0; j<2; j++){
                                    netNodeLineages[i][j] = new BitSet();

                                }
                            }
                        }

                        if(_printDetail){
                            System.out.print("AC1: {");
                            for(Configuration config: AC1){
                                System.out.print(config.toString(gtClusters)+"  ");
                            }
                            System.out.println("}");
                            System.out.print("AC2: {");
                            for(Configuration config: AC2){
                                System.out.print(config.toString(gtClusters)+"  ");
                            }
                            System.out.println("}");
                        }

                        Configuration optimalOne = null;
                        int xl1=0,xl2=0;

                        for(Configuration config1: AC1){
                            for(Configuration config2: AC2){
                                if(config1.isCompatible(config2)){
                                    Configuration mergedConfig = new Configuration(config1, config2);
                                    if(firstIndependent){
                                        if(optimalOne==null || optimalOne._xl>mergedConfig._xl){
                                            optimalOne = mergedConfig;
                                            xl1 = Math.max(config1.getLineageCount() - 1, 0);
                                            xl2 = Math.max(config2.getLineageCount() - 1, 0);
                                            for(int i=0; i<_netNodeNum; i++){
                                                for(int j=0; j<2; j++){
                                                    netNodeLineages[i][j] = (BitSet)(mergedConfig._netNodeLineages[i][j].clone());
                                                    //netNodeLineages[i][j] = mergedConfig._netNodeLineages[i][j].clone();
                                                }
                                            }

                                        }

                                        else if(optimalOne._xl == mergedConfig._xl){
                                            for(int i=0; i<_netNodeNum; i++){
                                                for(int j=0; j<2; j++){
                                                    netNodeLineages[i][j].and(mergedConfig._netNodeLineages[i][j]);

                                                    //for(int k=0; k<_gtClusters.size(); k++)
                                                    //netNodeLineages[i][j][k] = netNodeLineages[i][j][k] && mergedConfig._netNodeLineages[i][j][k];

                                                    //System.out._printDetail(netNodeLineages[i][0].cardinality() + "/" + netNodeLineages[i][1].cardinality() + "   ");
                                                }
                                            }
                                        }

                                    }
                                    else if(independent){
                                        xl1 = Math.max(config1.getLineageCount() - 1, 0);
                                        xl2 = Math.max(config2.getLineageCount() - 1, 0);
                                        List<Configuration> sameLineageConfigs = new ArrayList<Configuration>();
                                        sameLineageConfigs.add(mergedConfig);
                                        CACs.put(mergedConfig._lineages, sameLineageConfigs);
                                    }

                                    else{
                                        List<Configuration> sameLineageConfigs = CACs.get(mergedConfig._lineages);
                                        if(sameLineageConfigs==null){
                                            sameLineageConfigs = new ArrayList<Configuration>();
                                            CACs.put(mergedConfig._lineages, sameLineageConfigs);
                                        }
                                        sameLineageConfigs.add(mergedConfig);
                                    }

                                }

                            }
                        }


                        if(firstIndependent){
                            optimalOne.setNetNodeLineageNum(netNodeLineages);
                            List<Configuration> tempList = new ArrayList<Configuration>();
                            tempList.add(optimalOne);
                            CACs.put(optimalOne._lineages, tempList);
                        }

                    }
                }



                if(_printDetail){
                    System.out.print("AC: {");
                    for(List<Configuration> configList: CACs.values()){
                        for(Configuration config: configList)
                            System.out.print(config.toString(gtClusters)+"  ");
                    }
                    System.out.println("}");
                }

                //set AC- for a node
                if(node.isRoot()){

                    if(CACs.size()!=1){
                        System.err.println("Error");
                    }
                    Configuration optimalConfig = CACs.values().iterator().next().get(0);
                    xl = optimalConfig._xl;

                    for(int i=0; i< _netNodeNum; i++){
                        //System.out._printDetail(optimalConfig._netNodeLineages[i][0].cardinality() + "/" + optimalConfig._netNodeLineages[i][1].cardinality() + "   ");
                        //System.out.print("\n"+optimalConfig._netNodeLineages[i][0].cardinality() + " " + optimalConfig._netNodeLineages[i][1].cardinality());
                        _totalNetNodeLinNum[i][0] += optimalConfig._netNodeLineages[i][0].cardinality() * weight;
                        _totalNetNodeLinNum[i][1] += optimalConfig._netNodeLineages[i][1].cardinality() * weight;
                    }


                }
                else if(node.isTreeNode()){
                    List<Configuration> ACminus = new ArrayList<Configuration>();
                    //for (Map<Set<Integer>, List<Configuration>> lineages2configs : CACs)
                    for(List<Configuration> sameLineageConfigs: CACs.values()){
                        Iterator<Configuration> configIt = sameLineageConfigs.iterator();
                        Configuration config = configIt.next();
                        //System.out.print(config);
                        Map<Integer,Integer> parent2child = new HashMap<Integer, Integer>();
                        //boolean canMerge;
                        List<Integer> lineageList = new ArrayList<Integer>();
                        lineageList.addAll(config._lineages);
                        do{
                            List<Integer> newLineageList = new ArrayList<Integer>();
                            for(Integer lineage: lineageList){
                                Integer parent = child2parent.get(lineage);
                                Integer sibling = parent2child.get(parent);
                                if(sibling == null){
                                    sibling = lineage;
                                    parent2child.put(parent,sibling);
                                }else{
                                    if(!node2outdegree.containsKey(parent)){
                                        config.mergeCluster(parent, sibling, lineage);
                                        newLineageList.add(parent);
                                    }
                                    else{
                                        int totalChildren = node2outdegree.get(parent);
                                        Integer resolvedDegree1 = addedNode2resolvedDegree.get(lineage);
                                        if(resolvedDegree1==null){
                                            resolvedDegree1 = 1;
                                        }
                                        Integer resolvedDegree2 = addedNode2resolvedDegree.get(sibling);
                                        if(resolvedDegree2==null){
                                            resolvedDegree2 = 1;
                                        }
                                        int resolvedDegree = resolvedDegree1+resolvedDegree2;
                                        if(resolvedDegree==totalChildren){
                                            config.mergeCluster(parent, sibling, lineage);
                                            newLineageList.add(parent);
                                        }
                                        else{
                                            int newVirtualNode = numGTNode++;
                                            //child2parent.put(sibling, newVirtualNode);
                                            //child2parent.put(lineage, newVirtualNode);
                                            child2parent.put(newVirtualNode, parent);
                                            addedNode2resolvedDegree.put(newVirtualNode, resolvedDegree);
                                            config.mergeCluster(newVirtualNode, sibling, lineage);
                                            newLineageList.add(newVirtualNode);
                                            parent2child.remove(parent);
                                            STITreeCluster newCluster = gtClusters.get(lineage).merge(gtClusters.get(sibling));
                                            gtClusters.add(newCluster);
                                        }
                                    }
                                }
                            }
                            if(config.getLineageCount()==1)break;
                            lineageList = newLineageList;

                        }while(lineageList.size()>0);
                        int addXL = Math.max(0, config.getLineageCount()-1);
                        config.addExtraLineage(addXL);
                        ACminus.add(config);
                        //System.out.println(" -> "+config);
                        while(configIt.hasNext()){
                            Configuration sameConfig = configIt.next();
                            sameConfig.setLineages(config._lineages);
                            sameConfig.addExtraLineage(addXL);
                            ACminus.add(sameConfig);
                        }
                    }

                    BitSet newEdge = new BitSet();
                    newEdge.set(node.getData());
                    newEdge.set(node.getParents().iterator().next().getData());
                    edge2ACminus.put(newEdge, ACminus);

                    if(_printDetail){
                        System.out.print("ACminus: {");
                        for(Configuration config: ACminus){
                            System.out.print(config.toString(gtClusters)+"  ");
                        }
                        System.out.println("}");
                    }
                }
                else {
                    List<Configuration> ACminus1 = new ArrayList<Configuration>();
                    List<Configuration> ACminus2 = new ArrayList<Configuration>();
                    int configIndex = 1;
                    //for(Map<Set<Integer>,List<Configuration>> lineages2configs: CACs)
                    for(List<Configuration> configList: CACs.values()){
                        for(Configuration config: configList){
                            int numLineage = config.getLineageCount();
                            int[] lineageArray = new int[numLineage];
                            int index = 0;
                            for(int lineage: config._lineages){
                                lineageArray[index++] = lineage;
                            }

                            boolean fEven = (config.getLineageCount()%2)==0;
                            int upper = config.getLineageCount()/2;
                            Set<Set<Integer>> addedConfigs = null;
                            if(fEven){
                                addedConfigs = new HashSet<Set<Integer>>();
                            }

                            for(int i=0; i<=upper; i++){
                                for(boolean[] selectedLineages: getSelected(config.getLineageCount(),i)){
                                    for(int k=0; k<2; k++){
                                        Configuration newConfig = new Configuration();
                                        newConfig.setNetNodeLineageNum(config._netNodeLineages);
                                        newConfig.setNetNodeChoice(config._netNodeIndex);

                                        index = 0;
                                        for(int lin: lineageArray) {
                                            if(selectedLineages[index] && k==0){
                                                newConfig.addLineage(lin);
                                                newConfig.addNetNodeLineageNum(netNodeIndex, 0, lin);
                                            }
                                            else if(!selectedLineages[index] && k==1){
                                                newConfig.addLineage(lin);
                                                newConfig.addNetNodeLineageNum(netNodeIndex, 0, lin);
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


                                        newConfig.addNetNodeChoice(netNodeIndex, configIndex);
                                        Configuration newConfigCopy = new Configuration(newConfig);
                                        if(config.getLineageCount()==0){
                                            newConfigCopy.addNetNodeChoice(netNodeIndex, configIndex);
                                        }else{
                                            newConfigCopy.addNetNodeChoice(netNodeIndex, configIndex+1);
                                        }


                                        if(k==0){
                                            newConfig.setExtraLineage(0);
                                            newConfigCopy.setExtraLineage(0);
                                            newConfigCopy.switchNetNodeLineageNum(netNodeIndex);
                                            ACminus1.add(newConfig);
                                            ACminus2.add(newConfigCopy);

                                        }
                                        else{
                                            newConfigCopy.setExtraLineage(config._xl);
                                            newConfig.setExtraLineage(config._xl);
                                            newConfig.switchNetNodeLineageNum(netNodeIndex);
                                            ACminus1.add(newConfigCopy);
                                            ACminus2.add(newConfig);

                                        }

                                    }
                                    if(config.getLineageCount()==0){
                                        configIndex++;
                                    }else{
                                        configIndex = configIndex + 2;
                                    }
                                }
                            }
                        }
                    }
                    if(_printDetail){
                        System.out.print("CAC after 1: {");
                        for(Configuration config: ACminus1){
                            System.out.print(config.toString(gtClusters)+"  ");
                        }
                        System.out.println("}");
                        System.out.print("CAC after 2: {");
                        for(Configuration config: ACminus2){
                            System.out.print(config.toString(gtClusters)+"  ");
                        }
                        System.out.println("}");
                    }
                    Iterator<NetNode<Integer>> it = node.getParents().iterator();
                    for(int i=0; i<2; i++){
                        List<Configuration> ACminus;
                        if(i==0){
                            ACminus = ACminus1;
                        }
                        else{
                            ACminus = ACminus2;
                        }
                        NetNode<Integer> parentNode = it.next();
                        if(node.getParentDistance(parentNode) != 0){
                            for(Configuration config: ACminus){
                                config.addExtraLineage(Math.max(0, config.getLineageCount()-1));
                            }
                        }
                        BitSet newEdge = new BitSet();
                        newEdge.set(node.getData());
                        newEdge.set(parentNode.getData());
                        edge2ACminus.put(newEdge, ACminus);

                        if(_printDetail){
                            System.out.print("ACminus to " + parentNode.getName()+ ": {");
                            for(Configuration config: ACminus){
                                System.out.print(config.toString(gtClusters)+"  ");
                            }
                            System.out.println("}");

                        }
                    }
                    netNodeIndex ++;
                }

            }

            xls[treeID] = xl;
            if(_parallel){
                treeID = getNextTreeID();
            }
            else{
                treeID++;
            }


        }

        if(!_parallel){
            double[] probs = getHybridProbabilities();
            int index = 0;
            for(NetNode<Integer> node: Networks.postTraversal(network)){
                if(node.isNetworkNode()){
                    for(NetNode parent: node.getParents()){
                        node.setParentProbability(parent,probs[index]);
                        probs[index] = 1 - probs[index];
                    }
                    index++;
                }
            }
        }


    }

    public void setPrint(boolean toPrint){
        _printDetail = toPrint;
    }

    private boolean checkLeafAgreement(Map<String, List<String>> species2alleles, String[] gtTaxa){
        if(species2alleles==null){
            for(String alleleG: gtTaxa){
                boolean found = false;
                for(String alleleS: _netTaxa){
                    if(alleleG.equals(alleleS)){
                        found = true;
                        break;
                    }
                }
                if(!found){
                    return false;
                }
            }
        }
        else{
            //HashSet<String> speciesSet = new HashSet<String>();
            for(String alleleG: gtTaxa){
                boolean found = false;
                for(Map.Entry<String, List<String>> entry: species2alleles.entrySet()){
                    for(String alleleS: entry.getValue()){
                        if(alleleS.equals(alleleG)){
                            found = true;
                            //speciesSet.add(entry.getKey());
                            break;
                        }
                    }
                    if(found){
                        break;
                    }
                }
                if(!found){
                    return false;
                }
            }

        }
        return true;
    }


    public void processNetwork(Network<Integer> net){
        removeBinaryNodes(net);
        _netNodeNum = 0;
        _totalNodeNum = 0;
        List<String> taxa = new ArrayList<String>();
        for(NetNode<Integer> node: net.dfs()){
            node.setData(_totalNodeNum++);
            if(node.isLeaf()){
                taxa.add(node.getName());
            }else if(node.isNetworkNode()){
                _netNodeNum++;
            }
        }
        _netTaxa = taxa.toArray(new String[0]);
        _totalNetNodeLinNum = new double[_netNodeNum][2];
        _firstIndependentNodes = new HashSet<Integer>();
        _allIndependentNodes = new HashSet<Integer>();
        computeNodeCoverage(net, _netTaxa, _firstIndependentNodes, _allIndependentNodes);
        _preProcessed = true;
    }


    private void processGT(Tree gt, String[] gtTaxa, Map<Integer,Integer> child2parent, Map<Integer, Integer> node2outdegree, List<STITreeCluster> gtClusters){
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
                int childCount = 0;
                for (TNode child : node.getChildren()) {
                    bs.or(map.get(child));
                    child2parent.put(((STINode<Integer>)child).getData(), index);
                    childCount++;
                }
                if(childCount>2){
                    node2outdegree.put(index, childCount);
                }
            }
            map.put(node, bs);
            STITreeCluster cl = new STITreeCluster(gtTaxa);
            cl.setCluster(bs);
            gtClusters.add(cl);
            index++;
        }
    }

    public static void removeBinaryNodes(Network<Integer> net)
    {
        // Find all binary nodes.
        List<NetNode> binaryNodes = new LinkedList<NetNode>();
        for (NetNode node : net.bfs()) {
            if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                binaryNodes.add(node);
            }
        }

        // Remove them.
        for (NetNode<Integer> node : binaryNodes) {
            NetNode child = node.getChildren().iterator().next();	// Node's only child.
            if(child.getIndeg() != 1){
                continue;
            }
            NetNode parent = node.getParents().iterator().next();	// Node's only parent.
            double distance = node.getParentDistance(parent) + child.getParentDistance(node);
            double gamma = node.getParentProbability(parent) * child.getParentProbability(node);
            parent.removeChild(node);
            node.removeChild(child);
            parent.adoptChild(child, distance);
            child.setParentProbability(parent, gamma);
        }
    }



    public static void computeNodeCoverage(Network<Integer> net, String[] netTaxa, Set<Integer> firstIndependentNodes, Set<Integer> allIndependentNodes){
        //List<Integer> leaves = new ArrayList<Integer>();

        //_node2Coverage = new HashMap<Integer, STITreeCluster>();
        //HashMap<NetNode<Integer>, STITreeCluster> node2Cluster = new HashMap<NetNode<Integer>, STITreeCluster>();
        for(NetNode<Integer> node: Networks.postTraversal(net)){
            int id = node.getData();
            STITreeCluster cl = new STITreeCluster(netTaxa);
            //System.out.println(cl);
            if(node.isLeaf()){
                allIndependentNodes.add(id);
                //System.out.println("all:" + node.getName());
                cl.addLeaf(node.getName());
            }
            else if(node.isRoot()){
                boolean ftotal = true;
                for(NetNode<Integer> child: node.getChildren()){
                    //cl = cl.merge(node2Cluster.get(child));
                    if(!allIndependentNodes.contains(child.getData())){
                        ftotal = false;
                    }
                }
                if(!ftotal){
                    firstIndependentNodes.add(id);
                    //System.out.println("first:" + node.getName());
                }
                allIndependentNodes.add(id);
                //System.out.println("all:" + node.getName());
            }
            else if(node.isTreeNode()){
                boolean ftotal = true;
                for(NetNode<Integer> child: node.getChildren()){
                    //cl = cl.merge(node2Cluster.get(child));
                    if(!allIndependentNodes.contains(child.getData())){
                        ftotal = false;
                    }
                }
                if(ftotal){
                    allIndependentNodes.add(id);
                    //System.out.println("all:" + node.getName());
                }else{
                    NetNode parent = node.getParents().iterator().next();
                    double distance = node.getParentDistance(parent);
                    parent.removeChild(node);
                    boolean disconnect = isValidNetwork(net, parent);
                    parent.adoptChild(node, distance);
                    if (disconnect) {
                        firstIndependentNodes.add(id);
                        //System.out.println("first:" + node.getName());
                        allIndependentNodes.add(id);
                        //System.out.println("all:" + node.getName());
                    }
                }

            }
            /*
            else{
                for(NetNode<Integer> child: node.getChildren()){
                    cl = cl.merge(node2Cluster.get(child));
                }
            }
            node2Cluster.put(node, cl);
            */

        }
        //System.out.println(firstIndependentNodes.size() + " " + allIndependentNodes.size());

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


    private static boolean isValidNetwork(Network<Integer> net, NetNode ignoreNode){
        Set<NetNode> visited = new HashSet<NetNode>();
        Set<NetNode> seen = new HashSet<NetNode>();
        for(NetNode<Integer> node: net.bfs()){
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



    private class Configuration{
        private HashSet<Integer> _lineages;
        private int _xl;
        int[] _netNodeIndex;
        BitSet[][] _netNodeLineages;

        public Configuration(){
            _lineages = new HashSet();
            _netNodeIndex = new int[_netNodeNum];

            _netNodeLineages = new BitSet[_netNodeNum][2];
            for(int i=0; i<_netNodeNum; i++){
                for(int j=0; j<2; j++){
                    _netNodeLineages[i][j] = new BitSet();

                }
            }

            Arrays.fill(_netNodeIndex, 0);

        }

        public Configuration(Configuration config){
            _lineages = (HashSet)config._lineages.clone();
            _xl = config._xl;
            _netNodeIndex = config._netNodeIndex.clone();


            _netNodeLineages = new BitSet[_netNodeNum][2];
            for(int i=0; i<_netNodeNum; i++){
                for(int j=0; j<2; j++){
                    _netNodeLineages[i][j] = (BitSet)(config._netNodeLineages[i][j].clone());

                }
            }

        }


        public Configuration(Configuration config1, Configuration config2){
            _lineages = (HashSet)config1._lineages.clone();
            _lineages.addAll(config2._lineages);
            _xl = config1._xl + config2._xl;
            _netNodeIndex = new int[_netNodeNum];

            _netNodeLineages = new BitSet[_netNodeNum][2];
            for(int i=0; i< _netNodeNum; i++){
                for(int j=0; j<2; j++){
                    _netNodeLineages[i][j] = new BitSet();
                }
            }


            for(int i=0; i<_netNodeNum; i++){
                if(config1._netNodeIndex[i] == config2._netNodeIndex[i]){
                    _netNodeIndex[i] = config1._netNodeIndex[i];
                }
                else{
                    _netNodeIndex[i] = Math.max(config1._netNodeIndex[i], config2._netNodeIndex[i]);
                }

                _netNodeLineages[i][0] = (BitSet)(config1._netNodeLineages[i][0].clone());
                _netNodeLineages[i][0].or(config2._netNodeLineages[i][0]);
                _netNodeLineages[i][1] = (BitSet)(config1._netNodeLineages[i][1].clone());
                _netNodeLineages[i][1].or(config2._netNodeLineages[i][1]);


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

        public void setLineages(Set<Integer> lineages){
            _lineages.clear();
            _lineages.addAll(lineages);
        }

        public void mergeCluster(int coalTo, int[] coalFrom){
            /*
            if(_lineages.get(coalFrom[0]) && _lineages.get(coalFrom[1])){
                _lineages.set(coalFrom[0], false);
                _lineages.set(coalFrom[1], false);
                _lineages.set(coalTo, true);
            }
            */
            _lineages.remove(coalFrom[0]);
            _lineages.remove(coalFrom[1]);
            _lineages.add(coalTo);

        }


        public void mergeCluster(int coalTo, int coalFrom1, int coalFrom2){
            /*
            if(_lineages.get(coalFrom[0]) && _lineages.get(coalFrom[1])){
                _lineages.set(coalFrom[0], false);
                _lineages.set(coalFrom[1], false);
                _lineages.set(coalTo, true);
            }
            */
            _lineages.remove(coalFrom1);
            _lineages.remove(coalFrom2);
            _lineages.add(coalTo);

        }




        public void setExtraLineage(int xl){
            _xl = xl;
        }

        public void addExtraLineage(int xl){
            _xl += xl;
        }

        public int getLineageCount(){
            return _lineages.size();
        }


        public String toString(List<STITreeCluster> gtClusters){
            String exp = "";
            for(int id: _lineages) {
                exp = exp + gtClusters.get(id);
            }
            exp = exp + "/[";

            /*
            for(int i=0; i<_netNodeLineages.length; i++){
                exp = exp + _netNodeLineages[i][0].cardinality() + "/" + _netNodeLineages[i][1].cardinality();
                if(i!=_netNodeLineages.length-1){
                    exp = exp + ",";
                }
            }
            */
            exp = exp + "]:" + _xl;
            return exp;
        }

        public void addNetNodeChoice(int net, int index){
            _netNodeIndex[net] = index;
        }

        public void setNetNodeChoice(int[] choice){
            _netNodeIndex = choice.clone();
        }


        public void addNetNodeLineageNum(int net, int index, int lineage){
            _netNodeLineages[net][index].set(lineage);
        }

        public void switchNetNodeLineageNum(int net){
            BitSet temp = _netNodeLineages[net][0];
            _netNodeLineages[net][0] = _netNodeLineages[net][1];
            _netNodeLineages[net][1] = temp;

        }


        public void setNetNodeLineageNum(BitSet[][] lineageNum){

            for(int i=0; i<_netNodeNum; i++){
                for(int j=0; j<2; j++){
                    _netNodeLineages[i][j] = (BitSet)(lineageNum[i][j].clone());
                    //_netNodeLineages[i][j] = lineageNum[i][j].clone();
                }
            }
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
