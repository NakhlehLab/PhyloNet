package edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.genetreeprobability;

import edu.rice.cs.bioinfo.library.phylogenetics.*;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.STITreeCluster;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.rice.cs.bioinfo.library.programming.Predicate2;
import edu.rice.cs.bioinfo.library.programming.counters.Counter;
import edu.rice.cs.bioinfo.library.programming.counters.CounterInt;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/26/12
 * Time: 4:02 PM
 * To change this template use File | Settings | File Templates.
 */

public class GeneTreeProbabilityYF<NN,NE,TN,TE>
{
  //  Set<Integer> _totalCoverNodes;
 //   String[] _netTaxa;
    boolean[][] _R;
    boolean _printDetails = false;
 //   int _netNodeNum;
 //   int _totalNodeNum;
   // List<STITreeCluster> _gtClusters;

    private Predicate2<NN, GraphReadOnly<NN,NE>> _isNetworkNodeStrategy = new IsNetworkNode<NN, NE>();

    private Func2<GraphReadOnly<NN,NE>, Counter<Integer>, Integer> _getNetworkNodeCountStrategy = new GetNetworkNodeCount<NN, NE, Integer>();

    private Func1<GraphReadOnly<TN,?>, Iterable<TN>> _getLeafsStrategy = new GetLeafs<TN>();

    private Func1<GraphReadOnly<TN,TE>, Iterable<TN>> _getGeneTreeNodesPostOrderStrategy = new GetNodesPostOrder<TN, TE>();

    private Func2<GraphReadOnly<TN,TE>,TN,Iterable<TN>> _getGeneTreeDirectSuccessorsStrategy = new GetDirectSuccessors<TN, TE>();

    private Func1<GraphReadOnly<NN,?>, NN> _getNetworkRootStrategy = new FindRoot<NN>();

    private Func2<GraphReadOnly<NN,NE>,NN,Integer> _getNetworkNodeOutDegreeStrategy = new GetOutDegree<NN, NE>();

    private Func2<GraphReadOnly<NN,NE>,NN,Iterable<NN>> _getNetworkNodeDirectSucessorsStrategy = new GetDirectSuccessors<NN, NE>();

     private Func2<GraphReadOnly<NN,NE>,NN,Iterable<NN>> _getNetworkNodeDirectPredecessorsStrategy = new GetDirectPredecessors<NN, NE>();

    private Predicate2<GraphReadOnly, Object> _isNetworkLeafStrategy = new IsLeaf();

    private Func2<NN, GraphReadOnly<NN,NE>,Boolean> _isNetworkTreeNodeStrategy = new IsTreeNode<NN, NE>();

    private Func1<GraphReadOnly<NN,NE>, Iterable<NN>> _getNetworkNodesDFSOrderStrategy = new GetNodesDFSOrder<NN,NE>();

    private Predicate2<GraphReadOnly, Object> _isGeneTreeLeafStrategy = new IsLeaf();

    private Predicate1<GraphReadOnly> _isNetworkWeeklyConnectedStrategy = new IsConnected();


    public void setPrintDetails(boolean p){
        _printDetails = p;
    }


    public List<Double> calculateGTDistribution(GraphReadOnly<NN,NE> network, List<GraphReadOnly<TN,TE>> gts,
                                                Func2<TN, GraphReadOnly<TN,TE>, String> getGeneTreeNodeName,
                                                Func2<NN, GraphReadOnly<NN,NE>, String> getNetworkNodeName,
                                                Func2<NE, GraphReadOnly<NN,NE>, Double> getDistance,
                                                Func2<NE, GraphReadOnly<NN,NE>, Double> getProbability,
                                                Map<String, List<String>> species2alleles){

        new GraphValidator<NN,NE>().assertValidGraph(network);

        for(NE edge : network.getEdges())
        {
            if(getDistance.execute(edge, network) < 0)
            {
                throw new IllegalArgumentException("Given network contains an edge with distance less than zero. (" + edge + ")");
            }
        }


        NN netRoot = _getNetworkRootStrategy.execute(network);
        int totalNodeNum = IterableHelp.countInt(network.getNodes());
        int netNodeNum = _getNetworkNodeCountStrategy.execute(network, new CounterInt());

        List<Double> probList = new ArrayList<Double>();
        Map<NN,Integer> nodeToData = processNetwork(network);
        Set<Integer> totalCoverNodes = computeNodeCoverage(network, nodeToData, netRoot, getDistance);
        //System.out.println(gts.size());
        //System.exit(0);
        int tempCount = 0;
        List<STITreeCluster> gtClusters = new LinkedList<STITreeCluster>();
        for(GraphReadOnly<TN,TE> gt: gts){
            double gtProb = 0;

            List<String> gtTaxaNames = new LinkedList<String>();

            for(TN gtLeaf : _getLeafsStrategy.execute(gt))
            {
                gtTaxaNames.add(getGeneTreeNodeName.execute(gtLeaf,gt));
            }

            String[] gtTaxa = (String[]) gtTaxaNames.toArray(new String[0]);

            Map<Integer,Integer> child2parent = new HashMap<Integer, Integer>();
            gtClusters.addAll(processGT(gt, gtTaxa, child2parent, getGeneTreeNodeName));
            computeR(gtClusters);

            HashSet<String> gtTaxaSet = new HashSet<String>();
            for(String taxon: gtTaxa){
                gtTaxaSet.add(taxon);
            }

            HashMap<BitSet, List<Configuration>> edge2ACminus = new HashMap<BitSet, List<Configuration>>();
            int netNodeIndex = 0;
            //int nodeID = 0;
            for(NN node: walkNetwork(network)){

                //TODO
                //System.out.println();
                //System.out.println("On node #" + node.getData() + " " + node.getName());
            /*    if(_printDetails){
                    System.out.println();
                    System.out.println("On node #" + node.getData() + " " + node.getName());
                } */
                List<Map<Set<Integer>,List<Configuration>>> CACs = new ArrayList<Map<Set<Integer>,List<Configuration>>>();

                //long start = System.currentTimeMillis();
                //set AC for a node
                if(_isNetworkLeafStrategy.execute(network,  node)){
                    Map<Set<Integer>,List<Configuration>> sizeOneConfigs = new HashMap<Set<Integer>, List<Configuration>>();
                    Configuration config = new Configuration(netNodeNum);
                    String nodeName = getNetworkNodeName.execute(node, network);
                    if(species2alleles == null){
                        if(gtTaxaSet.contains(nodeName)){
                            STITreeCluster cl = new STITreeCluster(gtTaxa);
                            cl.addLeaf(nodeName);
                            config.addLineage(gtClusters.indexOf(cl));
                        }
                    }
                    else{
                        for(String allele: species2alleles.get(nodeName)){
                            if(gtTaxaSet.contains(allele)){
                                STITreeCluster cl = new STITreeCluster(gtTaxa);
                                cl.addLeaf(allele);
                                config.addLineage(gtClusters.indexOf(cl));
                            }
                        }
                    }
                    config.setProbability(1);
                    List<Configuration> tempList = new ArrayList<Configuration>();
                    tempList.add(config);
                    sizeOneConfigs.put(config._lineages, tempList);
                    CACs.add(sizeOneConfigs);
                }
                else{
                    if(_getNetworkNodeOutDegreeStrategy.execute(network, node) == 1){
                        Iterator<NN> childNode = _getNetworkNodeDirectSucessorsStrategy.execute(network,  node).iterator();
                        BitSet edge = new BitSet();
                        int nodeData = nodeToData.get(node);
                        edge.set(nodeData);
                        int childNodeNextData = nodeToData.get(childNode.next());
                        edge.set(childNodeNextData);
                        Map<Set<Integer>,List<Configuration>> temp = new HashMap<Set<Integer>, List<Configuration>>();
                        temp.put(null,edge2ACminus.remove(edge));
                        CACs.add(temp);
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
                        Iterator<NN> childNode = _getNetworkNodeDirectSucessorsStrategy.execute(network, node).iterator();
                        BitSet edge1 = new BitSet();
                        int nodeData = nodeToData.get(node);
                        edge1.set(nodeData);
                        edge1.set(nodeToData.get(childNode.next()));
                        List<Configuration> AC1 = edge2ACminus.remove(edge1);
                        BitSet edge2 = new BitSet();
                        edge2.set(nodeData);
                        edge2.set(nodeToData.get(childNode.next()));
                        List<Configuration> AC2 = edge2ACminus.remove(edge2);
                        //TODO
                        //System.out.println("AC1: " + AC1.size() + " & AC2:" + AC2.size());
                        //System.out.println("total");

                        //TODO see if needed to add more
                        List<Integer> configSizeList = new ArrayList<Integer>();
                        boolean total = totalCoverNodes.contains(nodeData);

                        for(Configuration config1: AC1){
                            for(Configuration config2: AC2){
                                if(config1.isCompatible(config2)){
                                    Configuration mergedConfig = new Configuration(config1, config2, netNodeNum);
                                    if(mergedConfig._prob <=0 ){
                                        continue;
                                    }
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
                                        }
                                        else{
                                            sameLineageConfigs.get(0).addProbability(mergedConfig._prob);
                                        }

                                    }
                                    else{
                                        sameLineageConfigs.add(mergedConfig);

                                    }
                                }
                            }
                        }
                        //System.out.println("Pooling time: " + (System.currentTimeMillis()-start)/1000.0);
                    }
                }

                /*
                int total = 0;
                for(Map<Set<Integer>,List<Configuration>> lineages2configs: CACs){
                    for(List<Configuration> configList: lineages2configs.values()){
                        total += configList.size();
                    }
                }
                System.out.println("CAC size: " + total);
                */

                if(_printDetails){
                    System.out.print("AC: {");
                    for(Map<Set<Integer>,List<Configuration>> lineages2configs: CACs){
                        for(List<Configuration> configList: lineages2configs.values()){
                            for(Configuration config: configList){
                                System.out.print(config.toString()+"  ");
                            }
                        }
                    }
                    System.out.println("}");
                }


                //set AC- for a node
                if(node == netRoot){
                    Configuration rootConfig = new Configuration(netNodeNum);
                    rootConfig.addLineage(gtClusters.size()-1);
                    for(Map<Set<Integer>,List<Configuration>> lineages2configs: CACs){
                        for(List<Configuration> configList: lineages2configs.values()){
                            for(Configuration preConfig: configList){
                                if(preConfig.getLineageCount()==1){
                                    gtProb += preConfig._prob;
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
                                    rootConfig.clearCoalEvents();
                                    rootConfig.addCoalEvents(events);
                                    double prob = computeProbability(preConfig, rootConfig, -1, 1);
                                    if(prob!=0){
                                        gtProb += prob*preConfig._prob;
                                    }
                                }

                            }
                        }

                    }
                    probList.add(gtProb);
                }
                else if(_isNetworkTreeNodeStrategy.execute(node, network) || _isNetworkLeafStrategy.execute(network, node) ){
                    NN parent = _getNetworkNodeDirectPredecessorsStrategy.execute(network, node).iterator().next();
                    NE parentEdge = network.getEdge(parent, node);
                    double distance = getDistance.execute(parentEdge, network);
                    List<Configuration> ACminus = new ArrayList<Configuration>();

                    if(distance==0){
                        for(Map<Set<Integer>,List<Configuration>> lineages2configs: CACs){
                            for(List configs: lineages2configs.values()){
                                ACminus.addAll(configs);
                            }
                        }
                    }
                    else{
                        computeACMinus(CACs, distance, 1, child2parent, ACminus);
                    }
                    BitSet newEdge = new BitSet();
                    newEdge.set(nodeToData.get(node));
                    newEdge.set(nodeToData.get(parent));
                    edge2ACminus.put(newEdge, ACminus);
                    //TODO
                    //System.out.println("ACminus: " + ACminus.size());
                    if(_printDetails){
                        System.out.print("ACminus: {");
                        for(Configuration config: ACminus){
                            System.out.print(config.toString()+"  ");
                        }
                        System.out.println("}");
                    }
                }
                else {
                    List<Integer> configSizeList = new ArrayList<Integer>();
                    List<Map<Set<Integer>,List<Configuration>>> newCACs1 = new ArrayList<Map<Set<Integer>,List<Configuration>>>();
                    List<Map<Set<Integer>,List<Configuration>>> newCACs2 = new ArrayList<Map<Set<Integer>,List<Configuration>>>();

                    int configIndex = 1;
                    for(Map<Set<Integer>,List<Configuration>> lineages2configs: CACs){
                        for(List<Configuration> configList: lineages2configs.values()){
                            for(Configuration config: configList){
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
                                            Configuration newConfig = new Configuration(netNodeNum);
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
                                            newConfig.setProbability(Math.sqrt(config._prob));

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
                                    System.out.print(config.toString()+"  ");
                                }
                            }
                        }
                        System.out.println("}");
                        System.out.print("CAC after: {");
                        for(Map<Set<Integer>,List<Configuration>> lineages2configs: newCACs2){
                            for(List<Configuration> configList: lineages2configs.values()){
                                for(Configuration config: configList){
                                    System.out.print(config.toString()+"  ");
                                }
                            }
                        }
                        System.out.println("}");
                    }
                    Iterator<NN> it = _getNetworkNodeDirectPredecessorsStrategy.execute(network, node).iterator();
                    NN parentNode1 = it.next();
                    NE parentNode1Edge = network.getEdge(parentNode1, node);
                    double distance1 = getDistance.execute(parentNode1Edge, network);
                    double hybridProb1 = getProbability.execute(parentNode1Edge, network);
                    NN parentNode2 = it.next();
                    NE parentNode2Edge = network.getEdge(parentNode2, node);
                    double distance2 = getDistance.execute(parentNode2Edge, network);
                    double hybridProb2 = getProbability.execute(parentNode2Edge, network);
                    List<Configuration> ACminus1 = new ArrayList<Configuration>();
                    List<Configuration> ACminus2 = new ArrayList<Configuration>();
                    computeTwoACMinus(newCACs1, distance1, hybridProb1, newCACs2, distance2, hybridProb2,child2parent, ACminus1, ACminus2);
                    BitSet newEdge1 = new BitSet();
                    newEdge1.set(nodeToData.get(node));
                    newEdge1.set(nodeToData.get(parentNode1));
                    edge2ACminus.put(newEdge1, ACminus1);
                  /*  if(_printDetails){
                        System.out.print("ACminus to " + parentNode1.getName()+ ":  {");
                        for(Configuration config: ACminus1){
                            System.out.print(config.toString()+"  ");
                        }
                        System.out.println("}");
                    } */
                    BitSet newEdge2 = new BitSet();
                    newEdge2.set(nodeToData.get(node));
                    newEdge2.set(nodeToData.get(parentNode2));
                    edge2ACminus.put(newEdge2, ACminus2);
                 /*   if(_printDetails){
                        System.out.print("ACminus to " + parentNode2.getName()+ ":  {");
                        for(Configuration config: ACminus2){
                            System.out.print(config.toString()+"  ");
                        }
                        System.out.println("}");
                    }   */
                    netNodeIndex ++;

                }


            }
            if(_printDetails){
                System.out.println("The probability of this gene tree is:" + gtProb);
            }
        }
        return probList;
    }


    private void computeACMinus(List<Map<Set<Integer>,List<Configuration>>> CACs, double distance, double hybridProb, Map<Integer, Integer> child2parent, List<Configuration> ACminus){
        Map<Integer, Map<Configuration, Configuration>> shape2ACminus = new HashMap<Integer, Map<Configuration, Configuration>>();
        for (Map<Set<Integer>, List<Configuration>> lineages2configs : CACs) {
            for (List<Configuration> sameLineageConfigs : lineages2configs.values()) {

                Configuration origConfig = sameLineageConfigs.get(0);
                if (origConfig._prob <= 0) {
                    throw new RuntimeException("Wrong probability!");
                }
                Stack<Configuration> configStack = new Stack<Configuration>();
                Configuration configCopy = new Configuration(origConfig);
                configCopy.clearCoalEvents();
                configStack.push(configCopy);
                Map<Integer, Set<Set<Integer>>> visitedACminus = new HashMap<Integer, Set<Set<Integer>>>();

                while (!configStack.empty()) {
                    Configuration cconfig = configStack.pop();
                    boolean allNegative = true;
                    double prob = computeProbability(origConfig, cconfig, distance, hybridProb);
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
                        double newProb = prob*config._prob;
                        if(newProb<=0){
                            continue;
                        }
                        allNegative = false;
                        Configuration cconfigCopy;
                        if(ffirst){
                            cconfigCopy = cconfig;
                            ffirst = false;
                        }
                        else{
                            cconfigCopy = new Configuration(cconfig);
                            cconfigCopy.setNetNodeChoice(config._netNodeIndex);
                        }
                        cconfigCopy.setProbability(newProb);
                        Map<Configuration, Configuration> cc = shape2ACminus.get(code);
                        if (cc != null) {
                            if (cc.containsKey(cconfigCopy)) {
                                cc.get(cconfigCopy).addProbability(cconfigCopy._prob);
                            } else {
                                cc.put(cconfigCopy, cconfigCopy);
                            }
                        } else {
                            cc = new HashMap<Configuration, Configuration>();
                            cc.put(cconfigCopy, cconfigCopy);
                            shape2ACminus.put(code, cc);
                        }
                    }

                    if(allNegative){
                        continue;
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


    private void computeTwoACMinus(List<Map<Set<Integer>,List<Configuration>>>CACs1, double distance1, double hybridProb1, List<Map<Set<Integer>,List<Configuration>>>CACs2, double distance2, double hybridProb2,Map<Integer, Integer> child2parent, List<Configuration> ACminus1, List<Configuration> ACminus2){
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
                Configuration origConfig2 = sameLineageConfigs2.get(0);
                if(origConfig1._prob<=0){
                    throw new RuntimeException("Wrong probability!");
                }
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
                    boolean allNegative = true;
                    //double prob = computeProbability(origConfig1, cconfig, distance1, hybridProb1);
                    String[] forPrint = new String[2];
                    double probCommon = computeProbabilityPart1(origConfig1, cconfig, forPrint);
                    double prob = computeProbabilityPart2(probCommon, origConfig1.getLineageCount(), cconfig.getLineageCount(), distance1, hybridProb1, forPrint);

                    boolean ffirst = true;
                    for(Configuration config: sameLineageConfigs1){
                        double newProb = prob*config._prob;
                        if(newProb<=0){
                            continue;
                        }
                        allNegative = false;
                        Configuration cconfigCopy;
                        if(ffirst){
                            cconfigCopy = cconfig;
                            ffirst = false;
                        }
                        else{
                            cconfigCopy = new Configuration(cconfig);
                            cconfigCopy.setNetNodeChoice(config._netNodeIndex);
                        }
                        cconfigCopy.setProbability(newProb);
                        Map<Configuration, Configuration> cc = shape2ACminus1.get(code);
                        if (cc != null) {
                            if (cc.containsKey(cconfigCopy)) {
                                cc.get(cconfigCopy).addProbability(cconfigCopy._prob);
                            } else {
                                cc.put(cconfigCopy, cconfigCopy);
                            }
                        } else {
                            cc = new HashMap<Configuration, Configuration>();
                            cc.put(cconfigCopy, cconfigCopy);
                            shape2ACminus1.put(code, cc);
                        }
                    }

                    prob = computeProbabilityPart2(probCommon, origConfig1.getLineageCount(), cconfig.getLineageCount(), distance2, hybridProb2, forPrint);

                    //prob = computeProbability(origConfig2, cconfig, distance2, hybridProb2);
                    for(Configuration config: sameLineageConfigs2){
                        double newProb = prob*config._prob;
                        if(newProb<=0){
                            continue;
                        }
                        allNegative = false;
                        Configuration cconfigCopy = new Configuration(cconfig);
                        cconfigCopy.setNetNodeChoice(config._netNodeIndex);
                        cconfigCopy.setProbability(newProb);
                        Map<Configuration, Configuration> cc = shape2ACminus2.get(code);
                        if (cc != null) {
                            if (cc.containsKey(cconfigCopy)) {
                                cc.get(cconfigCopy).addProbability(cconfigCopy._prob);
                            } else {
                                cc.put(cconfigCopy, cconfigCopy);
                            }
                        } else {
                            cc = new HashMap<Configuration, Configuration>();
                            cc.put(cconfigCopy, cconfigCopy);
                            shape2ACminus2.put(code, cc);
                        }
                    }

                    if(allNegative){
                        continue;
                    }


                    Map<Integer, List<Integer>> parent2children = new HashMap<Integer, List<Integer>>();
                    for(int i: cconfig._lineages)
                    {
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

    private Map<NN,Integer> processNetwork(GraphReadOnly<NN,NE> net){
        Map<NN,Integer> nodeToData = new HashMap<NN, Integer>();
    //    removeBinaryNodes(net);
      //  _netNodeNum = 0;
      int  totalNodeNum = 0;
    //    List<String> taxa = new ArrayList<String>();
        for(NN node: _getNetworkNodesDFSOrderStrategy.execute(net)){
            nodeToData.put(node, totalNodeNum++);
       //     if(node.isLeaf()){
      //          taxa.add(node.getName());
           /* }else */ if(_isNetworkNodeStrategy.execute(node, net)){
        //        _netNodeNum++;
            }
        }
    //    _netTaxa = taxa.toArray(new String[0]);

        nodeToData.clear();
        nodeToData.put((NN)"R", 0);
        nodeToData.put((NN)"K", 1);
        nodeToData.put((NN)"X", 2);
        nodeToData.put((NN)"L", 3);
        nodeToData.put((NN)"C", 4);
        nodeToData.put((NN)"B", 5);
        nodeToData.put((NN)"D", 6);
        nodeToData.put((NN)"J", 7);
        nodeToData.put((NN)"A", 8);

        return nodeToData;
    }


    private List<STITreeCluster> processGT(GraphReadOnly<TN,TE> gt, String[] gtTaxa, Map<Integer,Integer> child2parent, Func2<TN, GraphReadOnly<TN,TE>, String> getGeneTreeNodeName){

        List<STITreeCluster> gtClusters = new ArrayList<STITreeCluster>();
        Map<TN,Integer> nodeToData = new HashMap<TN, Integer>();

        Map<TN, BitSet> map = new HashMap<TN, BitSet>();
        int index = 0;
        for (TN node : _getGeneTreeNodesPostOrderStrategy.execute(gt)) {
            //((STINode<Integer>)node).setData(index);
            nodeToData.put(node, index);
            BitSet bs = new BitSet();

            if (_isGeneTreeLeafStrategy.execute(gt, node)) {

                for (int i = 0; i < gtTaxa.length; i++) {
                    String nodeName = getGeneTreeNodeName.execute(node, gt);
                    if (nodeName.equals(gtTaxa[i])) {
                        bs.set(i);
                        break;
                    }
                }
            }
            else {
                for (TN child : _getGeneTreeDirectSuccessorsStrategy.execute(gt, node)) {
                    bs.or(map.get(child));
                    child2parent.put(nodeToData.get(child), index);
                }
            }
            map.put(node, bs);
            STITreeCluster cl = new STITreeCluster(gtTaxa);
            cl.setCluster(bs);
            gtClusters.add(cl);
            index++;
        }

        return gtClusters;
    }
          /*
    private void removeBinaryNodes(Network<Integer> net)
    {
        // Find all binary nodes.
        List<NetNode<Integer>> binary_nodes = new LinkedList<NetNode<Integer>>();
        for (NetNode node : net.bfs()) {
            if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                binary_nodes.add(node);
            }
        }

        // Remove them.
        for (NetNode<Integer> node : binary_nodes) {
            NetNode<Integer> child = node.getChildren().iterator().next();	// Node's only child.
            if(child.getIndeg() != 1){
                continue;
            }
            NetNode<Integer> parent = node.getParents().iterator().next();	// Node's only parent.
            double distance = node.getParentDistance(parent) + child.getParentDistance(node);
            double gamma = node.getParentProbability(parent) * child.getParentProbability(node);
            parent.removeChild(node);
            node.removeChild(child);
            parent.adoptChild(child, distance);
            child.setParentProbability(parent, gamma);
        }
    } */

    private List<NN> walkNetwork(GraphReadOnly<NN,NE> net){
        Stack<NN> stack = new Stack<NN>();
        List<NN> searchedNodes = new ArrayList<NN>();
        NN netRoot = _getNetworkRootStrategy.execute(net);
        stack.push(netRoot);
        Map<NN, Integer> node2index = new HashMap<NN, Integer>();
        node2index.put(netRoot, 0);

        while(!stack.isEmpty()){
            NN topNode = stack.peek();
            int index = node2index.get(topNode);
            if(index == _getNetworkNodeOutDegreeStrategy.execute(net, topNode)){
                searchedNodes.add(stack.pop());
            }
            else{
                Iterator<NN> it = _getNetworkNodeDirectSucessorsStrategy.execute(net, topNode).iterator();
                for(int i=0; i<index; i++){
                    it.next();
                }
                NN child = it.next();
                if(searchedNodes.contains(child)){
                    node2index.put(topNode, index + 1);
                }
                else{
                    stack.push(child);
                    node2index.put(child, 0);
                }
            }
        }

        return (List<NN>) Arrays.asList("A", "B", "C", "L", "X", "J", "D", "K", "R");

    //    return searchedNodes;
    }


    private  HashSet<Integer> computeNodeCoverage(GraphReadOnly<NN,NE> net, Map<NN, Integer> nodeToData, NN netRoot, Func2<NE, GraphReadOnly<NN,NE>, Double> getDistance){
        //List<Integer> leaves = new ArrayList<Integer>();
        List<Integer> allTotalNodes = new ArrayList<Integer>();
        HashSet<Integer> totalCoverNodes = new HashSet<Integer>();
        for(NN node: walkNetwork(net)){
            int id = nodeToData.get(node);
            if(_isNetworkLeafStrategy.execute(net, node)){
                //leaves.add(id);
                allTotalNodes.add(id);
            }

            else if(netRoot == node){
                boolean ftotal = true;
                for(NN child: _getNetworkNodeDirectSucessorsStrategy.execute(net, node)){

                    if(!allTotalNodes.contains(nodeToData.get(child))){
                        ftotal = false;
                        break;
                    }
                }
                if(!ftotal){
                    totalCoverNodes.add(id);
                }

            }
            else if(_isNetworkTreeNodeStrategy.execute(node, net)){
                boolean ftotal = true;
                for(NN child: _getNetworkNodeDirectSucessorsStrategy.execute(net, node)){
                    if(!allTotalNodes.contains(nodeToData.get(child))){
                        ftotal = false;
                        break;
                    }
                }
                if(ftotal){
                    allTotalNodes.add(id);
                }else{
                    NN parent = _getNetworkNodeDirectPredecessorsStrategy.execute(net, node).iterator().next();
                    NE parentToNode = net.getEdge(parent, node);
                    double distance = getDistance.execute(parentToNode, net);
                    GraphMask<NN,NE,GraphReadOnly<NN,NE>> mask = new GraphMask<NN, NE, GraphReadOnly<NN, NE>>(net, new ArrayList<NN>(), Arrays.asList(parentToNode));
                  //  parent.removeChild(node);
                    boolean disconnect = !_isNetworkWeeklyConnectedStrategy.execute(mask); // isValidNetwork(net);
               //     parent.adoptChild(node, distance);
                    if (disconnect) {
                        totalCoverNodes.add(id);
                        allTotalNodes.add(id);
                    }
                }

            }
        }

        return totalCoverNodes;
    }
     /*
    private boolean isValidNetwork(Network<Integer> net){
        BitSet visited = new BitSet();
        BitSet seen = new BitSet();
        for(NetNode<Integer> node: net.bfs()){
            visited.set(node.getData(), true);
            for(NetNode<Integer> parent: node.getParents()){
                seen.set(parent.getData(), true);
            }
            for(NetNode<Integer> child: node.getChildren()){
                seen.set(child.getData(), true);
            }
        }
        return visited.cardinality()==seen.cardinality();
    } */




    /**
     * The function is to calculate the _R matrix for the given tree.
     * Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
     */
    private void computeR(List<STITreeCluster> gtClusters){
        _R = new boolean[gtClusters.size()][gtClusters.size()];
        for(int i=0; i<gtClusters.size(); i++){
            STITreeCluster cl1 = gtClusters.get(i);
            for(int j=i+1; j<gtClusters.size(); j++){
                STITreeCluster cl2 = gtClusters.get(j);
                if(cl1.containsCluster(cl2)){
                    _R[i][j] = true;
                }
                else if(cl2.containsCluster(cl1)){
                    _R[j][i] = true;
                }
            }
        }
    }



    private double computeProbability(Configuration preConfig, Configuration coalescedConfig, double distance, double portion){
        double prob;
        int u = preConfig.getLineageCount();
        int v = coalescedConfig.getLineageCount();
        prob = Math.pow(portion, u);
        if(distance == 0 && u != v){
            if(_printDetails){
                System.out.println(preConfig.toString()+"->"+coalescedConfig.toString()+ ": g"+u+v+"(0)=0");
            }
            return 0;
        }
        if(u == v && (u == 1 || u == 0)){
            if(_printDetails){
                if(portion==1 || u==0){
                    System.out.println(preConfig.toString()+"->"+coalescedConfig.toString()+ ": g"+u+v+"("+distance+")=1");
                }
                else{
                    System.out.println(preConfig.toString()+"->"+coalescedConfig.toString()+ ": g"+u+v+"("+distance+")*"+portion+"="+portion);
                }
            }
            return prob;
        }
        double gij = gij(distance, u, v);
        if(gij<=0){
            return -1;
        }
        double d = calculateD(u,u-v);
        double w = calculateW(coalescedConfig._coalEvents);

        if(_printDetails){
            if(portion!=1){
                System.out.println(preConfig.toString()+"->"+coalescedConfig.toString()+ ": g"+u+v+"("+distance+")*"+w+"/"+d+"*"+portion+"^"+u+"=" + gij*w/d*Math.pow(portion, u));
            }
            else{
                System.out.println(preConfig.toString()+"->"+coalescedConfig.toString()+ ": g"+u+v+"("+distance+")*"+w+"/"+d+"=" + gij*w/d);
            }
        }
        prob *= gij*w/d;
        return prob;
    }


    private double computeProbabilityPart1(Configuration preConfig, Configuration coalescedConfig, String[] forPrint){
        double prob;
        forPrint[0] = preConfig.toString()+"->"+coalescedConfig.toString()+ ": ";
        int u = preConfig.getLineageCount();
        int v = coalescedConfig.getLineageCount();
        if(u == v && (u == 1 || u == 0)){
            return 1;
        }
        double d = calculateD(u,u-v);
        double w = calculateW(coalescedConfig._coalEvents);
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
      /*  if(length == TNode.NO_DISTANCE || length == -1){
            if(j == 1){
                return 1;
            }else{
                return 0;
            }
        } */
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


    private double calculateW(Set<Integer> coalEvents){
        double w = 1.0;
        w = w*fact(1, coalEvents.size());
        //System.out.print(fact(1, coalEvents.cardinality()));
        for (int i: coalEvents) {
            int a = 0;
            for (int j: coalEvents) {
                if(i!=j && _R[i][j])
                    a++;
            }
            //System.out.print(" * " + 1.0/(1 + a));
            w = w*(1.0/(1 + a));
        }
        //System.out.println();
        return w;
    }


    private double calculateD(Set<Integer> coalEvents){
        double w = 1.0;
        w = w*fact(1, coalEvents.size());
        //System.out.print(fact(1, coalEvents.cardinality()));
        for (int i: coalEvents) {
            int a = 0;
            for (int j: coalEvents) {
                if(i!=j && _R[i][j])
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




    private static class Configuration{
        private HashSet<Integer> _lineages;
        private double _prob;
        int[] _netNodeIndex;
        private Set<Integer> _coalEvents;

        private int _netNodeCount;

        public Configuration(int netNodeCount){
            _netNodeCount = netNodeCount;
            _lineages = new HashSet<Integer>();
            _netNodeIndex = new int[_netNodeCount];
            Arrays.fill(_netNodeIndex, 0);
        }

        public Configuration(Configuration config){
            _lineages = (HashSet)config._lineages.clone();
            _prob = config._prob;
            _netNodeIndex = config._netNodeIndex.clone();
        }

        public Configuration(Configuration config1, Configuration config2, int netNodeCount){
            _netNodeCount = netNodeCount;
            _lineages = (HashSet)config1._lineages.clone();
            _lineages.addAll(config2._lineages);
            _prob = config1._prob*config2._prob;
            _netNodeIndex = new int[_netNodeCount];
            for(int i=0; i< _netNodeCount; i++){
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
            //TODO needed?
            for(int i=0; i< _netNodeCount; i++){
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

        public void setProbability(double prob){
            _prob = prob;
        }

        public int getLineageCount(){
            return _lineages.size();
        }

        public void addProbability(double adds){
            _prob += adds;
        }
        /*
        public String toString(){
            String exp = "";
            for(int id: _lineages) {
                exp = exp + _gtClusters.get(id);
            }
            exp = exp + "/[";
            for(int i=0; i<_netNodeIndex.length; i++){
                exp = exp + _netNodeIndex[i];
                if(i!=_netNodeIndex.length-1){
                    exp = exp + ",";
                }
            }
            exp = exp + "]:" + _prob;
            return exp;
        }     */

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

