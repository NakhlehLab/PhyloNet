package edu.rice.cs.bioinfo.library.phylogenetics.scoring.network;

import edu.rice.cs.bioinfo.library.phylogenetics.*;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func4;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/9/12
 * Time: 5:53 PM
 * To change this template use File | Settings | File Templates.
 */
public class MDCOnNetworkYF {
    List<STITreeCluster> _gtClusters;
    int[][] _gtclConstitution;
    // Iterable<String> gtTaxa;
    String[] _netTaxa;
    // Map<NetNode, Integer> netNode2id;
    BitSet _totalCoverNode;
    boolean _printDetail;
    int _netNodeNum;
    int _totalNodeNum;
    int[][] _totalNetNodeLinNum;


    public void setPrintDetails(boolean p){
        _printDetail = p;
    }

    public double[] getHybridProbabilities(){
        double[] probabilities = new double[_totalNetNodeLinNum.length];
        int index = 0;
        for(int[] lineageNum: _totalNetNodeLinNum){
            double total = lineageNum[0]+lineageNum[1];
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

    public <N1,E1,N2,E2> List<Integer> countExtraCoal(Graph<N1,E1> network, Iterable<Graph<N2,E2>> gts, Map<String, List<String>> species2alleles,
                                                      Func1<N1,String> getNetworkNodeLabel, Func1<N2,String> getGTNodeLabel,
                                                      Func2<GraphReadOnly<N1,E1>, E1, Double> getNetworkDistance,
                                                      Func2<GraphReadOnly<N1,E1>, E1, Double> getNetworkProbability,
                                                      Func2<GraphReadOnly<N2,E2>, E2, Double> getGTDistance,
                                                      Func2<GraphReadOnly<N2,E2>, E2, Double> getGTProbability,
                                                      Func4<N1,N1,Double,Double,E1> makeNetworkEdge,
                                                      Func4<N2,N2,Double,Double,E2> makeGTEdge){
        List<Integer> xlList = new ArrayList<Integer>();
        Map<N1, Integer> netNode2id = new HashMap<N1,Integer>();
        processNetwork(network, getNetworkNodeLabel, netNode2id, getNetworkDistance, getNetworkProbability, makeNetworkEdge);
        //System.out.println(gts.size());
        //System.exit(0);
        IsLeaf isLeaf = new IsLeaf();
        for(Graph<N2,E2> gt: gts){
            int xl = Integer.MAX_VALUE;
            ArrayList<String> gtTaxa = new ArrayList<String>();
            for(N2 node :  new GetLeafs<N2>().execute(gt))
            {
                gtTaxa.add(getGTNodeLabel.execute(node));
            }
            if(!checkLeafAgreement(species2alleles, gtTaxa)){
                throw new RuntimeException("Gene tree " + gt + " has leaf that the network doesn't have.");
            }
            _gtClusters = new ArrayList<STITreeCluster>();
            processGT(gt, gtTaxa, getGTNodeLabel, getGTDistance, getGTProbability, makeGTEdge);

            HashMap<BitSet, List<Configuration>> edge2ACminus = new HashMap<BitSet, List<Configuration>>();

            int netNodeIndex = 0;
            for(N1 node: walkNetwork(network,null)){
                //System.out.println(edge2ACminus);
                if(_printDetail){
                    System.out.println();
                    System.out.println("On node #" + netNode2id.get(node) + " " + node);
                }
                ArrayList<Configuration> CACs =  new ArrayList<Configuration>();

                //set AC for a node
                if(isLeaf.execute(network, node)){
                    Configuration config = new Configuration(gtTaxa);
                    String nodeLabel = getNetworkNodeLabel.execute(node);
                    if(species2alleles == null){

                        if(gtTaxa.contains(nodeLabel)){
                            STITreeCluster cl = new STITreeCluster(gtTaxa.toArray(new String[0]));
                            cl.addLeaf(nodeLabel);
                            config.addLineage(cl);
                        }
                    }
                    else{
                        for(String allele: species2alleles.get(nodeLabel)){
                            if(gtTaxa.contains(allele)){
                                STITreeCluster cl = new STITreeCluster(gtTaxa.toArray(new String[0]));
                                cl.addLeaf(allele);
                                config.addLineage(cl);
                            }
                        }
                    }
                    config.setExtraLineage(0);
                    CACs.add(config);
                }
                else{
                    if(new GetOutDegree().execute(network, node) == 1){
                        Iterator<N1> childNodes = new GetDirectSuccessors().execute(network, node).iterator();
                        BitSet edge = new BitSet();
                        edge.set(netNode2id.get(node));
                        edge.set(netNode2id.get(childNodes.next()));
                        CACs.addAll(edge2ACminus.remove(edge));
                    }
                    else{
                        //TODO only on binary nodes

                        Iterator<N1> childNodes = new GetDirectSuccessors().execute(network, node).iterator();
                        BitSet edge1 = new BitSet();
                        edge1.set(netNode2id.get(node));
                        edge1.set(netNode2id.get(childNodes.next()));
                        List<Configuration> AC1 = edge2ACminus.remove(edge1);
                        BitSet edge2 = new BitSet();
                        edge2.set(netNode2id.get(node));
                        edge2.set(netNode2id.get(childNodes.next()));
                        List<Configuration> AC2 = edge2ACminus.remove(edge2);
                        int minimum = Integer.MAX_VALUE;
                        boolean totalCover = _totalCoverNode.get(netNode2id.get(node));
                        BitSet[][] netNodeLineages = null;
                        if(totalCover){
                            netNodeLineages = new BitSet[_netNodeNum][2];
                            for(int i=0; i< _netNodeNum; i++){
                                for(int j=0; j<2; j++){
                                    netNodeLineages[i][j] = new BitSet(_gtClusters.size());
                                }
                            }
                        }
                        //System.out.println(AC1.size()+"	with	"+AC2.size());
                        for(Configuration config1: AC1){
                            for(Configuration config2: AC2){
                                if(config1.isCompatible(config2)){
                                    Configuration mergedConfig = new Configuration(config1, config2);
                                    if(totalCover){
                                        if(minimum > mergedConfig._xl){
                                            minimum = mergedConfig._xl;
                                            CACs.clear();
                                            CACs.add(mergedConfig);
                                            for(int i=0; i< _netNodeNum; i++){
                                                for(int j=0; j<2; j++){
                                                    netNodeLineages[i][j] = (BitSet)(mergedConfig._netNodeLineages[i][j].clone());
                                                }
                                            }
                                        }
                                        else if(minimum == mergedConfig._xl){
                                            for(Configuration optimal: CACs){
                                                if(optimal._coverage.equals(mergedConfig._coverage)){
                                                    if(optimal.getLineageCount() > mergedConfig.getLineageCount()){
                                                        CACs.remove(optimal);
                                                        CACs.add(mergedConfig);
                                                    }
                                                }
                                                else{
                                                    CACs.add(mergedConfig);
                                                }
                                            }

                                            for(int i=0; i< _netNodeNum; i++){
                                                for(int j=0; j<2; j++){
                                                    netNodeLineages[i][j].and(mergedConfig._netNodeLineages[i][j]);
                                                    //System.out._printDetail(netNodeLineages[i][0].cardinality() + "/" + netNodeLineages[i][1].cardinality() + "   ");
                                                }
                                            }
                                        }

                                        /*
                                                  else if(minimum == mergedConfig._xl){
                                                      CACs.get(0).addNetNodeLineageNum(mergedConfig);
                                                  }
                                                  */
                                    }
                                    else{

                                        CACs.add(mergedConfig);
                                        /*
                                                  int index = CACs.indexOf(mergedConfig);
                                                  if(index == -1){
                                                      CACs.add(mergedConfig);
                                                  }
                                                  else{

                                                      Configuration exist = CACs.get(index);
                                                      exist.setExtraLineage(Math.min(mergedConfig._xl, exist._xl));
                                                  }
                                                  */
                                    }
                                }
                            }
                        }
                        if(totalCover){
                            CACs.get(0).setNetNodeLineageNum(netNodeLineages);
                        }
                    }
                }
                if(_printDetail){
                    System.out.println("AC:" + CACs.toString());
                }
                //set AC- for a node
                int nodeInDegree = new GetInDegree<N1,E1>().execute(network, node);
                if(nodeInDegree == 0){
                    if(CACs.size()!=1){
                        System.err.println("Error");
                    }
                    Configuration optimalConfig = CACs.get(0);
                    xl = optimalConfig._xl;
                    xlList.add(xl);

                    for(int i=0; i< _netNodeNum; i++){
                        //System.out._printDetail(optimalConfig._netNodeLineages[i][0].cardinality() + "/" + optimalConfig._netNodeLineages[i][1].cardinality() + "   ");
                        _totalNetNodeLinNum[i][0] += optimalConfig._netNodeLineages[i][0].cardinality();
                        _totalNetNodeLinNum[i][1] += optimalConfig._netNodeLineages[i][1].cardinality();
                    }
                    //System.out.println();
                }
                else if(nodeInDegree == 1){
                    List<Configuration> ACminus = new ArrayList<Configuration>();
                    for(Configuration config: CACs){
                        for(int i=0; i< _gtclConstitution.length; i++){
                            config.mergeCluster(i, _gtclConstitution[i]);
                        }
                        config.addExtraLineage(Math.max(0, config.getLineageCount()-1));
                        ACminus.add(config);
                        /*
                        int index = ACminus.indexOf(config);
                        if(index == -1){
                            ACminus.add(config);
                        }
                        else{
                            //System.out.println(mergedConfig);
                            System.out.println("here");
                            Configuration exist = ACminus.get(index);
                            exist.setExtraLineage(Math.min(config._xl,exist._xl));
                        }
                        */
                    }

                    BitSet newEdge = new BitSet();
                    newEdge.set(netNode2id.get(node));
                    newEdge.set(netNode2id.get(new GetDirectPredecessors<N1, E1>().execute(network, node).iterator().next()));
                    edge2ACminus.put(newEdge, ACminus);

                    if(_printDetail){
                        System.out.println("ACminus: " + ACminus);
                    }
                }
                else {
                    List<Configuration> ACminus1 = new ArrayList<Configuration>();
                    List<Configuration> ACminus2 = new ArrayList<Configuration>();
                    int configIndex = 1;
                    for(Configuration config: CACs){
                        int numLineage = config.getLineageCount();
                        for(int i=0; i<=numLineage; i++){
                            for(BitSet selectedLineages: getSelected(numLineage,i)){
                                Configuration newConfig1 = new Configuration(gtTaxa);
                                Configuration newConfig2 = new Configuration(gtTaxa);
                                newConfig1.setNetNodeLineageNum(config._netNodeLineages);
                                newConfig2.setNetNodeLineageNum(config._netNodeLineages);

                                int index = 0;
                                for (int k = config._lineages.nextSetBit(0); k >= 0; k = config._lineages.nextSetBit(k+1)) {
                                    if(selectedLineages.get(index)){
                                        newConfig1.addLineage(k);
                                        newConfig1.addNetNodeLineageNum(netNodeIndex, 0, k);
                                    }
                                    else{
                                        newConfig2.addLineage(k);
                                        newConfig2.addNetNodeLineageNum(netNodeIndex, 1, k);
                                    }
                                    index ++;

                                }

                                newConfig1.setNetNodeChoice(config._netNodeIndex);
                                //newConfig1.setNetNodeLineageNum(config._netNodeLineages);
                                newConfig1.setExtraLineage(config._xl);
                                newConfig1.addNetNodeChoice(netNodeIndex, configIndex);
                                //newConfig1.addNetNodeLineageNum(netNodeIndex, i, 0);
                                ACminus1.add(newConfig1);

                                newConfig2.setNetNodeChoice(config._netNodeIndex);
                                //newConfig2.setNetNodeLineageNum(config._netNodeLinNum);
                                newConfig2.setExtraLineage(0);
                                newConfig2.addNetNodeChoice(netNodeIndex, configIndex);
                                //newConfig2.addNetNodeLineageNum(netNodeIndex, 0, numLineage-i);
                                ACminus2.add(newConfig2);
                                configIndex ++;
                            }
                        }
                    }
                    /*
                         System.out.println(ACminus1);
                         System.out.println(ACminus2);
                         System.out.println();
                         */
                    Iterator<N1> it = new GetDirectPredecessors<N1,E1>().execute(network, node).iterator();
                    for(int i=0; i<2; i++){
                        List<Configuration> ACminus;
                        if(i==0){
                            ACminus = ACminus1;
                        }
                        else{
                            ACminus = ACminus2;
                        }
                        N1 parentNode = it.next();
                        //TODO
                        //if(getNetworkDistance.execute(network, network.getEdge(parentNode,node)) != 0.0){
                            for(Configuration config: ACminus){
                                config.addExtraLineage(Math.max(0, config.getLineageCount()-1));
                            }
                        //}
                        BitSet newEdge = new BitSet();
                        newEdge.set(netNode2id.get(node));
                        newEdge.set(netNode2id.get(parentNode));
                        edge2ACminus.put(newEdge, ACminus);

                        if(_printDetail){
                            System.out.println("ACminus to " + parentNode + ": " + ACminus);
                        }
                    }
                    netNodeIndex ++;
                }

            }

        }
        return xlList;
    }


    private boolean checkLeafAgreement(Map<String, List<String>> species2alleles, Iterable<String> gtTaxa){
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
            /*
               if(speciesSet.size()!=species2alleles.size()){
                   return false;
               }
               */
        }
        return true;
    }


    private <N,E> void processNetwork(Graph<N,E> net, Func1<N,String> getNetworkNodeLabel, Map<N, Integer> netNode2id,
                                      Func2<GraphReadOnly<N,E>, E, Double> getDistance,
                                      Func2<GraphReadOnly<N,E>, E, Double> getProbability,
                                      Func4<N,N,Double,Double,E> makeEdge){
        removeBinaryNodes(net, getDistance, getProbability, makeEdge);
        _netNodeNum = 0;
        _totalNodeNum = 0;
        List<String> taxa = new ArrayList<String>();
        IsLeaf isLeaf = new IsLeaf();
        GetInDegree<N,E> inDegree = new GetInDegree<N, E>();
        for(N node: net.getNodes()){
            netNode2id.put(node, _totalNodeNum++);
            if(isLeaf.execute(net, node)){
                taxa.add(getNetworkNodeLabel.execute(node));
            }else if(inDegree.execute(net, node) == 2){
                _netNodeNum++;
            }
        }
        _netTaxa = (String[]) taxa.toArray(new String[0]);
        _totalNetNodeLinNum = new int[_netNodeNum][2];
        computeNodeCoverage(net, getDistance, netNode2id);
    }


    private <N,E> void processGT(Graph<N,E> gt, ArrayList<String> gtTaxa, Func1<N,String> getGTNodeLabel,
                                 Func2<GraphReadOnly<N,E>, E, Double> getDistance, Func2<GraphReadOnly<N,E>, E, Double> getProbability,
                                 Func4<N,N,Double,Double,E> makeEdge){
        removeBinaryNodes(gt, getDistance, getProbability, makeEdge);
        Map<N, STITreeCluster> map = new HashMap<N, STITreeCluster>();
        IsLeaf isLeaf = new IsLeaf();
        GetDirectSuccessors<N,E> getDirectSuccessors = new GetDirectSuccessors<N, E>();
        for (N node : new GetNodesPostOrder<N,E>().execute(gt)) {
            STITreeCluster cl = new STITreeCluster(gtTaxa.toArray(new String[0]));
            if (isLeaf.execute(gt, node)) {
                cl.addLeaf(getGTNodeLabel.execute(node));
                //_gtClusters.add(cl);
            }
            else {

                for(N child : getDirectSuccessors.execute(gt, node)) {
                    cl = cl.merge(map.get(child));
                }

                int i = 0;
                for(; i< _gtClusters.size(); i++){
                    if(_gtClusters.get(i).getClusterSize() > cl.getClusterSize()){
                        break;
                    }
                }
                _gtClusters.add(i, cl);
            }
            map.put(node, cl);
        }

        for(int i=0; i<gtTaxa.size(); i++){
            STITreeCluster cl = new STITreeCluster(gtTaxa.toArray(new String[0]));
            cl.addLeaf(gtTaxa.get(i));
            _gtClusters.add(cl);
        }

        _gtclConstitution = new int[_gtClusters.size()-gtTaxa.size()][2];
        for (Map.Entry<N, STITreeCluster> entry: map.entrySet()) {
            N pNode = entry.getKey();
            if(! isLeaf.execute(gt, entry.getKey())){
                int parent = _gtClusters.indexOf(map.get(pNode));
                int index = 0;
                for(N child : getDirectSuccessors.execute(gt, pNode)) {
                    _gtclConstitution[parent][index++] = _gtClusters.indexOf(map.get(child));
                }
            }
        }

    }


    private <N,E> void removeBinaryNodes(Graph<N,E> net, Func2<GraphReadOnly<N,E>, E, Double> getDistance,
                                         Func2<GraphReadOnly<N,E>, E, Double> getProbability,
                                         Func4<N,N,Double,Double,E> makeEdge)
    {
        // Find all binary nodes.
        List<N> binaryNodes = new LinkedList<N>();

        GetInDegree<N,E> gid = new GetInDegree<N, E>();
        GetOutDegree<N,E> god = new GetOutDegree<N, E>();
        for (N node : net.getNodes()) {
            if (god.execute(net,node) == 1 && gid.execute(net,node) == 1) {
                binaryNodes.add(node);
            }
        }

        GetDirectSuccessors<N,E> gs = new GetDirectSuccessors<N,E>();
        GetDirectPredecessors<N,E> gp = new GetDirectPredecessors<N, E>();
        // Remove them.
        for (N node : binaryNodes) {
            N child = gs.execute(net, node).iterator().next();	// Node's only child.
            if( gid.execute(net, child) != 1){
                continue;
            }
            N parent = gp.execute(net, node).iterator().next();	// Node's only parent.
            E parentEdge = net.getEdge(parent, node);
            E childEdge = net.getEdge(node, child);
            double distance = getDistance.execute(net, parentEdge) + getDistance.execute(net, childEdge);
            double gamma = getProbability.execute(net, parentEdge) + getProbability.execute(net, childEdge);
            net.removeNode(node);
            E edge = makeEdge.execute(parent,child,distance, gamma);
            net.addEdge(edge);
        }
    }


    private <N,E> List<N> walkNetwork(GraphReadOnly<N,E> net, N root){

        Stack<N> stack = new Stack<N>();
        List<N> searchedNodes = new ArrayList<N>();
        if(root == null){
            root = new FindRoot<N>().execute(net);
        }
        stack.push(root);
        Map<N, Integer> node2index = new HashMap<N, Integer>();
        node2index.put(root, 0);

        GetOutDegree<N,E> getOutDegree = new GetOutDegree<N, E>();
        GetDirectSuccessors<N,E> getDirectSuccessors = new GetDirectSuccessors<N, E>();

        while(!stack.isEmpty()){
            N topNode = stack.peek();
            int index = node2index.get(topNode);
            if(index == getOutDegree.execute(net,topNode)){
                searchedNodes.add(stack.pop());
            }
            else{
                Iterator<N> it = getDirectSuccessors.execute(net,topNode).iterator();
                for(int i=0; i<index; i++){
                    it.next();
                }
                N child = it.next();
                if(searchedNodes.contains(child)){
                    node2index.put(topNode, index + 1);
                }
                else{
                    stack.push(child);
                    node2index.put(child, 0);
                }
            }
        }

        return searchedNodes;
    }


    private <N,E> void computeNodeCoverage(Graph<N,E> net, Func2<GraphReadOnly<N,E>, E, Double> getDistance, Map<N, Integer> netNode2id){
        _totalCoverNode = new BitSet(_totalNodeNum);
        GetInDegree<N,E> getInDegree = new GetInDegree<N, E>();
        IsLeaf isLeaf = new IsLeaf();
        GetDirectPredecessors<N,E> getParents = new GetDirectPredecessors<N, E>();
        N root = new FindRoot<N>().execute(net);
        for(N trNode: net.getNodes()){
            int inDegree = getInDegree.execute(net, trNode);
            if(inDegree > 1){
                continue;
            }

            if (inDegree == 0) {
                _totalCoverNode.set(netNode2id.get(trNode), true);
            }
            else if(!isLeaf.execute(net, trNode)){
                N parent = getParents.execute(net, trNode).iterator().next();
                //double distance = getDistance.execute(net, parent);
                E edge = net.getEdge(parent, trNode);
                net.removeEdge(edge);
                boolean disconnect = isValidNetwork(net, netNode2id, root);
                net.addEdge(edge);
                if (disconnect){
                    _totalCoverNode.set(netNode2id.get(trNode), true);
                }
            }
        }
    }

    private List<BitSet> getSelected(int n, int m){
        List<BitSet> selectedList = new ArrayList<BitSet>();
        int[] order = new int[m+1];
        for(int i=0; i<=m; i++){
            order[i] = i-1;
        }
        int k = m;
        boolean flag = true;
        while(order[0] == -1){
            if(flag){
                BitSet bs = new BitSet(n);
                for(int i=1; i<=m; i++){
                    bs.set(order[i]);
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


    private <N,E> boolean isValidNetwork(GraphReadOnly<N,E> net, Map<N, Integer> netNode2id, N root){
        BitSet visited = new BitSet();
        BitSet seen = new BitSet();

        GetDirectPredecessors<N,E> getDirectPredecessors = new GetDirectPredecessors<N, E>();
        GetDirectSuccessors<N,E> getDirectSuccessors = new GetDirectSuccessors<N, E>();

        for(N node: walkNetwork(net, root)){
            //System.out.println("Visiting:" + node);
            visited.set(netNode2id.get(node), true);
            for(N parent: getDirectPredecessors.execute(net, node)){
                //System.out.println("Seeing parent:" + parent);
                seen.set(netNode2id.get(parent), true);
            }
            for(N child: getDirectSuccessors.execute(net, node)){
                //System.out.println("Seeing child:" + child);
                seen.set(netNode2id.get(child), true);
            }
        }
        return visited.cardinality()==seen.cardinality();
    }

    private class Configuration{
        private STITreeCluster _coverage;
        private BitSet _lineages;
        private int _xl;
        int[] _netNodeIndex;
        //int[][] _netNodeLinNum;
        BitSet[][] _netNodeLineages;

        public Configuration(List<String> gtTaxa){
            _lineages = new BitSet(_gtClusters.size());
            _netNodeIndex = new int[_netNodeNum];
            //_netNodeLinNum = new int[_netNodeNum][2];
            _netNodeLineages = new BitSet[_netNodeNum][2];
            for(int i=0; i< _netNodeNum; i++){
                for(int j=0; j<2; j++){
                    _netNodeLineages[i][j] = new BitSet(_gtClusters.size());
                }
            }
            _coverage = new STITreeCluster(gtTaxa.toArray(new String[0]));
            Arrays.fill(_netNodeIndex, 0);
        }


        public Configuration(Configuration config1, Configuration config2){
            _lineages = (BitSet)config1._lineages.clone();
            _lineages.or(config2._lineages);
            _xl = config1._xl + config2._xl;
            _coverage = new STITreeCluster(config1._coverage);
            _coverage = _coverage.merge(config2._coverage);
            _netNodeIndex = new int[_netNodeNum];
            //_netNodeLinNum = new int[_netNodeNum][2];
            _netNodeLineages = new BitSet[_netNodeNum][2];
            for(int i=0; i< _netNodeNum; i++){
                for(int j=0; j<2; j++){
                    _netNodeLineages[i][j] = new BitSet(_gtClusters.size());
                }
            }
            for(int i=0; i< _netNodeNum; i++){
                if(config1._netNodeIndex[i] == config2._netNodeIndex[i]){
                    _netNodeIndex[i] = config1._netNodeIndex[i];
                }
                else{
                    _netNodeIndex[i] = Math.max(config1._netNodeIndex[i], config2._netNodeIndex[i]);
                }
                //_netNodeLinNum[i][0] = config1._netNodeLinNum[i][0] + config2._netNodeLinNum[i][0];
                //_netNodeLinNum[i][1] = config1._netNodeLinNum[i][1] + config2._netNodeLinNum[i][1];
                _netNodeLineages[i][0] = (BitSet)(config1._netNodeLineages[i][0].clone());
                _netNodeLineages[i][0].or(config2._netNodeLineages[i][0]);
                _netNodeLineages[i][1] = (BitSet)(config1._netNodeLineages[i][1].clone());
                _netNodeLineages[i][1].or(config2._netNodeLineages[i][1]);
            }
        }


        public boolean isCompatible(Configuration config){
            boolean compatible = true;
            if(!_coverage.isDisjoint(config._coverage)){
                compatible = false;
            }
            for(int i=0; i< _netNodeNum; i++){
                if(_netNodeIndex[i] != config._netNodeIndex[i] && _netNodeIndex[i]!=0 && config._netNodeIndex[i]!=0){
                    compatible = false;
                }
            }
            return compatible;
        }

        public void addLineage(STITreeCluster cl){
            _lineages.set(_gtClusters.indexOf(cl));
            if(_coverage == null){
                _coverage = new STITreeCluster(cl);
            }
            else{
                _coverage = _coverage.merge(cl);
            }
        }

        public void addLineage(int index){
            _lineages.set(index);
            if(_coverage == null){
                _coverage = new STITreeCluster(_gtClusters.get(index));
            }
            else{
                _coverage = _coverage.merge(_gtClusters.get(index));
            }
        }


        public void mergeCluster(int coalTo, int[] coalFrom){
            if(_lineages.get(coalFrom[0]) && _lineages.get(coalFrom[1])){
                _lineages.set(coalFrom[0], false);
                _lineages.set(coalFrom[1], false);
                _lineages.set(coalTo, true);
            }
        }


        public void setExtraLineage(int xl){
            _xl = xl;
        }

        public void addExtraLineage(int xl){
            _xl += xl;
        }

        public int getLineageCount(){
            return _lineages.cardinality();
        }


        public String toString(){
            String exp = "";
            for (int i = _lineages.nextSetBit(0); i >= 0; i = _lineages.nextSetBit(i+1)) {
                exp = exp + _gtClusters.get(i);
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
            for(int i=0; i<_netNodeIndex.length; i++){
                exp = exp + _netNodeIndex[i];
                if(i!=_netNodeLineages.length-1){
                    exp = exp + ",";
                }
            }
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

        /*
          public void addNetNodeLineageNum1(Configuration config){
              for(int i=0; i<_netNodeNum; i++){
                  _netNodeLinNum[i][0] += config._netNodeLinNum[i][0];
                  _netNodeLinNum[i][1] += config._netNodeLinNum[i][1];
              }
          }
          */
        public void setNetNodeLineageNum(BitSet[][] lineageNum){
            for(int i=0; i< _netNodeNum; i++){
                for(int j=0; j<2; j++){
                    _netNodeLineages[i][j] = (BitSet)(lineageNum[i][j].clone());
                }
            }
            /*
               int index = 0;
               for(int[] num: linNum){
                   _netNodeLinNum[index++] = num.clone();
               }
               */
        }

        public boolean equals(Object o) {
            if(!(o instanceof Configuration)){
                return false;
            }

            Configuration config = (Configuration) o;
            if (config._coverage.equals(_coverage) && config._lineages.equals(_lineages) && Arrays.equals(config._netNodeIndex,_netNodeIndex)) {
                return true;
            }
            else {
                return false;

            }
        }

        public int hashCode(){
            return _coverage.hashCode()+_lineages.hashCode();
        }

    }

}
