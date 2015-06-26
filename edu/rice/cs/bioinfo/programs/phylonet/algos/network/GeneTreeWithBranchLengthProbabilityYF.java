package edu.rice.cs.bioinfo.programs.phylonet.algos.network;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.integration.GTBranchLengthsIntegrationForSpeciesPhylogeny;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;


/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 5/10/12
 * Time: 10:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneTreeWithBranchLengthProbabilityYF {
    Map<Integer, Double> _netNode2time;
    //String[] _netTaxa;
    Set<NetNode> _articulateNodes;
    boolean _printDetail = false;
    int _netNodeNum;
    int _totalNodeNum;

    boolean _parallel = false;
    int _currentTreeID = 0;
    int _totalTree;

    Network _network;
    List<Tree> _gts;
    Map<String, List<String>> _species2alleles;

    //int _binNumber = 1;
    //int _sampleSize = 10000;


    public GeneTreeWithBranchLengthProbabilityYF(Network network, List<Tree> gts, Map<String, List<String>> species2alleles){
        _network = network;
        _gts = gts;
        _species2alleles = species2alleles;
        _totalTree = gts.size();
        processNetwork();
    }

    public void setParallel(boolean parallel){
        _parallel = parallel;
    }

    public void setTotalNodes(Set<NetNode> articulateNodes){
        _articulateNodes = articulateNodes;
    }
    

    public synchronized int getNextTreeID(int batchSize){
        //System.out.println("In calculation:" + Thread.currentThread().getId() + ":");
        _currentTreeID = _currentTreeID + batchSize;
        //System.out.println("In calculation:" + Thread.currentThread().getId() + ":" + _currentTreeID);
        return _currentTreeID - batchSize;
    }

    public void setPrintDetails(boolean p){
        _printDetail = p;
    }

    /*
    public void setIntegralParameter(int binNumber, int sampleSize){
        _binNumber = binNumber;
        _sampleSize = sampleSize;
    }
    */

    public void calculateGTDistribution(double[] resultProbs){
        calculateGTDistribution(1, resultProbs);
    }

    public void calculateGTDistribution(int batchSize, double[] resultProbs){
        //System.out.println(gts.size());
        //System.exit(0);
        int treeID = 0;
        if(_parallel){
            treeID = getNextTreeID(batchSize);
        }

        int count = 0;
        while(treeID < _totalTree){
            Tree gt = _gts.get(treeID);
            /*
            for(TNode node: gt.postTraverse()){
                if(!node.isRoot() && node.getParentDistance()<0){
                    resultProbs[treeID] = 0;
                    break;
                }
            }
            */

            resultProbs[treeID] = computeBinaryGTProbability(gt);
                /*
                if (!Trees.isBinary(gt)) {
                    GTBranchLengthsIntegrationForSpeciesPhylogeny integral = new GTBranchLengthsIntegrationForSpeciesPhylogeny(_network, gt, _species2alleles);
                    integral.setArticulateNodes(_articulateNodes);
                    resultProbs[treeID] = integral.computeLikelihoodWithIntegral(_binNumber, _sampleSize);
                } else {
                    resultProbs[treeID] = computeBinaryGTProbability(gt);
                }
                */


            //System.out.println(gt + ": " + gtProb);
            if(_printDetail){
                System.out.println();
                System.out.println(gt + ": " + resultProbs[treeID]);
            }

            if(_parallel){
                count++;
                if(count == batchSize) {
                    treeID = getNextTreeID(batchSize);
                    count = 0;

                }
                else{
                    treeID ++;
                }
            }
            else{
                treeID++;
            }
        }
    }

    private List<Configuration> computeAC(NetNode<Integer> node, String[] gtTaxa, Set<String> gtTaxaSet, HashMap<Tuple<Integer,Integer>, List<Configuration>> edge2ACminus, HashSet<IntArray> invalidConfigs, List<STITreeCluster<Double>> gtClusters){

        List<Configuration> CACs =  new ArrayList<Configuration>();
        if(node.isLeaf()){
            Configuration config = new Configuration();
            if(_species2alleles == null){
                if(gtTaxaSet.contains(node.getName())){
                    STITreeCluster cl = new STITreeCluster(gtTaxa);
                    cl.addLeaf(node.getName());
                    config.addLineage(gtClusters.indexOf(cl), gtClusters);
                }
            }
            else{
                for(String allele: _species2alleles.get(node.getName())){
                    if(gtTaxaSet.contains(allele)){
                        STITreeCluster cl = new STITreeCluster(gtTaxa);
                        cl.addLeaf(allele);
                        config.addLineage(gtClusters.indexOf(cl), gtClusters);
                    }
                }
            }
            config.setTotalProbability(1.0);
            CACs.add(config);
        }
        else{
            if(node.getOutdeg() == 1){
                Iterator<NetNode<Integer>> childNode = node.getChildren().iterator();
                Tuple<Integer,Integer> edge = new Tuple<Integer, Integer>(node.getData(),childNode.next().getData());
                CACs.addAll(edge2ACminus.remove(edge));
            }
            else{
                boolean isArticulateNode = _articulateNodes.contains(node.getData());
                Iterator<NetNode<Integer>> childNode = node.getChildren().iterator();
                Tuple<Integer, Integer> edge1 = new Tuple<Integer, Integer>(node.getData(),childNode.next().getData());
                List<Configuration> AC1 = edge2ACminus.remove(edge1);
                Tuple<Integer, Integer> edge2 = new Tuple<Integer, Integer>(node.getData(),childNode.next().getData());
                List<Configuration> AC2 = edge2ACminus.remove(edge2);
                for(Configuration config1: AC1){
                    if(invalidConfigs.contains(new IntArray(config1._netNodeIndex))){
                        continue;
                    }
                    for(Configuration config2: AC2){
                        if(invalidConfigs.contains(new IntArray(config2._netNodeIndex))){
                            continue;
                        }
                        if(config1.isCompatible(config2)){
                            Configuration mergedConfig = new Configuration(config1, config2);
                            if(mergedConfig._totalProb == 0){
                                continue;
                            }
                            if(isArticulateNode){
                                mergedConfig.clearNetNodeChoice();
                                int index = CACs.indexOf(mergedConfig);
                                if(index == -1){
                                    CACs.add(mergedConfig);
                                }
                                else{
                                    //System.out.println(mergedConfig);
                                    Configuration exist = CACs.get(index);
                                    exist.addTotalProbability(mergedConfig._totalProb);
                                }
                            }
                            else{
                                CACs.add(mergedConfig);
                            }
                        }
                    }
                }
            }
        }
        if(_printDetail){
            System.out.print("AC:[");
            for(Configuration config: CACs){
                System.out.print(config.toString(gtClusters)+" ");
            };
            System.out.println("]");
        }

        return CACs;
    }

    private double computeACMinus(NetNode<Integer> node, int netNodeIndex, HashMap<Tuple<Integer,Integer>, List<Configuration>> edge2ACminus, HashSet<IntArray> invalidConfigs, List<STITreeCluster<Double>> gtClusters, List<Configuration> CACs, boolean[][] R){
        //set AC- for a node
        double lowTau = _netNode2time.get(node.getData());
        double gtProb = 0;
        if(node.isRoot()){
            Configuration rootConfig = new Configuration();
            rootConfig.addLineage(gtClusters.size()-1, gtClusters);
            for(Configuration preConfig: CACs){
                Configuration postConfig = new Configuration(preConfig);
                calculateProbability(lowTau, Double.MAX_VALUE, postConfig, 1, gtClusters, R);
                if(_printDetail){
                    System.out.println(preConfig.toString(gtClusters) + " --> " + postConfig.toString(gtClusters));
                }
                gtProb += postConfig._totalProb;
            }
            return gtProb;
        }
        else if(node.isTreeNode()){
            double highTau = _netNode2time.get(node.getParents().iterator().next().getData());
            List<Configuration> ACminus = new ArrayList<Configuration>();
            Iterator<Configuration> it = CACs.iterator();
            while(it.hasNext()){
                Configuration preConfig = it.next();
                Configuration postConfig = new Configuration(preConfig);
                calculateProbability(lowTau, highTau, postConfig, 1, gtClusters, R);
                if(_printDetail){
                    System.out.println(preConfig.toString(gtClusters) + " --> " + postConfig.toString(gtClusters));
                }
                if(postConfig._totalProb == 0){
                    it.remove();
                    invalidConfigs.add(new IntArray(postConfig._netNodeIndex));
                }
                else{
                    ACminus.add(postConfig);
                }
            }

            if(_printDetail){
                System.out.print("ACMinus:[");
                for(Configuration config: ACminus){
                    System.out.print(config.toString(gtClusters)+" ");
                }
                System.out.println("]");
            }
            if(ACminus.size() == 0 && CACs.size()!=0){
                return -1;
            }
            Tuple<Integer, Integer> newEdge = new Tuple<Integer, Integer>(node.getParents().iterator().next().getData(), node.getData());
            edge2ACminus.put(newEdge, ACminus);

        }
        else {
            List<Configuration> newCACs1 = new ArrayList<Configuration>();
            List<Configuration> newCACs2 = new ArrayList<Configuration>();
            int configIndex = 1;
            for(Configuration config: CACs){
                int[] lineageArray = new int[config.getLineageCount()];
                int index = 0;
                for(int lineage: config._lineages){
                    lineageArray[index++] = lineage;
                }
                for(int i=0; i<=config.getLineageCount(); i++){
                    for(boolean[] selectedLineages: getSelected(config.getLineageCount(),i)){
                        Configuration newConfig1 = new Configuration();
                        Configuration newConfig2 = new Configuration();
                        index = 0;
                        for(int lin: lineageArray) {
                            if(selectedLineages[index]){
                                newConfig1.addLineage(lin, gtClusters);
                            }
                            else{
                                newConfig2.addLineage(lin, gtClusters);
                            }
                            index ++;
                        }

                        newConfig1.setNetNodeChoice(config._netNodeIndex);
                        newConfig1.setTotalProbability(config._totalProb);
                        newConfig1.addNetNodeChoice(netNodeIndex, configIndex);
                        newCACs1.add(newConfig1);

                        newConfig2.setNetNodeChoice(config._netNodeIndex);
                        newConfig2.setTotalProbability(1);
                        newConfig2.addNetNodeChoice(netNodeIndex, configIndex);
                        newCACs2.add(newConfig2);
                        configIndex ++;
                    }
                }
            }
            Iterator<NetNode<Integer>> parentIt = node.getParents().iterator();
            for(int i=0; i<2; i++){
                List<Configuration> ACminus = new ArrayList<Configuration>();
                List<Configuration> newCACs;
                if(i==0){
                    newCACs = newCACs1;
                }
                else{
                    newCACs = newCACs2;
                }

                NetNode<Integer> parentNode = parentIt.next();
                double highTau = _netNode2time.get(parentNode.getData());
                double gamma = node.getParentProbability(parentNode);
                Iterator<Configuration> it = newCACs.iterator();
                while(it.hasNext()){
                    Configuration preConfig = it.next();
                    Configuration postConfig = new Configuration(preConfig);
                    calculateProbability(lowTau, highTau, postConfig, gamma, gtClusters, R);
                    if(_printDetail){
                        System.out.println(preConfig.toString(gtClusters) + " --> " + postConfig.toString(gtClusters));
                    }
                    if(postConfig._totalProb == 0){
                        it.remove();
                        invalidConfigs.add(new IntArray(postConfig._netNodeIndex));
                    }
                    else{
                        ACminus.add(postConfig);
                    }
                }

                if(_printDetail){
                    System.out.print("ACminus to " + parentNode.getName()+ ":[ ");
                    for(Configuration config: ACminus){
                        System.out.print(config.toString(gtClusters)+" ");
                    }
                    System.out.println("]");
                }

                if(ACminus.size()==0 && newCACs.size()!=0 && gamma!=0){
                    return -1;
                }

                Tuple<Integer, Integer> newEdge = new Tuple<Integer, Integer>(parentNode.getData(), node.getData());
                edge2ACminus.put(newEdge, ACminus);

            }
        }
        return gtProb;
    }

    private double computeBinaryGTProbability(Tree gt){
        double gtProb = 0;
        String[] gtTaxa = gt.getLeaves();
        List<STITreeCluster<Double>> gtClusters = new ArrayList<STITreeCluster<Double>>();
        boolean[][] R = processGT(gt, gtTaxa, gtClusters);

        //computeR();
        //HashMap<NetNode, BitSet> node2events = computeLowestEvents(gt, network, species2alleles);

        HashSet<String> gtTaxaSet = new HashSet<String>();
        Collections.addAll(gtTaxaSet, gtTaxa);

        HashMap<Tuple<Integer,Integer>, List<Configuration>> edge2ACminus = new HashMap<Tuple<Integer,Integer>, List<Configuration>>();
        HashSet<IntArray> invalidConfigs = new HashSet<IntArray>();
        int netNodeIndex = 0;
        for(Object nodeO: Networks.postTraversal(_network)){
            NetNode<Integer> node = (NetNode<Integer>)nodeO;
            if(_printDetail){
                System.out.println();
                System.out.println("On node #" + node.getData() + " " + node.getName());
            }

            List<Configuration> CACs = computeAC(node, gtTaxa, gtTaxaSet, edge2ACminus, invalidConfigs, gtClusters);
            gtProb = computeACMinus(node, netNodeIndex, edge2ACminus, invalidConfigs, gtClusters, CACs, R);
            if(node.isNetworkNode()){
                netNodeIndex++;
            }
            if(gtProb == -1){
                break;
            }
        }
        return gtProb;
    }

    private void calculateProbability(double lowTau, double highTau, Configuration config, double gamma, List<STITreeCluster<Double>> gtClusters, boolean R[][]){
        boolean hasPrint = false;
        int u = config.getLineageCount();
        config._totalProb = config._totalProb * Math.pow(gamma, u);
        if(_printDetail){
            if(gamma != 1){
                System.out.print("*("+gamma+")");
                if(u != 1){
                    System.out.print("^"+u);
                }
                hasPrint = true;
            }

        }
        if(u == 0){
            if(_printDetail){
                if(hasPrint)
                    System.out.print("	:");
                else
                    System.out.print("1	:");
            }
            return;
        }
        double prob = 1;
        int index = 0;
        for(STITreeCluster<Double> gtCl: gtClusters){
            double coalTime = gtCl.getData();
            if(coalTime>=highTau){
                break;
            }
            if(coalTime>=lowTau){
                BitSet temp = (BitSet) config._coverage.clone();
                temp.and(gtCl.getCluster());

                if(temp.cardinality() != 0){ // not disjoint
                    //temp = (BitSet) config._coverage.clone();
                    //temp.and(gtCl.getCluster());
                    if(temp.equals(gtCl.getCluster())){  //contains
                        config.mergeCluster(index, R);
                        prob *= Math.exp((-1)*(coalTime-lowTau)*u*(u-1)/2);
                        //System.out.println(Math.exp((-1)*(coalTime-lowTau)*u*(u-1)/2));
                        //prob *= Math.exp((-1)*(coalTime-lowTau)*u*(u-1));
                        if(_printDetail){
                            System.out.print("*exp(-" + u*(u-1)/2 + "*" + (coalTime-lowTau) + ")");
                            //System.out.print("*exp(-" + u*(u-1)/2 + "*" + (coalTime-lowTau) + ")");
                            hasPrint = true;
                        }
                        lowTau = coalTime;
                        u--;
                    }
                    else{
                        if(_printDetail){
                            System.out.print("*0	:");
                        }
                        config._totalProb = 0;
                        return;
                    }
                }
            }
            index ++;
        }

        if(u != 1){
            prob *= Math.exp((-1)*(highTau-lowTau)*u*(u-1)/2);
            //prob *= Math.exp((-1)*(highTau-lowTau)*u*(u-1));
            if(_printDetail){
                System.out.print("*exp(-" + u*(u-1)/2 + "*" + (highTau-lowTau) + ")");
                //System.out.print("*exp(-" + u*(u-1) + "*" + (highTau-lowTau) + ")");
                hasPrint = true;
            }
        }
        if(_printDetail){
            if(hasPrint)
                System.out.print("	:");
            else
                System.out.print("1	:");
        }
        prob = Math.max(0, prob);
        config._totalProb = config._totalProb * prob;
    }


    /**
     * The function is to calculate the _R matrix for the given tree.
     * Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
     */
    private boolean[][] computeR(List<STITreeCluster<Double>> gtClusters){
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


    private void processNetwork(){
        removeBinaryNodes(_network);
        _netNodeNum = 0;
        _totalNodeNum = 0;
        _netNode2time = new HashMap<Integer, Double>();
        List<String> taxa = new ArrayList<String>();
        for(Object nodeO: Networks.postTraversal(_network)){
            NetNode<Integer> node = (NetNode<Integer>)nodeO;
            double minTime = Double.MAX_VALUE;
            if(node.isLeaf()){
                taxa.add(node.getName());
                minTime = 0;

            }else{
                if(node.isNetworkNode()){
                    _netNodeNum++;
                }
                for (NetNode<Integer> child : node.getChildren()) {
                    minTime = Math.min(minTime, _netNode2time.get(child.getData()) + child.getParentDistance(node));
                }
            }
            _netNode2time.put(_totalNodeNum, minTime);
            node.setData(_totalNodeNum++);
        }
        if(_articulateNodes==null){
            _articulateNodes = Networks.getLowestArticulationNodes(_network);
        }

    }


    private boolean[][] processGT(Tree gt, String[] gtTaxa, List<STITreeCluster<Double>> gtClusters){
        Map<TNode, STITreeCluster> map = new HashMap<TNode, STITreeCluster>();
        for (TNode node : gt.postTraverse()) {
            STITreeCluster cl = new STITreeCluster(gtTaxa);
            double maxTime = Double.MIN_VALUE;
            int index = -1;
            if (node.isLeaf()) {
                cl.addLeaf(node.getName());
                ((STINode<Double>)node).setData(0.0);
                maxTime = 0;
            }
            else {
                for(TNode child : node.getChildren()) {
                    cl = cl.merge(map.get(child));
                    maxTime = Math.max(maxTime, ((STINode<Double>)child).getData()+child.getParentDistance());
                }
                ((STINode<Double>)node).setData(maxTime);
                for(index=0; index<gtClusters.size(); index++){
                    double clTime = gtClusters.get(index).getData();
                    if(clTime == 0 || clTime > maxTime){
                        break;
                    }
                }

            }
            map.put(node, cl);

            STITreeCluster<Double> newCl = new STITreeCluster<Double>(cl);
            newCl.setData(maxTime);
            if(index == -1){
                gtClusters.add(newCl);
            }
            else{
                gtClusters.add(index, newCl);
            }
        }
        boolean[][] R = computeR(gtClusters);
        return R;
    }



    private void removeBinaryNodes(Network<Integer> net)
    {
        // Find all binary nodes.
        List<NetNode<Integer>> binaryNodes = new LinkedList<NetNode<Integer>>();
        for (NetNode<Integer> node : net.bfs()) {
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




/*
    private void computeArticulateNodes(Network<Integer> net){
        //List<Integer> leaves = new ArrayList<Integer>();
        List<NetNode> allTotalNodes = new ArrayList<NetNode>();
        _articulateNodes = new HashSet<NetNode>();
        for(NetNode<Integer> node: Networks.postTraversal(net)){
            if(node.isLeaf()){
                //leaves.add(id);
                allTotalNodes.add(node);
            }

            else if(node.isRoot()){
                boolean ftotal = true;
                for(NetNode<Integer> child: node.getChildren()){
                    if(!allTotalNodes.contains(child.getData())){
                        ftotal = false;
                        break;
                    }
                }
                if(!ftotal){
                    _articulateNodes.add(node);
                }

            }
            else if(node.isTreeNode()){
                boolean ftotal = true;
                for(NetNode<Integer> child: node.getChildren()){
                    if(!allTotalNodes.contains(child.getData())){
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
                    boolean disconnect = isValidNetwork(net);
                    parent.adoptChild(node, distance);
                    if (disconnect) {
                        _articulateNodes.add(node);
                        allTotalNodes.add(node);
                    }
                }

            }
        }
    }

    private boolean isValidNetwork(Network<Integer> net){
        Set<NetNode> visited = new HashSet<NetNode>();
        Set<NetNode> seen = new HashSet<NetNode>();
        for(NetNode<Integer> node: net.bfs()){
            if(node.getIndeg()==1 && node.getOutdeg()==1) return false;
            visited.add(node);
            for(NetNode<Integer> parent: node.getParents()){
                seen.add(parent);
            }
            for(NetNode<Integer> child: node.getChildren()){
                seen.add(child);
            }
        }
        return visited.size()==seen.size();
    }


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

    private class Configuration{
        private HashSet<Integer> _lineages;
        private double _totalProb;
        private BitSet _coverage;
        int[] _netNodeIndex;
        //private Set<Integer> _coalEvents;

        public Configuration(){
            _lineages = new HashSet<Integer>();
            _netNodeIndex = new int[_netNodeNum];
            Arrays.fill(_netNodeIndex, 0);
            _coverage = new BitSet();
        }

        public Configuration(Configuration config){
            _lineages = (HashSet)config._lineages.clone();
            _totalProb = config._totalProb;
            _netNodeIndex = config._netNodeIndex.clone();
            _coverage = (BitSet)(config._coverage.clone());
        }

        public Configuration(Configuration config1, Configuration config2){
            _lineages = (HashSet)config1._lineages.clone();
            _lineages.addAll(config2._lineages);
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
            _coverage = (BitSet)(config1._coverage.clone());
            _coverage.or(config2._coverage);
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


        public void addLineage(int index, List<STITreeCluster<Double>> gtClusters){
            _lineages.add(index);
            _coverage.or(gtClusters.get(index).getCluster());
        }


        public void mergeCluster(int e, boolean[][] R){
            Iterator<Integer> lineageIt = _lineages.iterator();
            while(lineageIt.hasNext()) {
                int lin = lineageIt.next();
                if(R[e][lin]){
                    lineageIt.remove();
                }
            }
            _lineages.add(e);

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

        public String toString(List<STITreeCluster<Double>> gtClusters){
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

    class IntArray{
        int[] _data;
        int _code;

        public IntArray(int[] d){
            _data = d;
            for(int e: _data){
                _code += e;
            }
        }

        public boolean equals(Object o) {
            if(!(o instanceof IntArray)){
                return false;
            }

            IntArray ia = (IntArray) o;
            return Arrays.equals(ia._data, _data);
        }

        public int hashCode(){
            return _code;
        }
    }

}
