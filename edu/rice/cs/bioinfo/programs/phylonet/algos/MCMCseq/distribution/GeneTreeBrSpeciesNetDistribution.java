package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.distribution;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.TemporalConstraints;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;

import java.io.StringReader;
import java.util.*;

/**
 * Created by wendingqiao on 2/25/16.
 */

/**
 * create this class only when the network is changed
 * this class will support two kinds of likelihood calculation
 *      1. p(a list of gene trees | network) when network is changed
 *      2. p(one gene tree | network) when that gene tree is changed
 */
public class GeneTreeBrSpeciesNetDistribution {

    Network<NetNodeInfo> _network;
    int _netNodeNum; // number of reticulation nodes
    int _totalNodeNum;
    Set<NetNode<NetNodeInfo>> _articulateNodes;
    Map<String, List<String>> _species2alleles;
    private boolean _onePopSize;

    public GeneTreeBrSpeciesNetDistribution(Network network, Map<String, List<String>> species2alleles){
        this._network = network;
        this._species2alleles = species2alleles;
        removeBinaryNodes(_network);
        _netNodeNum = 0;
        _totalNodeNum = 0;

        for(NetNode<NetNodeInfo> node: Networks.postTraversal(_network)){
            if(node.isNetworkNode()) _netNodeNum++;
            node.getData().setIndex(_totalNodeNum++);
        }

        if(_articulateNodes==null) {
            _articulateNodes = Networks.getLowestArticulationNodes(_network);
        }
        _onePopSize = !Utils.varyPopSizeAcrossBranches();
    }

    private void removeBinaryNodes(Network<NetNodeInfo> net) {
        // Find all binary nodes.
        List<NetNode<NetNodeInfo>> binaryNodes = new LinkedList<NetNode<NetNodeInfo>>();
        for (NetNode<NetNodeInfo> node : net.bfs()) {
            if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                binaryNodes.add(node);
            }
        }
        // Remove them.
        for (NetNode<NetNodeInfo> node : binaryNodes) {
            NetNode<NetNodeInfo> child = node.getChildren().iterator().next();	// Node's only child.
            if(child.getIndeg() != 1){
                continue;
            }
            NetNode<NetNodeInfo> parent = node.getParents().iterator().next();	// Node's only parent.
            double gamma = node.getParentProbability(parent);
            double popSize = child.getParentSupport(node);
            parent.removeChild(node);
            node.removeChild(child);
            parent.adoptChild(child, parent.getData().getHeight() - child.getData().getHeight());
            if(gamma < 0 || gamma > 1) {
                System.err.println("GeneTreeBrSpeciesNetDistribution " + gamma);
            }
            child.setParentProbability(parent, gamma);
            child.setParentSupport(parent, popSize);
        }
    }

    public double calculateGTDistribution(List<UltrametricTree> geneTrees){
        double logL = 0.0;
        for(UltrametricTree gt : geneTrees) {
            logL += calculateGTDistribution(gt, null);
        }
        return logL;
    }

    public double calculateGTDistribution(UltrametricTree geneTree, TreeEmbedding embedding){
        double gtProb = 0;
        String[] gtTaxa = geneTree.getTree().getLeaves();
        List<STITreeCluster<Double>> gtClusters = new ArrayList<>();
        boolean[][] R = processGT(geneTree, gtTaxa, gtClusters);

        Set<String> gtTaxaSet = new HashSet<>();
        Collections.addAll(gtTaxaSet, gtTaxa);

        Map<Tuple<Integer,Integer>, List<Configuration>> edge2ACminus = new HashMap<>();
        Set<IntArray> invalidConfigs = new HashSet<>();
        int netNodeIndex = 0;
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network)){
            List<Configuration> CACs = computeAC(node, gtTaxa, gtTaxaSet, edge2ACminus, invalidConfigs, gtClusters);
            gtProb = computeACMinus(node, netNodeIndex, edge2ACminus, invalidConfigs, gtClusters, CACs, R, embedding);
            if(node.isNetworkNode()){
                netNodeIndex++;
            }
            if(gtProb == -1) {
                throw new IllegalArgumentException("gene tree proability is negative!!!! check!!!!!");
            }
        }
        return gtProb == 0 ? Utils.INVALID_MOVE : Math.log(gtProb);
    }


    private boolean[][] processGT(UltrametricTree geneTree,
                                  String[] gtTaxa,
                                  List<STITreeCluster<Double>> gtClusters){
        Map<TNode, STITreeCluster> map = new HashMap<>();
        for (TNode node : geneTree.getTree().postTraverse()) {
            STITreeCluster cl = new STITreeCluster(gtTaxa);
            double maxTime = 0.0;
            int index = -1;
            if (node.isLeaf()) {
                cl.addLeaf(node.getName());
                ((STINode<Double>)node).setData(maxTime);
            } else {
                for(TNode child : node.getChildren()) {
                    cl = cl.merge(map.get(child));
                }
                maxTime = geneTree.getNodeHeight(node);
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
        // The function is to calculate the _R matrix for the given tree.
        // Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
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

    private List<Configuration> computeAC(NetNode<NetNodeInfo> node,
                                          String[] gtTaxa,
                                          Set<String> gtTaxaSet,
                                          Map<Tuple<Integer,Integer>, List<Configuration>> edge2ACminus,
                                          Set<IntArray> invalidConfigs,
                                          List<STITreeCluster<Double>> gtClusters){

        List<Configuration> CACs =  new ArrayList<>();
        if(node.isLeaf()){
            Configuration config = new Configuration();
            if(_species2alleles == null){
                if(gtTaxaSet.contains(node.getName())){
                    STITreeCluster cl = new STITreeCluster(gtTaxa);
                    cl.addLeaf(node.getName());
                    config.addLineage(gtClusters.indexOf(cl), gtClusters);
                }
            } else{
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
        } else if(node.getOutdeg() == 1){
                Iterator<NetNode<NetNodeInfo>> childNode = node.getChildren().iterator();
                Tuple<Integer,Integer> edge = new Tuple<>(node.getData().getIndex(), childNode.next().getData().getIndex());
                CACs.addAll(edge2ACminus.remove(edge));
        } else{
            boolean isArticulateNode = _articulateNodes.contains(node);
            Iterator<NetNode<NetNodeInfo>> childNode = node.getChildren().iterator();
            Tuple<Integer, Integer> edge1 = new Tuple<>(node.getData().getIndex(), childNode.next().getData().getIndex());
            List<Configuration> AC1 = edge2ACminus.remove(edge1);
            Tuple<Integer, Integer> edge2 = new Tuple<>(node.getData().getIndex(), childNode.next().getData().getIndex());
            List<Configuration> AC2 = edge2ACminus.remove(edge2);
            for(Configuration config1: AC1){
                if(invalidConfigs.contains(new IntArray(config1._netNodeIndex))) continue;
                for(Configuration config2: AC2){
                    if(invalidConfigs.contains(new IntArray(config2._netNodeIndex))) continue;
                    if(config1.isCompatible(config2)){
                        Configuration mergedConfig = new Configuration(config1, config2);
                        if(mergedConfig._totalProb == 0) continue;
                        if(isArticulateNode){
                            mergedConfig.clearNetNodeChoice();
                            int index = CACs.indexOf(mergedConfig);
                            if(index == -1){
                                CACs.add(mergedConfig);
                            } else{
                                Configuration exist = CACs.get(index);
                                exist.addTotalProbability(mergedConfig._totalProb);
                            }
                        } else{
                            CACs.add(mergedConfig);
                        }
                    }
                }
            }
        }
        return CACs;
    }

    private double computeACMinus(NetNode<NetNodeInfo> node,
                                  int netNodeIndex,
                                  Map<Tuple<Integer,Integer>, List<Configuration>> edge2ACminus,
                                  Set<IntArray> invalidConfigs,
                                  List<STITreeCluster<Double>> gtClusters,
                                  List<Configuration> CACs,
                                  boolean[][] R,
                                  TreeEmbedding embedding){
        double gtProb = 0;
        double lowTau = node.getData().getHeight(); // node's height
        double rootPopSize = _network.getRoot().getRootPopSize();
        if(node.isRoot()){
            if(Double.isNaN(rootPopSize)) {
                throw new RuntimeException("Invalid root population size!!!");
            }
            Configuration rootConfig = new Configuration();
            rootConfig.addLineage(gtClusters.size()-1, gtClusters);
            for(Configuration preConfig: CACs){
                Configuration postConfig = new Configuration(preConfig);
                calculateProbability(lowTau, Double.MAX_VALUE, rootPopSize, postConfig, 1, gtClusters, R);
                gtProb += postConfig._totalProb;
            }
            return gtProb;

        } else if(node.isTreeNode()){
            NetNode<NetNodeInfo> parent = node.getParents().iterator().next();
            double highTau = parent.getData().getHeight(); // node's parent's height
            double branchPopSize = _onePopSize ? rootPopSize : node.getParentSupport(parent);
            if(Double.isNaN(branchPopSize)) {
                throw new RuntimeException("Invalid pop size of node " +
                        node.getName() + "<-" + parent.getName() + "!!! \n" + _network.toString());
            }
            List<Configuration> ACminus = new ArrayList<Configuration>();
            Iterator<Configuration> it = CACs.iterator();
            while(it.hasNext()){
                Configuration preConfig = it.next();
                Configuration postConfig = new Configuration(preConfig);
                calculateProbability(lowTau, highTau, branchPopSize, postConfig, 1, gtClusters, R);
                if(postConfig._totalProb == 0){
                    it.remove();
                    invalidConfigs.add(new IntArray(postConfig._netNodeIndex));
                } else{
                    ACminus.add(postConfig);
                }
            }
            if(ACminus.size() == 0 && CACs.size() != 0) return -1;
            Tuple<Integer, Integer> newEdge = new Tuple<>(parent.getData().getIndex(), node.getData().getIndex());
            edge2ACminus.put(newEdge, ACminus);
        } else { // network node
            List<Configuration> newCACs1 = new ArrayList<Configuration>();
            List<Configuration> newCACs2 = new ArrayList<Configuration>();
            int configIndex = 1;
            for(Configuration config: CACs){
                int[] lineageArray = new int[config.getLineageCount()];
                int index = 0;
                for(int lineage: config._lineages){
                    lineageArray[index++] = lineage;
                }

                if(embedding == null) {
                    for (int i = 0; i <= config.getLineageCount(); i++) {
                        for (boolean[] selectedLineages : getSelected(config.getLineageCount(), i)) {
                            Configuration newConfig1 = new Configuration();
                            Configuration newConfig2 = new Configuration();
                            index = 0;
                            for (int lin : lineageArray) {
                                if (selectedLineages[index]) {
                                    newConfig1.addLineage(lin, gtClusters);
                                } else {
                                    newConfig2.addLineage(lin, gtClusters);
                                }
                                index++;
                            }
                            newConfig1.setNetNodeChoice(config._netNodeIndex);
                            newConfig1.setTotalProbability(config._totalProb);
                            newConfig1.addNetNodeChoice(netNodeIndex, configIndex);
                            newCACs1.add(newConfig1);

                            newConfig2.setNetNodeChoice(config._netNodeIndex);
                            newConfig2.setTotalProbability(1);
                            newConfig2.addNetNodeChoice(netNodeIndex, configIndex);
                            newCACs2.add(newConfig2);
                            configIndex++;
                        }
                    }
                } else {
                    Iterator<NetNode<NetNodeInfo>> parentIt = node.getParents().iterator();
                    NetNode<NetNodeInfo> parent1 = parentIt.next();
                    NetNode<NetNodeInfo> parent2 = parentIt.next();
                    Configuration newConfig1 = new Configuration();
                    Configuration newConfig2 = new Configuration();
                    if(embedding.embedding.containsKey(new Tuple<>(parent1, node))) {
                        for (STITreeCluster cl : embedding.embedding.get(new Tuple<>(parent1, node))) {
                            newConfig1.addLineage(gtClusters.indexOf(cl), gtClusters);
                        }
                    }
                    if(embedding.embedding.containsKey(new Tuple<>(parent2, node))) {
                        for (STITreeCluster cl : embedding.embedding.get(new Tuple<>(parent2, node))) {
                            newConfig2.addLineage(gtClusters.indexOf(cl), gtClusters);
                        }
                    }

                    newConfig1.setNetNodeChoice(config._netNodeIndex);
                    newConfig1.setTotalProbability(config._totalProb);
                    newConfig1.addNetNodeChoice(netNodeIndex, configIndex);
                    newCACs1.add(newConfig1);

                    newConfig2.setNetNodeChoice(config._netNodeIndex);
                    newConfig2.setTotalProbability(1);
                    newConfig2.addNetNodeChoice(netNodeIndex, configIndex);
                    newCACs2.add(newConfig2);
                }
            }
            Iterator<NetNode<NetNodeInfo>> parentIt = node.getParents().iterator();
            for(int i=0; i<2; i++){
                List<Configuration> ACminus = new ArrayList<Configuration>();
                List<Configuration> newCACs = (i == 0) ? newCACs1 : newCACs2;

                NetNode<NetNodeInfo> parentNode = parentIt.next();
                double highTau = parentNode.getData().getHeight(); // one parent's height
                double gamma = node.getParentProbability(parentNode); // one parent's gamma
                // one parent's pop size
                double branchPopSize = _onePopSize ? rootPopSize : node.getParentSupport(parentNode);
                if(Double.isNaN(branchPopSize)) {
                    throw new RuntimeException("Invalid node " +
                            node.getName() + "<-" + parentNode.getName() + " population size!!!");
                }

                Iterator<Configuration> it = newCACs.iterator();
                while(it.hasNext()) {
                    Configuration preConfig = it.next();
                    Configuration postConfig = new Configuration(preConfig);
                    calculateProbability(lowTau, highTau, branchPopSize, postConfig, gamma, gtClusters, R);
                    if (postConfig._totalProb == 0) {
                        it.remove();
                        invalidConfigs.add(new IntArray(postConfig._netNodeIndex));
                    } else {
                        ACminus.add(postConfig);
                    }
                }
                if(ACminus.size() == 0 && newCACs.size() != 0 && gamma != 0) return -1;
                Tuple<Integer, Integer> newEdge = new Tuple<>(
                        parentNode.getData().getIndex(), node.getData().getIndex());
                edge2ACminus.put(newEdge, ACminus);

            }
        }
        return gtProb;
    }

    private void calculateProbability(double lowTau,
                                      double highTau,
                                      double popSize, // effective population size added
                                      Configuration config,
                                      double gamma,
                                      List<STITreeCluster<Double>> gtClusters,
                                      boolean R[][]){
        int u = config.getLineageCount();
        if(Double.isNaN(gamma) || gamma == Double.NEGATIVE_INFINITY) System.err.println("INVALID gamma");
        config._totalProb = config._totalProb * Math.pow(gamma, u);
        if(u == 0) return;

        double prob = 1;
        int index = 0;
        for(STITreeCluster<Double> gtCl: gtClusters){
            double coalTime = gtCl.getData();
            if(coalTime >= highTau) break;
            if(coalTime>=lowTau){
                BitSet temp = (BitSet) config._coverage.clone();
                temp.and(gtCl.getCluster());

                if(temp.cardinality() != 0){ // not disjoint
                    if(temp.equals(gtCl.getCluster())){  //contains
                        config.mergeCluster(index, R);
                        prob *= 2.0/popSize * Math.exp(-(coalTime-lowTau) * u*(u-1) / popSize);
                        lowTau = coalTime;
                        u--;
                    } else{
                        config._totalProb = 0;
                        return;
                    }
                }
            }
            index++;
        }
        if(u != 1){
            prob *= Math.exp(-(highTau-lowTau) * u*(u-1) / popSize);
        }
        prob = Math.max(0, prob);
        config._totalProb = config._totalProb * prob;

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

    private class Configuration{
        private HashSet<Integer> _lineages;
        private double _totalProb;
        private BitSet _coverage;
        int[] _netNodeIndex;

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
            if(index == -1) {
                throw new RuntimeException("!!!");
            }
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

    private class IntArray{
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
