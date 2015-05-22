package edu.rice.cs.bioinfo.programs.phylonet.algos.integration;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.*;

/**
 * Created by yunyu on 6/3/14.
 */
public class GTBranchLengthsIntegrationForSpeciesPhylogeny {
    Tree _gt;
    Network _network;
    String[] _networkTaxa;
    Map<String, List<String>> _species2alleles;
    Map<String,Double> _networkHeightConstraints;
    Set<NetNode> _articulateNodes;

    //int _batchSize = 1;


    public GTBranchLengthsIntegrationForSpeciesPhylogeny(Network network, Tree gt, Map<String, List<String>> species2alleles){
        _network = network;
        _gt = gt;
        _species2alleles = species2alleles;
        computeNetworkHeightConstraints();
    }

    public void setArticulateNodes(Set<NetNode> articulateNodes){
        _articulateNodes = articulateNodes;
    }

    public double computeLikelihoodWithIntegral(int sampleSize){
        return computeLikelihoodWithIntegral(sampleSize, 1);
    }


    public double computeLikelihoodWithIntegral(int sampleSize, int bins){
        double result = 0;

        nameAllInternalNodes(_gt);
        Map<String, Double> gtNodeName2height = getNodeHeight(_gt);

        for(Tree resolvedGT: Trees.getAllBinaryResolution(_gt)){
            List<TNode> nodesToIntegrate = new ArrayList<TNode>();
            List<Double> minList = new ArrayList<Double>();
            List<Double> maxList = new ArrayList<Double>();
            getNodesToIntegrate(resolvedGT, convertTaxaAssociationMap(_species2alleles), _networkHeightConstraints, gtNodeName2height, nodesToIntegrate, minList, maxList);
            //integrate(resolvedGT, nodesToIntegrate, minList, maxList, sampleSize, bins);
            result += integrate(resolvedGT, nodesToIntegrate, minList, maxList, sampleSize, bins);
        }
        return result;
    }

    private void computeNetworkHeightConstraints(){
        List<String> leafList = new ArrayList<String>();
        for(Object o: _network.getLeaves()){
            leafList.add(((NetNode)o).getName());
        }
        _networkTaxa = leafList.toArray(new String[leafList.size()]);
        Map<NetNode, STITreeCluster<Double>> node2cluster = new HashMap<NetNode, STITreeCluster<Double>>();
        _networkHeightConstraints = new HashMap<String, Double>();
        for(Object o: Networks.postTraversal(_network)){
            NetNode node = (NetNode)o;
            STITreeCluster<Double> cluster = new STITreeCluster<Double>(_networkTaxa);
            double clusterHeight = -1;
            if(node.isLeaf()){
                cluster.addLeaf(node.getName());
                clusterHeight = 0;

            }
            else{
                BitSet commonCluster = null;
                STITreeCluster childCluster1 = null;
                STITreeCluster childCluster2 = null;
                for(Object co: node.getChildren()){
                    NetNode childNode = (NetNode)co;
                    STITreeCluster<Double> childCluster = node2cluster.get(childNode).duplicate();
                    if(childCluster1==null){
                        childCluster1 = childCluster;
                    }
                    else{
                        childCluster2 = childCluster;
                    }
                    cluster = cluster.merge(childCluster);
                    if(clusterHeight==-1){
                        clusterHeight = childCluster.getData() + childNode.getParentDistance(node);
                    }
                    if(node.isTreeNode()) {
                        if (commonCluster == null) {
                            commonCluster = (BitSet)childCluster.getCluster().clone();
                        } else {
                            commonCluster.and(childCluster.getCluster());
                        }
                    }
                }
                if(node.isTreeNode()) {
                    childCluster1.getCluster().andNot(commonCluster);
                    childCluster2.getCluster().andNot(commonCluster);
                    for(String t1: childCluster1.getClusterLeaves()){
                        for(String t2: childCluster2.getClusterLeaves()){
                            String pair = t1.compareTo(t2)>0?t1+"|"+t2:t2+"|"+t1;
                            Double preDist = _networkHeightConstraints.get(pair);
                            if(preDist==null || preDist>clusterHeight){
                                _networkHeightConstraints.put(pair, clusterHeight);

                            }
                        }
                    }
                }
            }
            cluster.setData(clusterHeight);
            node2cluster.put(node, cluster);
        }

    }



    private Map<String, Double> getNodeHeight(Tree tree){
        Map<String, Double> nodeName2height = new HashMap<String, Double>();
        for(TNode node: tree.postTraverse()){
            double height = 0;
            if(!node.isLeaf()) {
                TNode child = node.getChildren().iterator().next();
                height = child.getParentDistance() + nodeName2height.get(child.getName());
            }
            nodeName2height.put(node.getName(),height);
        }
        return nodeName2height;
    }

    private Map<String,String> convertTaxaAssociationMap(Map<String, List<String>> species2alleles){
        if(species2alleles==null){
            return null;
        }
        Map<String, String> allele2species = new HashMap<String, String>();
        for(Map.Entry<String, List<String>> entry: species2alleles.entrySet()){
            for(String allele: entry.getValue()){
                allele2species.put(entry.getKey(), allele);
            }
        }
        return allele2species;
    }

    private void getNodesToIntegrate(Tree tree, Map<String, String> allele2species, Map<String,Double> networkHeightConstraints, Map<String, Double> gtNodeName2height, List<TNode> nodesToIntegrate, List<Double> mins, List<Double> maxs){
        Map<TNode,STITreeCluster<Double>> node2cluster = new HashMap<TNode, STITreeCluster<Double>>();
        for(TNode node: tree.postTraverse()){
            STITreeCluster<Double> cluster = new STITreeCluster<Double>(_networkTaxa);
            double minHeight = 0;
            if(node.isLeaf()){
                if(allele2species==null) {
                    cluster.addLeaf(node.getName());
                }else{
                    cluster.addLeaf(allele2species.get(node.getName()));
                }
            }
            else{
                STITreeCluster<Double> childCluster1 = null;
                STITreeCluster<Double> childCluster2 = null;
                for(TNode child: node.getChildren()){
                    STITreeCluster<Double> childCluster = node2cluster.get(child);
                    cluster = cluster.merge(childCluster);
                    if(childCluster1==null){
                        childCluster1 = childCluster;
                    }
                    else{
                        childCluster2 = childCluster;
                    }
                }
                if(node.getName().equals("")){
                    nodesToIntegrate.add(node);
                    for(TNode child: node.getChildren()){
                        minHeight = Math.max(node2cluster.get(child).getData(), minHeight);
                    }
                    minHeight = Math.max(getMinHeightFromSpeciesNetwork(childCluster1, childCluster2, networkHeightConstraints), minHeight);
                }
                else{
                    minHeight = gtNodeName2height.get(node.getName());
                }
            }

            cluster.setData(minHeight);
            node2cluster.put(node, cluster);

            if(!node.isRoot() && node.getParentDistance() == TNode.NO_DISTANCE){
                Double nodeHeight = gtNodeName2height.get(node.getName());
                Double parentHeight = gtNodeName2height.get(node.getParent().getName());
                if(nodeHeight!=null && parentHeight!=null){
                    node.setParentDistance(parentHeight - nodeHeight);
                }
            }
        }

        for(TNode node: nodesToIntegrate){
            mins.add(node2cluster.get(node).getData());
        }

        for(TNode node: nodesToIntegrate){
            Double max = null;
            do{
                node = node.getParent();
                max = gtNodeName2height.get(node.getName());
            }while(max==null);
            maxs.add(max);
        }
    }

    private double getMinHeightFromSpeciesNetwork(STITreeCluster cl1, STITreeCluster cl2, Map<String,Double> networkHeightConstraints){
        double minHeight = 0;
        for(String t1: cl1.getClusterLeaves()){
            for(String t2: cl2.getClusterLeaves()){
                String pair = t1.compareTo(t2)>0?t1+"|"+t2:t2+"|"+t1;
                minHeight = Math.max(minHeight, networkHeightConstraints.get(pair));
            }
        }
        return minHeight;
    }


    private void nameAllInternalNodes(Tree tree){
        int index = 1;
        for(TNode node: tree.postTraverse()){
            if(node.getName().equals("")){
                while(tree.getNode("i"+index)!=null){
                    index++;
                }
                ((STINode)node).setName("i"+index);
            }

        }

    }





    private double integrate(Tree resolvedGT, List<TNode> nodes, List<Double> minList, List<Double> maxList, int sampleSize, int bins){
        double[] mins = new double[minList.size()];
        int index = 0;
        for(double value: minList){
            mins[index++] = value;
        }
        double[] maxs = new double[maxList.size()];
        index = 0;
        for(double value: maxList){
            maxs[index++] = value;
        }

        MultivariateMonteCarloIntegral integration = new MultivariateMonteCarloIntegral(sampleSize, bins);
        GTProbInNetworkFunction function = new GTProbInNetworkFunction(_network, resolvedGT, _species2alleles, nodes, mins, maxs);
        double result = integration.integrate(function, nodes.size(), mins, maxs);
        return result;
    }


    class GTProbInNetworkFunction implements Func1<double[], Double> {
        Network _network;
        Tree _gt;
        Map<String, List<String>> _species2alleles;
        int _numArguments;
        double[] _lowerBound;
        double[] _upperBound;
        List<TNode> _nodes; //keep it in post order
        boolean[][] _R;
        Map<TNode, Double> _node2height;
        List<Tree> _gtList;


        public GTProbInNetworkFunction(Network network, Tree gt, Map<String, List<String>> species2alleles, List<TNode> nodes, double[] lowerBound, double[] upperBound){
            _network = network;
            _gt = gt;
            _species2alleles = species2alleles;
            _numArguments = nodes.size();
            _lowerBound = lowerBound;
            _upperBound = upperBound;
            _nodes = nodes;
            _gtList = new ArrayList<Tree>();
            _gtList.add(_gt);
            findGTNodesAncestralRelationship();
        }

        private void findGTNodesAncestralRelationship(){
            _R = new boolean[_numArguments][_numArguments];
            _node2height = new HashMap<TNode, Double>();
            for(TNode node: _gt.postTraverse()){
                double height = 0;
                if(!node.isLeaf()){
                    for(TNode child: node.getChildren()){
                        if(_node2height.containsKey(child) && child.getParentDistance()!=TNode.NO_DISTANCE){
                            height = child.getParentDistance() + _node2height.get(child);
                        }
                    }
                }
                int nodeID = _nodes.indexOf(node);
                int parentID = _nodes.indexOf(node.getParent());
                if(nodeID!=-1 && parentID!=-1){
                    _R[parentID][nodeID] = true;
                }
                _node2height.put(node, height);
            }
        }



        public Double execute(double[] argument){

            for(int i=0; i<_numArguments; i++){
                for(int j=0; j<_numArguments; j++){
                    if(_R[i][j]){
                        if(argument[i] <= argument[j]){
                            return 0.0;
                        }
                    }
                }
            }
            int index = 0;
            for(TNode node: _nodes){
                for(TNode child: node.getChildren()){
                    child.setParentDistance(argument[index] - _node2height.get(child));
                }
                TNode parent = node.getParent();
                if(!_nodes.contains(parent)) {
                    node.setParentDistance(_node2height.get(parent) - argument[index]);
                }
                _node2height.put(node, argument[index]);
                index++;
            }

            if(!Trees.isUltrametric(_gt, 0.001)){
                throw new RuntimeException("Not Ultrametric");
            }

            GeneTreeWithBranchLengthProbabilityYF gtp = new GeneTreeWithBranchLengthProbabilityYF(_network, _gtList, _species2alleles);
            if(_articulateNodes!=null){
                gtp.setTotalNodes(_articulateNodes);
            }
            double[] probs = new double[1];
            gtp.calculateGTDistribution(probs);
            //System.out.println(_gt + ": " + probs[0]);
            //if(probs[0]!=0)System.out.println(argument[0] + ": " + probs[0]);
            return probs[0];
        }

    }

    public static void main(String[] args){
        //Network net = Networks.string2network("((A:2,(B:1)#H1:1::0.3):1,(C:2,#H1:1::0.7):1);");


        try{
            Network net = Networks.readNetwork("(((((D:1)#H1:1::0.3,C:2):1,(B:2)#H2:1::0.7):1,(A:3,#H2:1::0.3):1):1,(E:2,#H1:1::0.7):3);");
            String gtString = "(A:4,B:4,C:4);";
            NewickReader nr = new NewickReader(new StringReader(gtString));
            Tree gt = nr.readTree();
            double[] minHeight = new double[2];
            Arrays.fill(minHeight, 0);
            double[] maxHeight = new double[2];
            Arrays.fill(maxHeight, 10);
            GTBranchLengthsIntegrationForSpeciesPhylogeny integrator = new GTBranchLengthsIntegrationForSpeciesPhylogeny(net, gt, null);
            System.out.println(integrator.computeLikelihoodWithIntegral(50000));

        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
    }
}
