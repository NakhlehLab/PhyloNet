package edu.rice.cs.bioinfo.programs.phylonet.algos.integration;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import edu.rice.cs.bioinfo.programs.phmm.src.phylogeny.Felsenstein;
import edu.rice.cs.bioinfo.programs.phmm.src.substitutionModel.GTRSubstitutionModel;
import edu.rice.cs.bioinfo.programs.phmm.src.substitutionModel.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;

import java.io.StringReader;
import java.util.*;

/**
 * Created by yunyu on 6/3/14.
 */
public class CompleteLikelihoodFromSequence {
    Tree _gt;
    double _minBL;
    double _maxBL;
    double _theta;
    SubstitutionModel _model;
    Network _network;
    String[] _networkTaxa;
    Map<String, List<String>> _species2alleles;
    Map<String,char[]> _sequences;
    Map<String,Double> _networkHeightConstraints;

    //int _counter = 0;


    public CompleteLikelihoodFromSequence(Network network, Tree gt, Map<String, List<String>> species2alleles, Map<String,char[]> sequences, SubstitutionModel model, double theta, double minbl, double maxbl){
        _network = network;
        _gt = gt;
        _species2alleles = species2alleles;
        _minBL = minbl;
        _maxBL = maxbl;
        _model = model;
        _theta = theta;
        _sequences = sequences;
        computeNetworkHeightConstraints();
    }

    public double computeLikelihoodWithIntegral(int sampleSize){
        return computeLikelihoodWithIntegral(sampleSize, 1);
    }


    public double computeLikelihoodWithIntegral(int sampleSize, int bins){
        CompleteProbabilityFunction function = new CompleteProbabilityFunction(_network, _gt, _species2alleles, _sequences, _model, _minBL, _maxBL);
        int dimension = function.getNumArguments();
        double[] mins = new double[dimension];
        double[] maxs = new double[dimension];
        computeGTNodeHeightBound(_gt, convertTaxaAssociationMap(_species2alleles), _networkHeightConstraints, _minBL, _maxBL, mins, maxs);

        MultivariateMonteCarloIntegral integration = new MultivariateMonteCarloIntegral(sampleSize, bins);
        double result = integration.integrate(function, dimension, mins, maxs);
        return result;
    }


    private void computeGTNodeHeightBound(Tree tree, Map<String, String> allele2species, Map<String,Double> networkHeightConstraints, double minBL,  double maxBL, double[] mins, double[] maxs){
        Map<TNode,STITreeCluster<Double>> node2cluster = new HashMap<TNode, STITreeCluster<Double>>();
        int index = 0;
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
                    minHeight = Math.max(childCluster.getData(), minHeight);
                    cluster = cluster.merge(childCluster);
                    if(childCluster1==null){
                        childCluster1 = childCluster;
                    }
                    else{
                        childCluster2 = childCluster;
                    }
                }
                minHeight = Math.max(getMinHeightFromSpeciesNetwork(childCluster1, childCluster2, networkHeightConstraints), minHeight);
                mins[index] = minHeight;
                maxs[index++] = minHeight + maxBL;
            }

            cluster.setData(minHeight);
            node2cluster.put(node, cluster);

        }
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


    class CompleteProbabilityFunction implements Func1<double[],Double>{
        Network _network;
        Tree _gt;
        List<MutableTuple<Tree,Double>> _gtList;
        Map<String, List<String>> _species2alleles;
        Map<String,char[]> _sequences;
        int _numArguments;
        double _lowerBound;
        double _upperBound;
        SubstitutionModel _model;
        boolean[][] _R;

        public CompleteProbabilityFunction(Network network, Tree gt, Map<String, List<String>> species2alleles, Map<String,char[]> sequences, SubstitutionModel model, double lowerBound, double upperBound){
            _network = network;
            _gt = gt;
            _gtList = new ArrayList<MutableTuple<Tree,Double>>();
            _gtList.add(new MutableTuple<Tree, Double>(_gt,1.0));
            _species2alleles = species2alleles;
            _numArguments = _gt.getNodeCount() - _gt.getLeafCount();
            _model = model;
            _sequences = sequences;
            _lowerBound = lowerBound;
            _upperBound = upperBound;
            findGTNodesAncestralRelationship();
        }

        public void findGTNodesAncestralRelationship(){
            _R = new boolean[_numArguments][_numArguments];
            Map<TNode, Integer> node2index = new HashMap<TNode, Integer>();
            int index = 0;
            for(TNode node: _gt.postTraverse()){
                for(TNode child: node.getChildren()){
                    if(!child.isLeaf()) {
                        _R[index][node2index.get(child)] = true;
                    }
                }
                if(!node.isLeaf()) {
                    node2index.put(node, index++);
                }
            }
        }

        public Double execute(double[] argument){
            /*
            if(_counter++ % 1000 == 0){
                System.out.println(_counter);
            }
            */

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
            Map<TNode, Double> node2height = new HashMap<TNode, Double>();
            for(TNode node: _gt.postTraverse()){
                double height = 0;
                if(!node.isLeaf()){
                    height = argument[index++];
                }
                for(TNode child: node.getChildren()){
                    child.setParentDistance((height - node2height.get(child))*_theta);
                }
                node2height.put(node, height);
            }

            double result = 1;


            int sequenceLength = 1;
            index = 0;
            do{
                Map<String,Character> omap = new HashMap<String, Character>();
                for(Map.Entry<String,char[]> entry: _sequences.entrySet()){
                    if(index == 0){
                        sequenceLength = entry.getValue().length;
                    }
                    omap.put(entry.getKey(), entry.getValue()[index]);
                }

                ObservationMap observationMap = new ObservationMap(omap);
                Felsenstein fcalc = new Felsenstein(_model);
                result *= fcalc.getLikelihoodtree(_gt, observationMap);

                index++;
            }while(index < sequenceLength);


            for(TNode node: _gt.postTraverse()){
                if(node.isRoot()){
                    node.setParentDistance(node.getParentDistance()/_theta);
                }
            }
            GeneTreeWithBranchLengthProbabilityYF gtp = new GeneTreeWithBranchLengthProbabilityYF();
            double[] probs = new double[1];
            gtp.calculateGTDistribution(_network, _gtList, _species2alleles, probs);
            result *= probs[0];

            //System.out.println(_gt.toString()+": " + result);

            return result;
        }


        public int getNumArguments(){
            return _numArguments;
        }
    }

    public static void main(String[] args){
        GTRSubstitutionModel gtrsm = new GTRSubstitutionModel();
        double[] rates = {1.0, 1.0, 1.0, 1.0, 1.0};
        double[] freqs = {0.25, 0.25, 0.25, 0.25};
        gtrsm.setSubstitutionRates(rates, freqs);
        NewickReader nr = new NewickReader(new StringReader("((B:2.5,C:2.5):3,A:5.5);"));
        STITree<Double> tree = new STITree<Double>(true);
        try {
            nr.readTree(tree);
        }
        catch(Exception e) {
            System.err.println(e);
            e.printStackTrace();
            return;
        }
        //Network net = Networks.string2network("((A:2,(B:1)#H1:1::0.3):1,(C:2,#H1:1::0.7):1);");
        Network net = Networks.string2network("(A:2,(B:1,C:1):1):1;");

        Map<String,char[]> omap = new Hashtable<String,char[]>();
        String s1 = "ATAGT";
        omap.put("A", s1.toCharArray());
        String s2 = "ATTGT";
        omap.put("B", s2.toCharArray());
        String s3 = "TCTAG";
        omap.put("C", s3.toCharArray());

        CompleteLikelihoodFromSequence integrator = new CompleteLikelihoodFromSequence(net, tree, null, omap, gtrsm, 1, 0, 10);
        System.out.println(integrator.computeLikelihoodWithIntegral(1000, 2));
    }
}
