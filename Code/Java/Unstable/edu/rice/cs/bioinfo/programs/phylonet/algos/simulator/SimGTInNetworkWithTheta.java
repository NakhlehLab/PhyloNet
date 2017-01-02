package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.InferTreeWrapper;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 12/26/16
 * Time: 12:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimGTInNetworkWithTheta {
    private boolean _printDetails = false;
    private Random _random;
    private Long _seed = null;

    /**
     * This function is to set seed to control the randomness
     */
    public void setSeed(Long seed){
        _seed = seed;
        if(_seed != null)
            _random = new Random(_seed);
        else
            _random = new Random();
    }

    //parent distance - tau - divergence time in units of expected number of mutations per site
    //parent support - theta - population size parameter in units of population mutation rate per site

    public List<Tree> generateGTs(Network network, Map<String, List<String>> species2alleles, double mu, int numGTs){
        if(_seed == null){
            _random = new Random();
        }

        List<Tree> gts = new ArrayList<Tree>();
        if(species2alleles==null){
            species2alleles = new HashMap<String, List<String>>();
            for(Object leaf: network.getLeaves()){
                String species = ((NetNode)leaf).getName();
                List<String> alleles = new ArrayList<String>();
                alleles.add(species);
                species2alleles.put(species, alleles);
            }
        }
        for(int i=0; i<numGTs; i++){
            if(_printDetails){
                System.out.println("\n\nConstructing gt#" + (i+1) +" ...");
            }

            STITree gt = new STITree();
            STINode root = gt.getRoot();

            Map<Tuple<NetNode,NetNode>, List<TNode>> netEdge2geneLineages = new HashMap<Tuple<NetNode,NetNode>, List<TNode>>();
            for(Object nodeObject: Networks.postTraversal(network)){
                NetNode node = (NetNode)nodeObject;
                if(_printDetails){
                    System.out.println("\nNetNode " + node.getName());
                }
                List<TNode> geneLineages = new ArrayList<TNode>();
                if(node.isLeaf()){
                    for(String allele: species2alleles.get(node.getName())){
                        TNode newNode = root.createChild(allele);
                        newNode.setParentDistance(0.0);
                        geneLineages.add(newNode);
                    }
                }
                else{
                    for(Object childObject: node.getChildren()){
                        NetNode child = (NetNode)childObject;
                        geneLineages.addAll(netEdge2geneLineages.get(new Tuple<NetNode, NetNode>(node, child)));
                    }
                }

                if(_printDetails){
                    System.out.println(geneLineages);
                }



                if(node.isRoot()){
                    randomlyCoalGeneLineages(geneLineages, node, null, mu, root);
                }
                else if(node.isTreeNode()){
                    NetNode parent = (NetNode)node.getParents().iterator().next();
                    randomlyCoalGeneLineages(geneLineages, node, parent, mu, root);
                    netEdge2geneLineages.put(new Tuple<NetNode, NetNode>(parent, node), geneLineages);
                }
                else{
                    NetNode parent1 = (NetNode)node.getParents().iterator().next();
                    double inheritanceProb = node.getParentProbability(parent1);
                    if(inheritanceProb == NetNode.NO_PROBABILITY){
                        throw new RuntimeException("The network has network node that doesn't have inheritance probability.");

                    }
                    List<TNode> geneLineages1 = new ArrayList<TNode>();
                    List<TNode> geneLineages2 = new ArrayList<TNode>();
                    for(TNode gl: geneLineages){
                        double random = _random.nextDouble();
                        if(random < inheritanceProb){
                            geneLineages1.add(gl);
                        }
                        else{
                            geneLineages2.add(gl);
                        }
                    }
                    if(_printDetails){
                        System.out.println("Dividing into :" + geneLineages1);
                        System.out.println("Dividing into :" + geneLineages2);
                    }

                    int index = 0;
                    for(Object parentObject: node.getParents()){
                        NetNode parent = (NetNode)parentObject;
                        if(index == 0){
                            geneLineages = geneLineages1;
                        }
                        else{
                            geneLineages = geneLineages2;
                        }
                        randomlyCoalGeneLineages(geneLineages, node, parent, mu, root);
                        netEdge2geneLineages.put(new Tuple<NetNode, NetNode>(parent, node), geneLineages);
                        index++;
                    }
                }

            }
            Trees.removeBinaryNodes(gt);
            gts.add(gt);
        }

        return gts;
    }

    private void randomlyCoalGeneLineages(List<TNode> geneLineages, NetNode node, NetNode parent, double mu, STINode root) {
        int k = geneLineages.size();
        if(k == 0) return;
        double height = 0.0;
        double length = 0.0;
        double theta = 0.0;
        if(parent != null) {
            length = node.getParentDistance(parent);
            theta = node.getParentSupport(parent);
        } else {
            theta = node.getRootPopSize();
            length = 1e8;
        }
        int numCoal = 0;

        while(true) {
            if(k == 1) {
                if(node.isRoot()) {
                    return ;
                } else {
                    for(TNode lineage : geneLineages) {
                        lineage.setParentDistance(lineage.getParentDistance() + length - height);
                    }
                    break;
                }
            }

            double waittime = randomExp(2.0 / (1.0 * k * (k - 1) * 2.0 * mu / theta));

            if(!node.isRoot() && height + waittime >= length) {
                for(TNode lineage : geneLineages) {
                    lineage.setParentDistance(lineage.getParentDistance() + length - height);
                }
                break;
            }

            numCoal++;
            for(TNode lineage : geneLineages) {
                lineage.setParentDistance(lineage.getParentDistance() + waittime);
            }

            height += waittime;

            int a = _random.nextInt(k);
            int b = _random.nextInt(k - 1);
            if(b >= a) b++;
            else {
                int c = a;
                a = b;
                b = c;
            }

            STINode newNode = root.createChild();
            newNode.setParentDistance(0.0);
            STINode child1 = (STINode) geneLineages.get(a);
            STINode child2 = (STINode) geneLineages.get(b);
            newNode.adoptChild(child1);
            newNode.adoptChild(child2);
            geneLineages.remove(child1);
            geneLineages.remove(child2);
            geneLineages.add(newNode);
            k--;
        }

    }

    private double randomExp(double mean) {
        return -mean * Math.log(_random.nextDouble());
    }

    public static Map<String, Double> comparingTest() {
        double gamma = 1 - 0.35;

        double y = 0.5;
        double x = 0.5;

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork("((((b:1000.0,c:1000.0)I4:" + y + ")I3#H1:0.0::" + (1 - gamma) + ",a:" + (1000 + y) + ")I1:" + x + ",(I3#H1:0.0::" + gamma + ",d:" + (1000 + y) + ")I2:"+ x +")I0;");
        trueNetwork = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");

        SimGTInNetwork simulator = new SimGTInNetwork();
        List<Tree> simulatedTrees = simulator.generateGTs(trueNetwork, null, 100000);
        List<List<MutableTuple<Tree, Double>>> treeList = new ArrayList<>();
        for(Tree tree : simulatedTrees) {
            treeList.add(Arrays.asList(new MutableTuple<Tree, Double>(tree, 1.0)));
        }

        InferTreeWrapper inferTreeWrapper = new InferTreeWrapper();
        List<MutableTuple<Tree, Double>> distinctTrees = new ArrayList<>();
        inferTreeWrapper.summarizeData(treeList, null, distinctTrees);
        Map<String, Double> stat = new TreeMap<>();
        for(MutableTuple<Tree, Double> tuple: distinctTrees) {
            stat.put(Trees.getLexicographicNewickString(tuple.Item1, null), tuple.Item2 );
        }
        for(String s : stat.keySet()) {
            System.out.print(s + '\t');
        }
        System.out.println();
        for(String s : stat.keySet()) {
            System.out.print(stat.get(s).toString() + '\t');
        }
        System.out.println();
        return stat;
    }

    public static void main(String []args) {
        Map<String, Double> stat1 = comparingTest();

        double gamma = 1 - 0.35;

        double y = 0.5;
        double x = 0.5;

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork("((((b:1000.0,c:1000.0)I4:" + y + ")I3#H1:0.0::" + (1 - gamma) + ",a:" + (1000 + y) + ")I1:" + x + ",(I3#H1:0.0::" + gamma + ",d:" + (1000 + y) + ")I2:"+ x +")I0;");
        trueNetwork = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");

        double constTheta = 0.036;
        double mu = 1.0;
        for(NetNode node : Networks.postTraversal(trueNetwork)) {
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, constTheta);
                node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
            }
        }
        trueNetwork.getRoot().setRootPopSize(constTheta);

        SimGTInNetworkWithTheta simulator0 = new SimGTInNetworkWithTheta();
        List<Tree> simulatedTrees = simulator0.generateGTs(trueNetwork, null, mu, 100000);
        List<List<MutableTuple<Tree, Double>>> treeList = new ArrayList<>();
        for(Tree tree : simulatedTrees) {
            treeList.add(Arrays.asList(new MutableTuple<Tree, Double>(tree, 1.0)));
        }

        InferTreeWrapper inferTreeWrapper = new InferTreeWrapper();
        List<MutableTuple<Tree, Double>> distinctTrees = new ArrayList<>();
        inferTreeWrapper.summarizeData(treeList, null, distinctTrees);
        Map<String, Double> stat = new TreeMap<>();
        for(MutableTuple<Tree, Double> tuple: distinctTrees) {
            stat.put(Trees.getLexicographicNewickString(tuple.Item1, null), tuple.Item2 );
        }
        for(String s : stat.keySet()) {
            System.out.print(s + '\t');
        }
        System.out.println();
        for(String s : stat.keySet()) {
            System.out.print(stat.get(s).toString() + '\t');
        }
        System.out.println();

        for(String s : stat1.keySet()) {
            if(!stat.containsKey(s))
                stat.put(s, 0.0);
        }
        for(String s : stat.keySet()) {
            if(!stat1.containsKey(s))
                stat1.put(s, 0.0);
            System.out.println(s + '\t' + stat.get(s) + '\t' + stat1.get(s));
        }
    }
}
