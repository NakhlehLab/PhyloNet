package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.InferTreeWrapper;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 4/19/17
 * Time: 2:33 PM
 * To change this template use File | Settings | File Templates.
 */

// TODO: Not finish yet!!!

public class SimGTInNetworkByMSWithRecombination {
    private static String _msdir = "/Users/zhujiafan/Documents/Luay/msdir/ms";

    /**
     * This is the main function for using ms to generate gene trees given a species network
     *
     * @param network             the species network (need to be ultrametric) that the gene trees are generated from
     * @param species2alleles     mapping from species to alleles sampled from it
     * @param numGTs              the number of gene trees need to be generated
     * @param recombinationRate   \rho=4 N_0 r
     * @param recombinationNSites the number of sites between which recombination can occur
     * @param MSPath              local ms path
     *
     * @return   the list of simulated gene trees
     */
    public List<Tree> generateGTs(Network network, Map<String, List<String>> species2alleles, int numGTs, double recombinationRate, int recombinationNSites, String MSPath){
        double epsilon = 0.000001;
        Map<String,String> MSName2AlleleName = new Hashtable<String, String>();

        String MSCommand = generateMSCommand(network, species2alleles,numGTs, epsilon, MSName2AlleleName, recombinationRate, recombinationNSites);
        System.out.println(MSCommand);
        List<Tree> gts = new ArrayList<Tree>();
        try{
            Process proc = Runtime.getRuntime().exec(MSPath+MSCommand,null,null);
            BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            String line;
            while((line=stdout.readLine())!=null){
                line = line.trim();
                if(line.startsWith("//")){
                    line=stdout.readLine();
                    NewickReader nr = new NewickReader(new StringReader(line));
                    Tree newtr = nr.readTree();

                    for(TNode node: newtr.postTraverse()){
                        if(node.isLeaf()){
                            ((STINode)node).setName("MSTemp"+ MSName2AlleleName.get(node.getName()));
                        }
                        node.setParentDistance(node.getParentDistance()*2);
                    }
                    for(TNode node: newtr.postTraverse()){
                        if(node.isLeaf()){
                            ((STINode)node).setName(node.getName().substring(6));
                        }
                    }
                    gts.add(newtr);
                }
            }
            proc.waitFor();
            stdout.close();
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return gts;
    }



    /**
     * This function is to generate the ms command
     *
     * @param net                 the species network (need to be ultrametric) that the gene trees are generated from
     * @param species2alleles     mapping from species to alleles sampled from it
     * @param numGTs              the number of gene trees need to be generated
     * @param epsilon             the threshold for checking ultrametricy of a network
     * @param MSName2AlleleName   mapping from the taxa of gene trees generated by ms to original names
     * @param recombinationRate   \rho=4 N_0 r
     * @param recombinationNSites the number of sites between which recombination can occur
     *
     * @return   the ms command
     */
    public String generateMSCommand(Network net, Map<String, List<String>> species2alleles, int numGTs, double epsilon, Map<String,String> MSName2AlleleName, double recombinationRate, int recombinationNSites){
        Map<String, Integer> speciesName2MSName = new Hashtable<String, Integer>();
        List<Integer> numAllelesList = new ArrayList<Integer>();
        List<Tuple<NetNode<Integer>,Double>> nodeTimeList = new ArrayList<Tuple<NetNode<Integer>,Double>>();
        processNetwork(net, epsilon, species2alleles, speciesName2MSName, numAllelesList, MSName2AlleleName, nodeTimeList);
        int numPopulations = speciesName2MSName.size();
        int totalAlleles = 0;
        String outputCommand = "";
        for(int numAllele: numAllelesList){
            outputCommand += " "+numAllele;
            totalAlleles += numAllele;
        }
        outputCommand = " " + totalAlleles + " " + numGTs + " -T -r " + recombinationRate + " " + recombinationNSites + " -I " + numPopulations + outputCommand;

        Map<BitSet, Integer> edge2population = new HashMap<BitSet, Integer>();

        for(Tuple<NetNode<Integer>,Double> bundle: nodeTimeList){
            NetNode<Integer> node = bundle.Item1;
            double height = bundle.Item2;
            if(node.isLeaf()){
                NetNode<Integer> parent = node.getParents().iterator().next();
                BitSet bs = new BitSet();
                bs.set(parent.getData());
                bs.set(node.getData());
                edge2population.put(bs, speciesName2MSName.get(node.getName()));
            }
            else if(node.isTreeNode()){
                Iterator<NetNode<Integer>> children = node.getChildren().iterator();
                NetNode<Integer> child1 = children.next();
                NetNode<Integer> child2 = children.next();
                BitSet edge1 = new BitSet();
                edge1.set(child1.getData());
                edge1.set(node.getData());
                int population1 = edge2population.get(edge1);
                BitSet edge2 = new BitSet();
                edge2.set(node.getData());
                edge2.set(child2.getData());
                int population2 = edge2population.get(edge2);
                outputCommand += " -ej " + height + " " + population1 + " " + population2;

                if(node.isRoot()){
                    break;
                }
                NetNode<Integer> parent = node.getParents().iterator().next();
                BitSet bs = new BitSet();
                bs.set(parent.getData());
                bs.set(node.getData());
                edge2population.put(bs, population2);
            }
            else if(node.isNetworkNode()){
                NetNode<Integer> child = node.getChildren().iterator().next();
                BitSet childEdge = new BitSet();
                childEdge.set(child.getData());
                childEdge.set(node.getData());
                int population = edge2population.get(childEdge);
                Iterator<NetNode<Integer>> parents = node.getParents().iterator();
                NetNode<Integer> parent1 = parents.next();
                double hybridProb = node.getParentProbability(parent1);
                outputCommand += " -es " + height + " " + population + " " + hybridProb;
                BitSet parentEdge1 = new BitSet();
                parentEdge1.set(parent1.getData());
                parentEdge1.set(node.getData());
                edge2population.put(parentEdge1, population);
                NetNode<Integer> parent2 = parents.next();
                BitSet parentEdge2 = new BitSet();
                parentEdge2.set(parent2.getData());
                parentEdge2.set(node.getData());
                numPopulations++;
                edge2population.put(parentEdge2, numPopulations);
            }
        }

        return outputCommand;

    }



    /**
     * This function is to gather information of the species network for generating ms command
     *
     * @param net                   the species network (need to be ultrametric) that the gene trees are generated from
     * @param epsilon               the threshold for checking ultrametricy of a network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param speciesName2MSName    mapping from the population id to the taxa of gene trees generated by ms
     * @param MSName2AlleleName     mapping from the taxa of gene trees generated by ms to original names
     * @param numAllelesList        the number of alleles sampled from every corresponding species
     * @param nodeTimeList          the heights of every internal node of the species network
     */
    private void processNetwork(Network<Integer> net, double epsilon, Map<String, List<String>> species2alleles, Map<String, Integer> speciesName2MSName, List<Integer> numAllelesList, Map<String,String> MSName2AlleleName, List<Tuple<NetNode<Integer>,Double>> nodeTimeList){
        int nodeIndex = 0;
        int leafIndex = 1;
        int populationIndex = 1;
        if(species2alleles == null){
            species2alleles = new Hashtable<String, List<String>>();
            for(NetNode<Integer> node: net.dfs()){
                node.setData(nodeIndex++);
                if(node.isLeaf()){
                    List<String> alleles = new ArrayList<String>();
                    alleles.add(node.getName());
                    species2alleles.put(node.getName(), alleles);
                    speciesName2MSName.put(node.getName(), leafIndex++);
                    numAllelesList.add(1);
                    MSName2AlleleName.put((populationIndex++)+"", node.getName());
                }
            }
        }
        else{
            for(NetNode<Integer> node: net.dfs()){
                node.setData(nodeIndex++);
                if(node.isLeaf()){
                    if(!species2alleles.containsKey(node.getName())){
                        throw new RuntimeException("The species " + node.getName() + " doesn't specify the number of alleles sampled.");
                    }
                    speciesName2MSName.put(node.getName(), leafIndex++);
                    List<String> alleles = species2alleles.get(node.getName());
                    int numAlleles = alleles.size();
                    numAllelesList.add(numAlleles);
                    for(String allele: alleles){
                        MSName2AlleleName.put((populationIndex++)+"", allele);
                    }
                }
            }
        }
        if(species2alleles.size()!=leafIndex-1){
            throw new RuntimeException("The mapping for alleles doesn't match the network.");
        }
        Map<Integer, Double> node2height = new HashMap<Integer, Double>();
        for(NetNode<Integer> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(net)){
            double height;
            if(node.isLeaf()){
                node2height.put(node.getData(), 0.0);
                height = 0;
            }
            else if(node.isTreeNode()){
                Iterator<NetNode<Integer>> children = node.getChildren().iterator();
                NetNode<Integer> child1 = children.next();
                NetNode<Integer> child2 = children.next();
                Double height1 = node2height.get(child1.getData()) + child1.getParentDistance(node)*0.5;
                Double height2 = node2height.get(child2.getData()) + child2.getParentDistance(node)*0.5;
                if(Math.abs(height1-height2)>epsilon){
                    throw new RuntimeException("The network is not ultrametric! (" + node.getName() + ")");
                }
                height = Math.max(height1,height2);
                node2height.put(node.getData(), height);
            }
            else{
                NetNode<Integer> child = node.getChildren().iterator().next();
                height = node2height.get(child.getData()) + child.getParentDistance(node)*0.5;
                node2height.put(node.getData(),height);
            }
            Tuple<NetNode<Integer>,Double> bundle = new Tuple<NetNode<Integer>,Double>(node, height);
            int i = 0;
            for(; i<nodeTimeList.size(); i++){
                if(nodeTimeList.get(i).Item2>height){
                    break;
                }
            }
            nodeTimeList.add(i, bundle);
        }
    }

    public static void main(String []args) {
        double pi0 = 0.5;
        double pi1 = 1.0 - pi0;
        double u = 1.0 / (2.0 * pi0);
        double v = 1.0 / (2.0 * pi1);
        double mu = 2.0 * u * v / (u + v);

        double gamma = 1 - 0.3;

        double y = 1.5;
        double x = 0.5;

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork("((((b:1.0,c:1.0)I4:" + y + ")I3#H1:0.1::" + (1 - gamma) + ",a:" + (1.1 + y) + ")I1:" + x + ",(I3#H1:0.1::" + gamma + ",d:" + (1.1 + y) + ")I2:"+ x +")I0;");

        SimGTInNetworkWithTheta sim1 = new SimGTInNetworkWithTheta();
        SimGTInNetworkByMS simMS = new SimGTInNetworkByMS();
        SimGTInNetworkByMSWithRecombination simMSRe = new SimGTInNetworkByMSWithRecombination();

        int numGTs = 10;

        GeneTreeProbability geneTreeProbability = new GeneTreeProbability();
        List<Double> gtprob = new ArrayList<>();
        List<Tree> gts = new ArrayList<>();
        gts.add(Trees.readTree("(((b, c), a), d);"));
        gts.add(Trees.readTree("(((b, c), d), a);"));
        gts.add(Trees.readTree("((a, b), (c, d));"));
        gts.add(Trees.readTree("((a, c), (b, d));"));
        gts.add(Trees.readTree("(((a, b), c), d);"));
        gts.add(Trees.readTree("(((a, c), b), d);"));
        gts.add(Trees.readTree("(a, (b, (c, d)));"));
        gts.add(Trees.readTree("(((b, d), c), a);"));
        gts.add(Trees.readTree("((a, d), (b, c));"));
        gts.add(Trees.readTree("(((a, b), d), c);"));
        gts.add(Trees.readTree("(b, (a, (c, d)));"));
        gts.add(Trees.readTree("(((a, d), b), c);"));
        gts.add(Trees.readTree("(((b, d), a), c);"));
        gts.add(Trees.readTree("(((a, c), d), b);"));
        gts.add(Trees.readTree("(((a, d), c), b);"));

        for(Tree gt : gts) {
            Trees.convertToLexicographicTree(gt);
        }


        gtprob = geneTreeProbability.calculateGTDistribution((Network)trueNetwork, gts, null, false);
        //System.out.println(gts);
        for(double prob : gtprob){
            System.out.print('\t');
            System.out.print(prob);
        }
        System.out.println("");

        List<Tree> treesFromMS = simMS.generateGTs(trueNetwork, null, numGTs, _msdir);
        List<List<MutableTuple>> treesFromMS0 = new ArrayList<List<MutableTuple>>();
        for (Tree tr : treesFromMS) {
            //for (Tree tr : simulator.generateGTs(network, species2alleles, numGTs)) {
            treesFromMS0.add(Arrays.asList(new MutableTuple(tr, 1.0)));
        }

        List<Tree> treesFromMSRe = simMSRe.generateGTs(trueNetwork, null, numGTs, 0.0075, 500000, _msdir);
        List<List<MutableTuple>> treesFromMSRe0 = new ArrayList<List<MutableTuple>>();
        for (Tree tr : treesFromMS) {
            //for (Tree tr : simulator.generateGTs(network, species2alleles, numGTs)) {
            treesFromMSRe0.add(Arrays.asList(new MutableTuple(tr, 1.0)));
        }

        double constTheta = 0.036;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, constTheta);
                node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
            }
        }
        trueNetwork.getRoot().setRootPopSize(constTheta);
        List<Tree> trees = sim1.generateGTs(trueNetwork, null, mu, numGTs);
        List<List<MutableTuple>> trees0 = new ArrayList<List<MutableTuple>>();
        for (Tree tr : trees) {
            //for (Tree tr : simulator.generateGTs(network, species2alleles, numGTs)) {
            trees0.add(Arrays.asList(new MutableTuple(tr, 1.0)));
        }

        InferTreeWrapper wrapper = new InferTreeWrapper();
        List<List<MutableTuple<Tree, Double>>> dataForStartingNetwork = new ArrayList();
        wrapper.summarizeData(treesFromMS0, null, dataForStartingNetwork);
        Map<String, Double> treesFromMS1 = new HashMap<>();
        for(Object object :  dataForStartingNetwork) {
            MutableTuple<Tree, Double> tuple = (MutableTuple<Tree, Double>) object;
            //MutableTuple<Tree, Double> tuple = (MutableTuple<Tree, Double>) (list.get(0));
            Trees.convertToLexicographicTree(tuple.Item1);
            treesFromMS1.put(tuple.Item1.toString(), tuple.Item2);
        }
        for(Tree gt : gts){
            System.out.print('\t');
            System.out.print(treesFromMS1.get(gt.toString()));
        }
        System.out.println("");

        dataForStartingNetwork.clear();
        wrapper.summarizeData(trees0, null, dataForStartingNetwork);
        Map<String, Double> trees1 = new HashMap<>();
        for(Object object :  dataForStartingNetwork) {
            MutableTuple<Tree, Double> tuple = (MutableTuple<Tree, Double>) object;
            //MutableTuple<Tree, Double> tuple = (MutableTuple<Tree, Double>) (list.get(0));
            Trees.convertToLexicographicTree(tuple.Item1);
            trees1.put(tuple.Item1.toString(), tuple.Item2);
        }
        for(Tree gt : gts){
            System.out.print('\t');
            System.out.print(trees1.get(gt.toString()));
        }
        System.out.println("");

        dataForStartingNetwork.clear();
        wrapper.summarizeData(treesFromMSRe0, null, dataForStartingNetwork);
        Map<String, Double> treesFromMSRe1 = new HashMap<>();
        for(Object object :  dataForStartingNetwork) {
            MutableTuple<Tree, Double> tuple = (MutableTuple<Tree, Double>) object;
            //MutableTuple<Tree, Double> tuple = (MutableTuple<Tree, Double>) (list.get(0));
            Trees.convertToLexicographicTree(tuple.Item1);
            treesFromMSRe1.put(tuple.Item1.toString(), tuple.Item2);
        }
        for(Tree gt : gts){
            System.out.print('\t');
            System.out.print(treesFromMSRe1.get(gt.toString()));
        }
        System.out.println("");

    }
}
