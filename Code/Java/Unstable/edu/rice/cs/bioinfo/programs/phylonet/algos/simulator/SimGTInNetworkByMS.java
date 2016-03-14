package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 3/28/13
 * Time: 2:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimGTInNetworkByMS {


    public List<Tree> generateGTs(Network network, Map<String, List<String>> species2alleles, int numGTs, String MSPath){
        double epsilon = 0.000001;
        Map<String,String> MSName2AlleleName = new Hashtable<String, String>();

        String MSCommand = generateMSCommand(network, species2alleles,numGTs, epsilon, MSName2AlleleName);
        //System.out.println(MSCommand);
        List<Tree> gts = new ArrayList<Tree>();
        try{
            Process proc = Runtime.getRuntime().exec(MSPath+MSCommand,null,null);
            BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            String line;
            while((line=stdout.readLine())!=null){
                //System.out.println(line);
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


    public String generateMSCommand(Network net, Map<String, List<String>> species2alleles, int numGTs, double epsilon, Map<String,String> MSName2AlleleName){
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
        outputCommand = " " + totalAlleles + " " + numGTs + " -T -I " + numPopulations + outputCommand;

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
                //System.out.println(node.getName()+":"+height1);
                NetNode<Integer> parent = node.getParents().iterator().next();
                BitSet bs = new BitSet();
                bs.set(parent.getData());
                bs.set(node.getData());
                edge2population.put(bs, population2);
            }

            else if(node.isNetworkNode()){
                NetNode<Integer> child = node.getChildren().iterator().next();

                //System.out.println(node.getName()+":"+height);
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
                    //System.out.println(node.getName());
                    //System.out.println(MSName2AlleleName);
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

    private List<Tree> summarizeGTs(List<Tree> originalGTs){
        List<Tree> distinctGTs = new ArrayList<Tree>();
        for(Tree gt: originalGTs){
            boolean exist = false;
            for(Tree exgt: distinctGTs){
                if(Trees.haveSameRootedTopology(gt, exgt)){
                    ((STINode<Integer>)exgt.getRoot()).setData(((STINode<Integer>)exgt.getRoot()).getData()+1);
                    exist = true;
                    break;
                }
            }
            if(!exist){
                ((STINode<Integer>)gt.getRoot()).setData(1);
                for(TNode node: gt.getNodes()){
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                distinctGTs.add(gt);
            }
        }
        for(int i=0; i<distinctGTs.size(); i++){
            Tree tr1 = distinctGTs.get(i);
            int count1 = ((STINode<Integer>)tr1.getRoot()).getData();
            for(int j=0; j<i; j++){
                int count2 = ((STINode<Integer>)distinctGTs.get(j).getRoot()).getData();
                if(count2 < count1){
                    distinctGTs.remove(i);
                    distinctGTs.add(j,tr1);
                    break;
                }
            }
        }
        return distinctGTs;
    }
}
