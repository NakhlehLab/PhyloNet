package edu.rice.cs.bioinfo.programs.phylonet.structs.network.util;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import org.apache.commons.math3.util.ArithmeticUtils;
import java.util.*;

/**
 * Created by yunyu on 9/16/14.
 */
public class NetworkLoad {
    
    /*
    public static  int computeMaxLoad(Network network, Map<String, List<String>> species2alleles){
        int maxLoad = 0;
        Set<NetNode> lowestArticulateNodes = Networks.getLowestArticulationNodes(network);
        Map<NetNode, Integer> node2load = new HashMap<>();
        Map<NetNode, Set<String>> node2taxaUnder = new HashMap<>();
        for(NetNode node: Networks.postTraversal(network)){
            if(node.isLeaf()){
                Set<String> taxaUnder = new HashSet<>();
                int load = 1;
                if(species2alleles == null){
                    taxaUnder.add(node.getName());
                }
                else{
                    List<String> alleles = species2alleles.get(node.getName());
                    taxaUnder.addAll(alleles);
                    load = alleles.size();

                }
                node2taxaUnder.put(node, taxaUnder);
                node2load.put(node, load);
            }
            else if(node.isNetworkNode()){
                NetNode childNode = node.getChildren().iterator().next();
                Set<String> taxaUnder = node2taxaUnder.get(childNode);
                node2taxaUnder.put(node, taxaUnder);
                int load = node2load.get(childNode) + taxaUnder.size();
                node.setName(node.getName()+"_"+load);
                node2load.put(node, load);
            }
            else{
                boolean isArticulate = lowestArticulateNodes.contains(node);
                Set<String> taxaUnder = new HashSet<>();
                int load = 0;
                for(NetNode childNode: node.getChildren()){
                    taxaUnder.addAll(node2taxaUnder.get(childNode));
                    if(!isArticulate){
                        load += node2load.get(childNode);
                    }
                }
                if(isArticulate){
                    load = taxaUnder.size();
                }
                node2taxaUnder.put(node, taxaUnder);
                node2load.put(node, load);
                //System.out.println(node.getName() + ": " + load);
                node.setName(node.getName()+"_"+load);
                maxLoad = Math.max(maxLoad, load);
            }
        }
        return maxLoad;
    }



    public static  double computeMaxLoad(Network network, Map<String, List<String>> species2alleles){
        double maxLoad = 0;
        Set<NetNode> lowestArticulateNodes = Networks.getLowestArticulationNodes(network);
        Map<NetNode, Double> node2load = new HashMap<>();
        Map<NetNode, Set<String>> node2taxaUnder = new HashMap<>();
        for(NetNode node: Networks.postTraversal(network)){
            if(node.isLeaf()){
                Set<String> taxaUnder = new HashSet<>();
                double load = 1;
                if(species2alleles == null){
                    taxaUnder.add(node.getName());
                }
                else{
                    List<String> alleles = species2alleles.get(node.getName());
                    taxaUnder.addAll(alleles);
                    load = Math.pow(2, alleles.size()/2);

                }
                node2taxaUnder.put(node, taxaUnder);
                node2load.put(node, load);
            }
            else if(node.isNetworkNode()){
                NetNode childNode = node.getChildren().iterator().next();
                Set<String> taxaUnder = node2taxaUnder.get(childNode);
                node2taxaUnder.put(node, taxaUnder);
                double load = Math.pow(2, taxaUnder.size()+1);
                node.setName(node.getName()+"_"+load);
                node2load.put(node, load);
            }
            else{
                boolean isArticulate = lowestArticulateNodes.contains(node);
                Set<String> taxaUnder = new HashSet<>();
                double load = 1;
                for(NetNode childNode: node.getChildren()){
                    taxaUnder.addAll(node2taxaUnder.get(childNode));
                    if(!isArticulate){
                        double childLoad = node2load.get(childNode);
                        load *= childLoad==1&&childNode.isLeaf()? 1: childLoad*1.5;
                    }
                }
                if(isArticulate){
                    load = Math.pow(2, taxaUnder.size()/2);
                }
                node2taxaUnder.put(node, taxaUnder);
                node2load.put(node, load);
                //System.out.println(node.getName() + ": " + load);
                node.setName(node.getName()+"_"+load);
                maxLoad = Math.max(maxLoad, load);
            }
        }
        return maxLoad;
    }




    public static  double computeMaxLoad(Network network, Map<String, List<String>> species2alleles){
        double maxLoad = 0;
        Set<NetNode> lowestArticulateNodes = Networks.getLowestArticulationNodes(network);
        Map<NetNode, Double> node2load = new HashMap<>();
        Map<NetNode, Set<String>> node2taxaUnder = new HashMap<>();
        for(NetNode node: Networks.postTraversal(network)){
            if(node.isLeaf()){
                Set<String> taxaUnder = new HashSet<>();
                double load = 1;
                if(species2alleles == null){
                    taxaUnder.add(node.getName());
                }
                else{
                    List<String> alleles = species2alleles.get(node.getName());
                    taxaUnder.addAll(alleles);
                    load = Math.pow(2, alleles.size()/2);

                }
                node2taxaUnder.put(node, taxaUnder);
                node2load.put(node, load);
            }
            else if(node.isNetworkNode()){
                NetNode childNode = node.getChildren().iterator().next();
                Set<String> taxaUnder = node2taxaUnder.get(childNode);
                node2taxaUnder.put(node, taxaUnder);
                double load = Math.pow(2, taxaUnder.size());
                node.setName(node.getName()+"_"+load);
                node2load.put(node, load);
            }
            else{
                boolean isArticulate = lowestArticulateNodes.contains(node);
                Set<String> taxaUnder = new HashSet<>();
                double load = 1;
                for(NetNode childNode: node.getChildren()){
                    taxaUnder.addAll(node2taxaUnder.get(childNode));
                    if(!isArticulate){
                        double childLoad = node2load.get(childNode);
                        //load *= childLoad==1&&childNode.isLeaf()? 1: childLoad*1.5;
                        load *= childLoad==1&&childNode.isLeaf()? 1: childLoad + 1;
                    }
                }
                if(isArticulate){
                    load = Math.pow(2, taxaUnder.size()/2);
                }
                node2taxaUnder.put(node, taxaUnder);
                node2load.put(node, load);
                //System.out.println(node.getName() + ": " + load);
                node.setName(node.getName()+"_"+load);
                maxLoad = Math.max(maxLoad, load);
            }
        }
        return maxLoad;
    }


    public double computeMaxLoad2(Network network, Map<String, List<String>> species2alleles){
        Network copy = Networks.readNetwork(network.toString());
        Networks.removeAllParameters(copy);
        Networks.autoLabelNodes(copy);
        System.out.println(copy.toString());
        double maxLoad = 0;
        Set<NetNode> lowestArticulateNodes = Networks.getLowestArticulationNodes(network);
        Map<NetNode, Double> node2load = new HashMap<>();
        Map<NetNode, Set<String>> node2taxaUnder = new HashMap<>();
        for(Object nodeO: Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeO;
            if(node.isLeaf()){
                Set<String> taxaUnder = new HashSet<>();
                double load = 1;
                if(species2alleles == null){
                    taxaUnder.add(node.getName());
                }
                else{
                    List<String> alleles = species2alleles.get(node.getName());
                    taxaUnder.addAll(alleles);
                    load = Math.pow(2, alleles.size()/2)/2;

                }
                node2taxaUnder.put(node, taxaUnder);
                node2load.put(node, load);
            }
            else if(node.isNetworkNode()){
                NetNode childNode = (NetNode)node.getChildren().iterator().next();
                Set<String> taxaUnder = node2taxaUnder.get(childNode);
                node2taxaUnder.put(node, taxaUnder);
                double load = Math.pow(2, taxaUnder.size()/2) * Math.max(1,Math.log(node2load.get(childNode)));
                node.setName(node.getName()+"_"+load);
                node2load.put(node, load);
            }
            else{
                boolean isArticulate = lowestArticulateNodes.contains(node);
                Set<String> taxaUnder = new HashSet<>();
                double load = 1;
                for(Object o: node.getChildren()){
                    NetNode childNode = (NetNode)o;
                    taxaUnder.addAll(node2taxaUnder.get(childNode));
                    if(!isArticulate){
                        double childLoad = node2load.get(childNode);
                        //load *= childLoad==1&&childNode.isLeaf()? 1: childLoad*1.5;
                        load *= childLoad==1&&childNode.isLeaf()? 1: childLoad + 1;
                    }
                }
                if(isArticulate){
                    load = Math.pow(2, taxaUnder.size()/2)/2;
                }
                node2taxaUnder.put(node, taxaUnder);
                node2load.put(node, load);
                //System.out.println(node.getName() + ": " + load);
                node.setName(node.getName()+"_"+load);
                maxLoad = Math.max(maxLoad, load);
            }
        }
        return maxLoad;
    }
*/

    public int computeMaxLoad(Network network, Map<String, List<String>> species2alleles){
        Network copy = Networks.readNetwork(network.toString());
        Networks.removeAllParameters(copy);
        Networks.autoLabelNodes(copy);
        System.out.println(copy.toString());
        Set<NetNode> articulateNodes = Networks.getLowestArticulationNodes(copy);

        Map<NetNode, Set<String>> node2taxaUnder = new HashMap<>();
        Map<NetNode, List<Tuple<Integer,Integer>>> node2containedArticulateNodeInfo = new HashMap<>();
        for(Object o: Networks.postTraversal(copy)){
            NetNode node = (NetNode)o;
            List<Tuple<Integer,Integer>> containedArticulateNodeInfo = new ArrayList<Tuple<Integer, Integer>>();
            node2containedArticulateNodeInfo.put(node, containedArticulateNodeInfo);
            Set<String> taxaUnder = new HashSet<>();
            node2taxaUnder.put(node, taxaUnder);
            if(node.isLeaf()){
                if(species2alleles==null){
                    taxaUnder.add(node.getName());
                }
                else{
                    taxaUnder.addAll(species2alleles.get(node.getName()));
                }
            }
            else{
                for(Object child: node.getChildren()){
                    taxaUnder.addAll(node2taxaUnder.get(child));
                    containedArticulateNodeInfo.addAll(node2containedArticulateNodeInfo.get(child));
                }

                if(articulateNodes.contains(node)){
                    boolean isRoot = node.isRoot();
                    NetNode parent = null;
                    if(!isRoot) {
                        parent = (NetNode) node.getParents().iterator().next();
                        parent.removeChild(node);
                    }
                    Network subNetwork = new BniNetwork((BniNetNode)node);
                    int numAlleleMappings = computeNumAlleleMappings1(subNetwork, node2taxaUnder);
                    int load = computeLoad(taxaUnder.size(), containedArticulateNodeInfo, numAlleleMappings);
                    containedArticulateNodeInfo.clear();
                    containedArticulateNodeInfo.add(new Tuple(taxaUnder.size(), load));
                    if(!isRoot){
                        pruneToRandomContainingTree(subNetwork);
                        parent.adoptChild(node, NetNode.NO_DISTANCE);
                    }
                    else{
                        return load;
                    }
                }
                else if(node.isRoot()){
                    int load = computeLoad(taxaUnder.size(), containedArticulateNodeInfo, 1);
                    return load;
                }

            }

        }
        return 0;
    }

    private int computeLoad(int numTaxa, List<Tuple<Integer,Integer>> containedArticulateNodeInfo, int numAlleleMappings){
        int load = getMaxACsCount(numTaxa);
        for(Tuple<Integer,Integer> info: containedArticulateNodeInfo){
            load -= getMaxACsCount(info.Item1);
            load += info.Item2;
        }
        load *= numAlleleMappings;

        //load = numAlleleMappings * numTaxa;
        return load;
    }

    private int getMaxACsCount(int numTaxa){
        //return (int)(numTaxa*(numTaxa-1)/2 + Math.pow(2, numTaxa-1)/2);
        return numTaxa*numTaxa;
    }



    private void pruneToRandomContainingTree(Network network){

        for(Object o: network.getNetworkNodes()){
            NetNode node = (NetNode)o;
            NetNode parent = (NetNode) node.getParents().iterator().next();
            parent.removeChild(node);
        }
        Networks.removeBinaryNodes(network);
    }


    private long[] computeCombinatorial(int n){
        long[] result = new long[n];
        for(int i=1; i<=n; i++){
            result[i-1] = ArithmeticUtils.binomialCoefficient(n, i);
        }
        return result;
    }


    private int computeNumAlleleMappings1(Network net, Map<NetNode, Set<String>> node2taxaUnder){
        Map<String, Integer> leaf2count = new HashMap<>();
        STITree mulTree = new STITree<Double>();
        Queue<NetNode<Double>> source = new LinkedList<NetNode<Double>>();
        Queue<TMutableNode> dest = new LinkedList<TMutableNode>();
        source.offer(net.getRoot());
        dest.offer((TMutableNode) mulTree.getRoot());
        while(!source.isEmpty()){
            NetNode<Double> parent = source.poll();
            TMutableNode peer = dest.poll();
            for (NetNode<Double> child : parent.getChildren()) {
                TMutableNode copy = peer.createChild();
                if(child.isLeaf()){
                    String name = child.getName();
                    Integer amount = leaf2count.get(name);
                    if(amount==null){
                        amount = 0;
                    }
                    leaf2count.put(name, ++amount);
                    copy.setName(name + "_" + amount);
                }
                // Continue to iterate over the children of nn and tn.
                source.offer(child);
                dest.offer(copy);
                //index ++;
            }
        }
        System.out.println(mulTree);
        int load = 1;
        for(Integer count: leaf2count.values()){
            //load *= Math.pow(2,count-1);
            load *= count;
        }
        return load;
    }


    private int computeNumAlleleMappings2(Network net, Map<NetNode, Set<String>> node2taxaUnder){
        Map<NetNode, Integer> node2count = new HashMap<>();
        STITree mulTree = new STITree<Double>();
        Queue<NetNode<Double>> source = new LinkedList<NetNode<Double>>();
        Queue<TMutableNode> dest = new LinkedList<TMutableNode>();
        source.offer(net.getRoot());
        dest.offer((TMutableNode) mulTree.getRoot());
        while(!source.isEmpty()){
            NetNode<Double> parent = source.poll();
            TMutableNode peer = dest.poll();
            for (NetNode<Double> child : parent.getChildren()) {
                TMutableNode copy = peer.createChild();
                String name = child.getName();
                Integer amount = node2count.get(child);
                if(amount==null){
                    amount = 0;
                }
                node2count.put(child, ++amount);
                copy.setName(name + "_" + amount);
                
                // Continue to iterate over the children of nn and tn.
                source.offer(child);
                dest.offer(copy);
                //index ++;
            }
        }
        System.out.println(mulTree);



        int load = 1;
        for(Map.Entry<NetNode,Integer> entry: node2count.entrySet()){
            if(entry.getValue()>1) {
                NetNode node = entry.getKey();
                if(node.isLeaf()){
                    //load *= entry.getValue();
                    load *= Math.pow(2,entry.getValue()-1);
                }
                else if(node.isNetworkNode()){

                    //load *= node2taxaUnder.get(node).size();//Math.max(1,Math.log(node2taxaUnder.get(node).size()));
                }
            }
        }
        return load;
    }

}
