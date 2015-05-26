package edu.rice.cs.bioinfo.programs.phylonet.structs.network;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.ByteArrayInputStream;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 3/19/13
 * Time: 2:02 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkMetricNakhleh<T> {
    Network<T> _network1;
    Network<T> _network2;


    public NetworkMetricNakhleh(){}



    public NetworkMetricNakhleh(Network<T> net1, Network<T> net2){
        _network1 = net1;
        _network2 = net2;
    }


    private void computeConvergentSet(Network<T> network, Map<NetNode<T>, Set<NetNode<T>>> node2reachableLeaves){
        Map<Set<NetNode<T>>, List<NetNode<T>>> reachableLeaves2nodes = new HashMap<Set<NetNode<T>>, List<NetNode<T>>>();
        for(NetNode<T> node: Networks.postTraversal(network)){
            Set<NetNode<T>> reachableLeaves = new HashSet<NetNode<T>>();
            if(node.isLeaf()){
                reachableLeaves.add(node);
            }
            else{
                for(NetNode<T> child: node.getChildren()){
                    reachableLeaves.addAll(node2reachableLeaves.get(child));
                }
            }
            node2reachableLeaves.put(node, reachableLeaves);
            List<NetNode<T>> equivalantNodes = reachableLeaves2nodes.get(reachableLeaves);
            if(equivalantNodes == null){
                equivalantNodes = new ArrayList<NetNode<T>>();
                reachableLeaves2nodes.put(reachableLeaves, equivalantNodes);
            }
            equivalantNodes.add(node);
        }


        for(Map.Entry<Set<NetNode<T>>, List<NetNode<T>>> entry: reachableLeaves2nodes.entrySet()){
            if(entry.getValue().size() == 1){
                node2reachableLeaves.remove(entry.getValue().iterator().next());
            }
        }
    }

    public double computeDistanceBetweenTwoNetworks(Network<T> network1, Network<T> network2){
        Map<NetNode<T>, MutableTuple<Integer, NetNode<T>>> KInfoInNet1 = new HashMap<NetNode<T>, MutableTuple<Integer, NetNode<T>>>();
        Map<NetNode<T>, MutableTuple<Integer, NetNode<T>>> KInfoInNet2 = new HashMap<NetNode<T>, MutableTuple<Integer, NetNode<T>>>();
        computeEquivalentNodeSets(network1, network2, KInfoInNet1, KInfoInNet2);
        return computeDistance(KInfoInNet1, KInfoInNet2);
    }

    /*
    private void computeConvergentSet(Network<T> network, List<Tuple<Set<String>, List<NetNode>>> convergentSetList){
        Map<Set<String>, List<NetNode>> reachableLeaves2nodes = new HashMap<Set<String>, List<NetNode>>();
        Map<NetNode, Set<String>> node2reachableLeaves = new HashMap<NetNode, Set<String>>();
        for(NetNode<T> node: Networks.postTraversal(network)){
            Set<String> reachableLeaves = new HashSet<String>();
            if(node.isLeaf()){
                reachableLeaves.add(node.getName());
            }
            else{
                for(NetNode<T> child: node.getChildren()){
                    reachableLeaves.addAll(node2reachableLeaves.get(child));
                }
            }
            node2reachableLeaves.put(node, reachableLeaves);
            List<NetNode> equivalantNodes = reachableLeaves2nodes.get(reachableLeaves);
            if(equivalantNodes == null){
                equivalantNodes = new ArrayList<NetNode>();
                reachableLeaves2nodes.put(reachableLeaves, equivalantNodes);
            }
            equivalantNodes.add(node);
        }

        for(Map.Entry<Set<String>, List<NetNode>> entry: reachableLeaves2nodes.entrySet()){
            if(entry.getValue().size() == 1){
                continue;
            }
            Set<String> reachableLeaves = entry.getKey();
            int size = reachableLeaves.size();
            int index = 0;
            for(Tuple<Set<String>, List<NetNode>> tuple: convergentSetList){
                if(size > tuple.Item1.size()){
                    break;
                }
                index ++;
            }
            convergentSetList.add(index, new Tuple<Set<String>, List<NetNode>>(reachableLeaves, entry.getValue()));
        }
    }
    */



    public Network<T> computeReducedNetwork(Network<T> network){
        Network<T> reducedNetwork = cloneNetwork(network);
        List<Tuple<NetNode<T>, NetNode<T>>> symbolicLeaf2subtree = new ArrayList<Tuple<NetNode<T>, NetNode<T>>>();
        condenseMaximalSubtrees(reducedNetwork, symbolicLeaf2subtree);
        reduceNetwork(reducedNetwork);
        deCondenseMaximalSubtrees(symbolicLeaf2subtree);
        return reducedNetwork;
    }


    private String network2string(Network speciesNetwork){
        RnNewickPrinter<Double> rnNewickPrinter = new RnNewickPrinter<Double>();
        StringWriter sw = new StringWriter();
        rnNewickPrinter.print(speciesNetwork, sw);
        return sw.toString();
    }


    private void reduceNetwork(Network<T> network){
        Map<NetNode<T>, Set<NetNode<T>>> node2ReachableLeaves = new HashMap<NetNode<T>, Set<NetNode<T>>>();
        computeConvergentSet(network, node2ReachableLeaves);
        System.out.println(network2string(network)+"\n");
        for(NetNode<T> node: Networks.postTraversal(network)){
            Set<NetNode<T>> reachableLeaves = node2ReachableLeaves.get(node);
            if(reachableLeaves!=null && !node.isLeaf()){
                System.out.println("after handling " + node.getName());
                //node.removeAllChild();
                List<NetNode<T>> children = new ArrayList<NetNode<T>>();
                for(NetNode<T> child: node.getChildren()){
                    children.add(child);
                }
                for(NetNode<T> child: children){
                    System.out.println("removing " + child.getName());
                    child.removeItself();
                }
                for(NetNode<T> reachableLeaf: reachableLeaves){
                    node.adoptChild(reachableLeaf, NetNode.NO_DISTANCE);
                }
                Networks.removeBinaryNodes(network);
                System.out.println(network2string(network) + "\n");
            }

        }
        Networks.removeBinaryNodes(network);
    }

    private void condenseMaximalSubtrees(Network<T> network, List<Tuple<NetNode<T>, NetNode<T>>> symbolicLeaf2subtree){
        Set<NetNode<T>> treeNodes = new HashSet<NetNode<T>>();
        Set<NetNode<T>> maxSubtrees = new HashSet<NetNode<T>>();
        for(NetNode<T> node: Networks.postTraversal(network)){
            if(node.isLeaf()){
                treeNodes.add(node);
            }
            else{
                boolean allTreeNodes = node.isTreeNode();
                boolean containsTreeNodes = false;
                for(NetNode<T> child: node.getChildren()){
                    if(!treeNodes.contains(child)){
                        allTreeNodes = false;
                    }
                    else {
                        containsTreeNodes = true;
                    }
                }
                if(allTreeNodes){
                    if(!node.isNetworkNode()){
                        treeNodes.add(node);
                    }
                }else if(containsTreeNodes){
                    for(NetNode<T> child: node.getChildren()){
                        if(treeNodes.contains(child)){
                            maxSubtrees.add(child);
                        }
                    }
                }
            }
        }
        int index = 0;
        for(NetNode<T> subtreeRoot: maxSubtrees){
            NetNode<T> parent = subtreeRoot.getParents().iterator().next();
            double distance = subtreeRoot.getParentDistance(parent);
            parent.removeChild(subtreeRoot);
            BniNetNode<T> symbolicNode = new BniNetNode<T>();
            symbolicNode.setName("SymbolicNode" + (index++));
            parent.adoptChild(symbolicNode, distance);
            symbolicLeaf2subtree.add(new Tuple<NetNode<T>, NetNode<T>>(symbolicNode, subtreeRoot));
        }
    }



    private void deCondenseMaximalSubtrees(List<Tuple<NetNode<T>, NetNode<T>>> symbolicLeaf2subtree){
        for(Tuple<NetNode<T>, NetNode<T>> tuple: symbolicLeaf2subtree){
            NetNode<T> symbolicNode = tuple.Item1;
            NetNode<T> subtreeRoot = tuple.Item2;
            //NetNode<T> parent = symbolicNode.getParents().iterator().next();
            //double distance = symbolicNode.getParentDistance(parent);
            //parent.removeChild(symbolicNode);
            symbolicNode.adoptChild(subtreeRoot,NetNode.NO_DISTANCE);
        }
    }



    private Network<T> cloneNetwork(Network<T> network){
        Network<T> clonedNetwork = null;
        try{
            StringWriter writer = new StringWriter();
            RnNewickPrinter printer = new RnNewickPrinter();
            printer.print(network, writer);
            RichNewickReaderAST reader = new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
            reader.setHybridSumTolerance(BigDecimal.valueOf(0.0001));
            NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
            RichNewickReadResult<edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks> readResult = reader.read(new ByteArrayInputStream(writer.toString().getBytes()));
            clonedNetwork = transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return clonedNetwork;
    }


/*
    private void computeEquivalentNodeSetsWithinANetwork(Network<T> network, Map<String, NetNode<T>> leafName2node, Map<Set<Tuple<NetNode<T>,Integer>>, NetNode> children2node, Map<NetNode<T>, NetNode<T>> node2repNode){
        //Map<Integer, List<NetNode>> childNumber2Nodes = new HashMap<Integer, List<NetNode>>();
        //Map<Set<Tuple<NetNode<T>,Integer>>, NetNode> children2node = new HashMap<Set<Tuple<NetNode<T>,Integer>>, NetNode>();
        //Map<NetNode<T>, NetNode<T>> node2repNode = new HashMap<NetNode<T>, NetNode<T>>();
        
        for(NetNode<T> node: Networks.postTraversal(network)){
            Set<Tuple<NetNode<T>,Integer>> childrenSet = new HashSet<Tuple<NetNode<T>, Integer>>();
            if(node.isLeaf()){
                leafName2node.put(node.getName(), node);
                node2repNode.put(node, node);
                childrenSet.add(new Tuple<NetNode<T>,Integer>(node, 0)); 
                children2node.put(childrenSet, node);
            }else{               
                Map<NetNode<T>, Integer> childNode2count = new HashMap<NetNode<T>, Integer>();
                for(NetNode<T> child: node.getChildren()){
                    NetNode<T> equivalentNode = node2repNode.get(child);
                    Integer count = childNode2count.get(equivalentNode);
                    if(count == null){
                        count = 0;
                    }
                    childNode2count.put(equivalentNode, count + 1);
                }
                for(Map.Entry<NetNode<T>, Integer> entry: childNode2count.entrySet()){
                    childrenSet.add(new Tuple<NetNode<T>,Integer>(entry.getKey(), entry.getValue()));
                }
                NetNode<T> equivalentNode = children2node.get(childrenSet);
                if(equivalentNode == null){
                    node2repNode.put(node, node);
                    children2node.put(childrenSet, node);
                }
                else{
                    node2repNode.put(node, equivalentNode);
                }
            }
            
        }
    }
*/

    //
    private void computeEquivalentNodeSets(Network<T> network1, Network<T> network2, Map<NetNode<T>, MutableTuple<Integer, NetNode<T>>> KInfoInNet1, Map<NetNode<T>, MutableTuple<Integer, NetNode<T>>> KInfoInNet2){
        //Map<Integer, List<NetNode>> childNumber2Nodes = new HashMap<Integer, List<NetNode>>();
        
        Map<Set<Tuple<NetNode<T>,Integer>>, NetNode> children2nodeInNet1 = new HashMap<Set<Tuple<NetNode<T>,Integer>>, NetNode>();
        Map<NetNode<T>, NetNode<T>> node2repNodeInNet1 = new HashMap<NetNode<T>, NetNode<T>>();
        Map<String, NetNode<T>> leafName2nodeInNet1 = new HashMap<String, NetNode<T>>();
        for(NetNode<T> node: Networks.postTraversal(network1)){
            Set<Tuple<NetNode<T>,Integer>> childrenSet = new HashSet<Tuple<NetNode<T>, Integer>>();
            if(node.isLeaf()){
                leafName2nodeInNet1.put(node.getName(), node);
                node2repNodeInNet1.put(node, node);
                childrenSet.add(new Tuple<NetNode<T>,Integer>(node, 0));
                children2nodeInNet1.put(childrenSet, node);
            }else{
                Map<NetNode<T>, Integer> childNode2count = new HashMap<NetNode<T>, Integer>();
                for(NetNode<T> child: node.getChildren()){
                    NetNode<T> equivalentNode = node2repNodeInNet1.get(child);
                    Integer count = childNode2count.get(equivalentNode);
                    if(count == null){
                        count = 0;
                    }
                    childNode2count.put(equivalentNode, count + 1);
                }
                for(Map.Entry<NetNode<T>, Integer> entry: childNode2count.entrySet()){
                    childrenSet.add(new Tuple<NetNode<T>,Integer>(entry.getKey(), entry.getValue()));
                }
                NetNode<T> equivalentNode = children2nodeInNet1.get(childrenSet);
                if(equivalentNode == null){
                    node2repNodeInNet1.put(node, node);
                    children2nodeInNet1.put(childrenSet, node);
                }
                else{
                    node2repNodeInNet1.put(node, equivalentNode);
                }
            }
        }

        /*
        System.out.print(children2nodeInNet1.values().size() + " : ");
        for(NetNode<T> node: children2nodeInNet1.values()){
            System.out.print(node.getName() + " ");
        }
        System.out.println();
        */

        Map<Set<Tuple<NetNode<T>,Integer>>, NetNode> children2nodeInNet2 = new HashMap<Set<Tuple<NetNode<T>,Integer>>, NetNode>();
        Map<NetNode<T>, NetNode<T>> node2repNodeInNet2 = new HashMap<NetNode<T>, NetNode<T>>();
        Map<NetNode<T>, NetNode<T>> node2TOnode1 = new HashMap<NetNode<T>, NetNode<T>>();
        Map<NetNode<T>, NetNode<T>> node1TOnode2 = new HashMap<NetNode<T>, NetNode<T>>();
        for(NetNode<T> node: Networks.postTraversal(network2)){
            Set<Tuple<NetNode<T>,Integer>> childrenSetInNet2 = new HashSet<Tuple<NetNode<T>, Integer>>();
            if(node.isLeaf()){
                node2repNodeInNet2.put(node, node);
                childrenSetInNet2.add(new Tuple<NetNode<T>,Integer>(node, 0));
                children2nodeInNet2.put(childrenSetInNet2, node);

                //TO build the common set
                NetNode<T> equivalentNodeInNet1 = leafName2nodeInNet1.get(node.getName());
                if(equivalentNodeInNet1 != null){
                    node2TOnode1.put(node, equivalentNodeInNet1);
                    node1TOnode2.put(equivalentNodeInNet1, node);
                }
            }else{
                Map<NetNode<T>, Integer> childNode2countInNet2 = new HashMap<NetNode<T>, Integer>();
                Map<NetNode<T>, Integer> childNode2countInNet1 = new HashMap<NetNode<T>, Integer>();
                for(NetNode<T> child: node.getChildren()){
                    NetNode<T> equivalentNodeInNet2 = node2repNodeInNet2.get(child);
                    Integer countInNet2 = childNode2countInNet2.get(equivalentNodeInNet2);
                    if(countInNet2 == null){
                        countInNet2 = 0;
                    }
                    childNode2countInNet2.put(equivalentNodeInNet2, countInNet2 + 1);

                    //TO build the common set
                    if(childNode2countInNet1!=null){
                        NetNode<T> equivalentNodeInNet1 = node2TOnode1.get(equivalentNodeInNet2);
                        if(equivalentNodeInNet1 == null){
                            childNode2countInNet1 = null;
                        }
                        else{
                            Integer countInNet1 = childNode2countInNet1.get(equivalentNodeInNet1);
                            if(countInNet1 == null){
                                countInNet1 = 0;
                            }
                            childNode2countInNet1.put(equivalentNodeInNet1, countInNet1 + 1);
                        }
                    }
                }
                for(Map.Entry<NetNode<T>, Integer> entry: childNode2countInNet2.entrySet()){
                    childrenSetInNet2.add(new Tuple<NetNode<T>,Integer>(entry.getKey(), entry.getValue()));
                }
                NetNode<T> equivalentNodeInNet2 = children2nodeInNet2.get(childrenSetInNet2);
                if(equivalentNodeInNet2 == null){
                    node2repNodeInNet2.put(node, node);
                    children2nodeInNet2.put(childrenSetInNet2, node);

                    //TO build the common set
                    if(childNode2countInNet1 != null){
                        Set<Tuple<NetNode<T>,Integer>> childrenSetInNet1 = new HashSet<Tuple<NetNode<T>, Integer>>();
                        for(Map.Entry<NetNode<T>, Integer> entry: childNode2countInNet1.entrySet()){
                            childrenSetInNet1.add(new Tuple<NetNode<T>,Integer>(entry.getKey(), entry.getValue()));
                        }
                        NetNode<T> equivalentNodeInNet1 = children2nodeInNet1.get(childrenSetInNet1);
                        if(equivalentNodeInNet1 != null){
                            node2TOnode1.put(node, equivalentNodeInNet1);
                            node1TOnode2.put(equivalentNodeInNet1, node);
                        }
                    }
                }
                else{
                    node2repNodeInNet2.put(node, equivalentNodeInNet2);
                }
            }

        }
/*
        System.out.print(children2nodeInNet2.values().size() + " : ");
        for(NetNode<T> node: children2nodeInNet2.values()){
            System.out.print(node.getName() + " ");
        }
        System.out.println();
        System.out.println(node1TOnode2.size() + " : ");
        for(Map.Entry<NetNode<T>, NetNode<T>> entry: node1TOnode2.entrySet()){
            System.out.println(entry.getKey().getName() + " = " + entry.getValue().getName());
        }
        System.out.println();
*/
        //Summarize Result
        for(Map.Entry<NetNode<T>, NetNode<T>> entry: node2repNodeInNet1.entrySet()){
            NetNode<T> uniqueNodeInNet1 = entry.getValue();
            MutableTuple<Integer, NetNode<T>> KInfo = KInfoInNet1.get(uniqueNodeInNet1);
            if(KInfo == null){
                KInfo = new MutableTuple<Integer, NetNode<T>>(1, node1TOnode2.get(uniqueNodeInNet1));
                KInfoInNet1.put(uniqueNodeInNet1, KInfo);
            }
            else{
                KInfo.Item1++;
            }
        }

        for(Map.Entry<NetNode<T>, NetNode<T>> entry: node2repNodeInNet2.entrySet()){
            NetNode<T> uniqueNodeInNet2 = entry.getValue();
            MutableTuple<Integer, NetNode<T>> KInfo = KInfoInNet2.get(uniqueNodeInNet2);
            if(KInfo == null){
                KInfo = new MutableTuple<Integer, NetNode<T>>(1, node2TOnode1.get(uniqueNodeInNet2));
                KInfoInNet2.put(uniqueNodeInNet2, KInfo);
            }
            else{
                KInfo.Item1++;
            }
        }

    }


    private double computeDistance(Map<NetNode<T>, MutableTuple<Integer, NetNode<T>>> KInfoInNet1, Map<NetNode<T>, MutableTuple<Integer, NetNode<T>>> KInfoInNet2){
        int totalDiff1 = 0;
        for(MutableTuple<Integer, NetNode<T>> kInfo: KInfoInNet1.values()){
            int KV = kInfo.Item1;
            int KVPrime = 0;
            if(kInfo.Item2 != null){
                KVPrime = KInfoInNet2.get(kInfo.Item2).Item1;
            }
            totalDiff1 += KV-KVPrime>0 ? KV-KVPrime:0;
        }

        int totalDiff2 = 0;
        for(MutableTuple<Integer, NetNode<T>> kInfo: KInfoInNet2.values()){
            int KV = kInfo.Item1;
            int KVPrime = 0;
            if(kInfo.Item2 != null){
                KVPrime = KInfoInNet1.get(kInfo.Item2).Item1;
            }
            totalDiff2 += KV-KVPrime>0 ? KV-KVPrime:0;
        }

        return (totalDiff1+totalDiff2) / 2;
    }


}