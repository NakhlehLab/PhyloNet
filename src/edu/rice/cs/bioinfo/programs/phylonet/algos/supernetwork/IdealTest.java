package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import com.google.common.collect.Lists;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.reducing.GenerateTrinets3Sets;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 10/4/18
 * Time: 10:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class IdealTest {
    public static Network PrepareNetworkFromBeastString(String input) {
        Network<NetNodeInfo> net = null;
        try {
            String netstring = Pipeline.convertBeastNetworkString(input);
            net = Networks.readNetwork(netstring);
            if(net == null)
                return null;
        } catch (Exception e) {
            return null;
        }

        NetNode rootChild = net.getRoot().getChildren().iterator().next();
        NetNode root = net.getRoot();
        // Reset height
        rootChild.setParentDistance(root, rootChild.getParentDistance(root) + 0.25);

//        // Scale branch lengths
//        for(NetNode<NetNodeInfo> node : net.dfs()) {
//            for(NetNode parent : node.getParents()) {
//                node.setParentDistance(parent, node.getParentDistance(parent) * 200.0);
//            }
//        }

        // Reset inheritance probabilities
        double prob = 0.6;
        for(NetNode<NetNodeInfo> node : net.dfs()) {
            if(node.isNetworkNode()) {
                Iterator<NetNode<NetNodeInfo>> it = node.getParents().iterator();
                NetNode parent1 = it.next();
                NetNode parent2 = it.next();
                node.setParentProbability(parent1, prob);
                node.setParentProbability(parent2, 1.0 - prob);
                prob += 0.01;
            }
        }

        //Pipeline.initNetHeights(net);


        BniNetNode<NetNodeInfo> outgroup = new BniNetNode<>();
        outgroup.setData(new NetNodeInfo(0.0));
        outgroup.setName("Z");
        net.getRoot().adoptChild(outgroup, 2000);

        // Reset heights
        List<NetNode<NetNodeInfo>> internalNodes = new ArrayList<>();
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(net)) {
            if(!node.isLeaf()) {
                internalNodes.add(node);
            } else {
                node.setData(new NetNodeInfo(0.0));
            }
        }
        //Collections.sort(internalNodes, (NetNode<NetNodeInfo> a, NetNode<NetNodeInfo> b)->Double.compare(a.getData().getHeight(), b.getData().getHeight()));

        double curheight = 1.0;
        double increment = 1.0;
        for(NetNode<NetNodeInfo> node : internalNodes) {
            if(!node.isLeaf()) {
                node.setData(new NetNodeInfo(curheight));
                curheight *= 1.2;//+= increment;
                increment += 0.2;
            }
        }

        net.getRoot().getData().setHeight(3000.0);

        // Reset branch lengths
        for(NetNode<NetNodeInfo> node : net.dfs()) {
            for(NetNode<NetNodeInfo> parent : node.getParents()) {
                node.setParentDistance(parent, parent.getData().getHeight() - node.getData().getHeight());
            }
        }

        System.out.println("Full string: " + Networks.getFullString(net));
        System.out.println("Dendroscope: " + Networks.getDendroscopeCompatibleString(net));

        return net;
    }

    public static boolean isD3Net(Network network) {
        for(Object nodeObj : network.dfs()) {
            NetNode node = (NetNode) nodeObj;
            if(node.isNetworkNode()) {
                if(SuperNetwork3.getReticulationNodeDiameter(node) <= 4) {
                    return true;
                }
            }
        }
        return false;
    }

    static void test(String args[]) {
        String filename = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/no_parallel_edges.trees";

        if(args.length > 1) {
            filename = args[1];
        }
        long starttime = System.currentTimeMillis();

        Map<Integer, Integer> countByReti = new TreeMap<>();
        Map<Integer, Integer> correctnessByReti = new TreeMap<>();
        int count = 0;
        int total = 0;
        int exceptions = 0;
        int readfailed = 0;
        int d3net = 0;
        List<String> problematic = new ArrayList<>();
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String s;
            int index = 0;
            while((s = in.readLine()) != null) {
                //if(index < 8500) {index++;continue;}
                //if(index == 682) {index++;continue;}

                try {
                    Network trueNetwork = PrepareNetworkFromBeastString(s);
                    total++;
                    if(trueNetwork == null) {
                        index++;readfailed++;continue;
                    }

                    int numReti = trueNetwork.getReticulationCount();
                    if(!countByReti.containsKey(numReti)) {
                        countByReti.put(numReti, 0);
                        correctnessByReti.put(numReti, 0);
                    }

                    countByReti.put(numReti, countByReti.get(numReti) + 1);

                    //if(true) continue;
                    //if(isD3Net(trueNetwork)) {index++;d3net++;continue;}

                    //if (trueNetwork.getReticulationCount() > 5) {
                    //    index++;
                    //    continue;
                    //}

                    List<Network> subnetworks = NetworkUtils.genAllSubNetworks(trueNetwork, 3);
                    List<Tuple3<Network, String, Double>> newnetlist = new ArrayList<>();

                    for (Network net : subnetworks) {
                        newnetlist.add(new Tuple3<>(net, "", 1.0));
                        //if(net.findNode("I") == null) continue;
                        //if(net.findNode("T") == null) continue;
                        //if(net.findNode("B") == null) continue;
                        //System.out.println(net.toString());
                    }

                    SuperNetwork3.reconcileHeights = false;
                    SuperNetwork3.printDetails_ = true;
                    // Compute best score - only one
                    // component size <= 5
                    SuperNetwork3 superNetwork3 = new SuperNetwork3(newnetlist);
                    Network result = superNetwork3.compute();

                    boolean correctness = Networks.hasTheSameTopology(result, trueNetwork);




                    System.out.println(correctness);
                    System.out.println(index);
                    if (correctness) {
                        count++;
                        correctnessByReti.put(numReti, correctnessByReti.get(numReti) + 1);
                    }
                    else problematic.add(trueNetwork.toString());
                    System.out.println("Correct: " + count);
                    System.out.println("Exceptions: " + exceptions);
                    System.out.println("Read failed: " + readfailed);
                    System.out.println("D3Net: " + d3net);
                    //if(!correctness) break;

                } catch (Exception e) {
                    e.printStackTrace();
                    problematic.add(s);
                    System.out.println("!!!!! " + total);
                    exceptions++;
                    //break;
                }
                index++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("Problematic networks:");
        for(String s : problematic) {
            System.out.println(s);
        }
        System.out.println("Correct: " + count);
        System.out.println("Exceptions: " + exceptions);
        System.out.println("Read failed: " + readfailed);

        for(Integer numReti : countByReti.keySet()) {
            System.out.println("NumReti: " + numReti);
            System.out.println("Count: " + countByReti.get(numReti));
            System.out.println("Correctness: " + correctnessByReti.get(numReti));
        }

        System.out.println("Time (s): " + (System.currentTimeMillis()-starttime)/1000.0);
    }

    static void testOneNetwork(String args[]) {
        String netstring = "((((((C:1.0)#H1:1.0::0.5,B:2.0):2.0)#H2:1.0::0.5,A:5.0):3.0,((D:3.0,#H1:2.0::0.5):3.0,#H2:2.0::0.5):2.0):2.0,E:10.0);";
        Network trueNetwork = Networks.readNetwork(netstring);
        List<Network> subnetworks = NetworkUtils.genAllSubNetworks(trueNetwork, 3);
        List<Tuple3<Network, String, Double>> newnetlist = new ArrayList<>();

        Random random = new Random(12345);
        for (Network net : subnetworks) {
            NetworkUtils.alterHeights(net, random, 0.1);
            if(net.findNode("A") != null && net.findNode("B") != null && net.findNode("C")!= null) {
                NetNode node = ((NetNode)net.getNetworkNodes().iterator().next());
                ((NetNode)node.getParents().iterator().next()).removeChild(node);
                Networks.removeBinaryNodes(net);
            }

            if(net.findNode("A") != null && net.findNode("C") != null && net.findNode("E")!= null) {
                NetNode node = ((NetNode)net.getNetworkNodes().iterator().next());
                ((NetNode)node.getParents().iterator().next()).removeChild(node);
                Networks.removeBinaryNodes(net);
            }

            if(net.findNode("B") != null && net.findNode("C") != null && net.findNode("E")!= null) {
                NetNode node = ((NetNode)net.getNetworkNodes().iterator().next());
                ((NetNode)node.getParents().iterator().next()).removeChild(node);
                Networks.removeBinaryNodes(net);
            }

            newnetlist.add(new Tuple3<>(net, "", 1.0));
            System.out.println(net.toString());
        }

        SuperNetwork3.reconcileHeights = true;
        SuperNetwork3.printDetails_ = true;
        SuperNetwork3.DemoMode = true;
        SuperNetwork3.outgroup = "E";

        SuperNetwork3 superNetwork3 = new SuperNetwork3(newnetlist);
        Network result = superNetwork3.compute();

        boolean correctness = Networks.hasTheSameTopology(result, trueNetwork);


        System.out.println(correctness);

    }

    static void testBigNets(String args[]) {
        String filename = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/bignetworks/no_parallel_edges_80.trees";
        String gtDir = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/bignetworks/80_gts";
        String msPath = null;//"/Users/zhujiafan/Documents/Luay/msdir/ms";
        if(args.length > 1) {
            filename = args[1];
        }

        if(args.length > 2) {
            gtDir = args[2].equals("null") ? null : args[2];
        }

        if(args.length > 3) {
            msPath = args[3].equals("null") ? null : args[3];
        }

        long starttime = System.currentTimeMillis();

        Map<Integer, Integer> countByReti = new TreeMap<>();
        Map<Integer, Integer> correctnessByReti = new TreeMap<>();
        int count = 0;
        int total = 0;
        int exceptions = 0;
        int readfailed = 0;
        int d3net = 0;
        List<String> problematic = new ArrayList<>();
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String s;
            int index = 0;
            while((s = in.readLine()) != null) {

                try {
                    Network trueNetwork = Networks.readNetwork(s);
                    total++;
                    if(trueNetwork == null) {
                        index++;readfailed++;continue;
                    }
                    //if(index < 24) {index++;continue;}

                    NetNode rootChild = (NetNode) trueNetwork.getRoot().getChildren().iterator().next();
                    NetNode root = trueNetwork.getRoot();

                    BniNetNode<NetNodeInfo> outgroup = new BniNetNode<>();
                    outgroup.setData(new NetNodeInfo(0.0));
                    outgroup.setName("Z");
                    trueNetwork.getRoot().adoptChild(outgroup, 0.2);

                    for(Object nodeObj : trueNetwork.dfs()) {
                        NetNode node = (NetNode) nodeObj;
                        for(Object parentObj : node.getParents()) {
                            NetNode parent = (NetNode) parentObj;
                            node.setParentDistance(parent, node.getParentDistance(parent) * 200.0);
                            if(node.isNetworkNode()) {
                                node.setParentProbability(parent, 0.5);
                            }
                        }
                    }

                    List<Network> subnetworks = new ArrayList<>();
                    // Generate gene trees and triplets
                    if(msPath != null) {
                        GenerateTrinets3Sets gen = new GenerateTrinets3Sets();

                        SimGTInNetworkByMS simulator = new SimGTInNetworkByMS();
                        List<Tree> gts = simulator.generateGTs(trueNetwork, null, 1000, msPath);
                        PrintWriter writer = new PrintWriter(gtDir + "/genetree_" + total + ".gts");
                        for(Tree tree : gts) {
                            gen.addTree(tree.toNewick());
                            writer.println(tree.toNewick());
                        }
                        writer.close();

                        writer = new PrintWriter(gtDir + "/triplets_" + total + ".txt");
                        Set<String[]> sets = gen.generate();
                        for(String[] set : sets) {
                            writer.println(set[0] + " " + set[1] + " " + set[2]);
                        }
                        writer.close();
                        continue;
                    }

                    // Read triplets
                    if(gtDir != null) {


                        BufferedReader reader = new BufferedReader(new FileReader(gtDir + "/triplets_" + total + ".txt"));
                        List<List<String>> triplets = new ArrayList<>();
                        while((s = reader.readLine()) != null) {
                            String[] ss = s.split(" ");
                            triplets.add(Arrays.asList(ss));
                        }

                        System.out.println("Reduced set: " + triplets.size());
                        for(List<String> triplet : triplets) {
                            System.out.println(triplet.get(0) + " " + triplet.get(1) + " " + triplet.get(2));
                            subnetworks.add(NetworkUtils.getSubNetwork(trueNetwork, triplet, true));
                        }
                    } else {
                        subnetworks = NetworkUtils.genAllSubNetworks(trueNetwork, 3);
                    }

                    int numReti = trueNetwork.getReticulationCount();
                    if(!countByReti.containsKey(numReti)) {
                        countByReti.put(numReti, 0);
                        correctnessByReti.put(numReti, 0);
                    }

                    countByReti.put(numReti, countByReti.get(numReti) + 1);


                    List<Tuple3<Network, String, Double>> newnetlist = new ArrayList<>();

                    for (Network net : subnetworks) {
                        newnetlist.add(new Tuple3<>(net.clone(), "", 1.0));
                    }

                    SuperNetwork3.reconcileHeights = false;
                    SuperNetwork3.printDetails_ = true;
                    SuperNetwork3.eps = 0.0001;
                    // Compute best score - only one
                    // border size <= 5
                    SuperNetwork3 superNetwork3 = new SuperNetwork3(newnetlist);
                    if(gtDir != null) {
                        List<Set<String>> requiredTrinets = superNetwork3.FindMoreRequiredTrinets();
                        List<String> buildOrder = superNetwork3.GetBuildOrder();
                        Set<String> backboneLeaves = superNetwork3.GetBackboneLeaves();

                        while(requiredTrinets.size() > 0) {

                            System.out.println("Need more trinets: " + requiredTrinets.size());
                            for (Set<String> requiredTaxa : requiredTrinets) {
                                List<String> taxa = new ArrayList<>(requiredTaxa);
                                System.out.println(taxa.get(0) + " " + taxa.get(1) + " " + taxa.get(2));
                                subnetworks.add(NetworkUtils.getSubNetwork(trueNetwork, taxa, true));
                            }
                            newnetlist.clear();
                            for (Network net : subnetworks) {
                                newnetlist.add(new Tuple3<>(net.clone(), "", 1.0));
                            }
                            superNetwork3 = new SuperNetwork3(newnetlist);
                            requiredTrinets = superNetwork3.FindMoreRequiredTrinets();
                            buildOrder = superNetwork3.GetBuildOrder();
                            backboneLeaves = superNetwork3.GetBackboneLeaves();
                        }
                        superNetwork3 = new SuperNetwork3(newnetlist);
                        superNetwork3.SetBuildOrder(buildOrder);
                        superNetwork3.SetBackboneLeaves(backboneLeaves);
                    }
                    //SuperNetwork3.DemoMode = true;
                    Network result = superNetwork3.compute();

                    boolean correctness = Networks.hasTheSameTopology(result, trueNetwork);




                    System.out.println(correctness);
                    System.out.println(index);
                    if (correctness) {
                        count++;
                        correctnessByReti.put(numReti, correctnessByReti.get(numReti) + 1);
                    }
                    else problematic.add(trueNetwork.toString());
                    System.out.println("Correct: " + count);
                    System.out.println("Exceptions: " + exceptions);
                    System.out.println("Read failed: " + readfailed);
                    System.out.println("D3Net: " + d3net);
                    //if(!correctness) break;
                    //if(true) break;
                } catch (Exception e) {
                    e.printStackTrace();
                    problematic.add(s);
                    System.out.println("!!!!! " + total + " " + index);
                    exceptions++;
                    break;
                }
                index++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("Problematic networks:");
        for(String s : problematic) {
            System.out.println(s);
        }
        System.out.println("Correct: " + count);
        System.out.println("Exceptions: " + exceptions);
        System.out.println("Read failed: " + readfailed);

        for(Integer numReti : countByReti.keySet()) {
            System.out.println("NumReti: " + numReti);
            System.out.println("Count: " + countByReti.get(numReti));
            System.out.println("Correctness: " + correctnessByReti.get(numReti));
        }

        System.out.println("Time (s): " + (System.currentTimeMillis()-starttime)/1000.0);
    }

    public static void main(String[] args) {
        System.out.println(args.length);
        for(int i = 0 ; i < args.length ; i++) {
            System.out.println(args[i]);
        }

        if(args.length > 0) {
            if(args[0].equals("BigNets")) {
                testBigNets(args);
            } else if(args[0].equals("SmallNets")) {
                test(args);
            } else if(args[0].equals("Demo")) {
                testOneNetwork(args);
            }
        }

    }
}
