package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import javafx.beans.binding.ObjectExpression;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import sun.nio.ch.Net;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 9/22/18
 * Time: 11:30 AM
 * To change this template use File | Settings | File Templates.
 */
public class SuperNetwork3 {
    class NetworkWithInfo {

        Network<NetNodeInfo> network;
        double PairwiseDistanceSum;
        List<String> taxa;
        double distMat[][];
        boolean dirty = false;
        String backup;
        Map<String, Set<String>> dependencies;
        Map<Tuple<String, String>, List<Double>> ehm; // Extended height matrix
        List<Tuple<Network, Map<NetNode, NetNode>>> binets;
        String filename;
        boolean trustTime = true;
        double percentage = 0.0;

        NetworkWithInfo(Network<NetNodeInfo> net, String name, double percentage) {
            network = net;
            backup = net.toString();
            filename = name;
            this.percentage = percentage;

            taxa = new ArrayList<>();
            for(Object leafObj : network.getLeaves()) {
                NetNode leaf = (NetNode) leafObj;
                taxa.add(leaf.getName());
            }
            Collections.sort(taxa);

            // Compute pairwise height
            distMat = new double[taxa.size()][taxa.size()];
            PairwiseDistanceSum = 0;
            for(int i = 0 ; i < taxa.size() ; i++) {
                for(int j = 0 ; j < taxa.size() ; j++) {
                    if(i == j) {
                        distMat[i][j] = 0;
                    } else {
                        distMat[i][j] = Double.MAX_VALUE;
                        Set<NetNode> MCRAs = getMCRAs(network, taxa.get(i), taxa.get(j));
                        for(NetNode node : MCRAs) {
                            double value = ((NetNodeInfo) node.getData()).getHeight();
                            distMat[i][j] = Math.min(value, distMat[i][j]);
                        }
                        if(i < j)
                            PairwiseDistanceSum += distMat[i][j];
                    }
                }
            }

            ehm = getEHMFromNetwork(network);
            binets = new ArrayList<>();
            for(int i = 0 ; i < taxa.size() ; i++) {
                for(int j = i + 1 ; j < taxa.size() ; j++) {
                    List<String> selected = new ArrayList<>();
                    selected.add(taxa.get(i));
                    selected.add(taxa.get(j));
                    binets.add(getSubNetwork(net, selected, true));
                }
            }

            dependencies = new HashMap<>();
            Map<NetNode, Set<String>> children = new HashMap<>();
            for(int i = 0 ; i < taxa.size() ; i++) {
                NetNode ancenstor = network.findNode(taxa.get(i)).getParents().iterator().next();
                while(!ancenstor.isNetworkNode()) {
                    if(!children.containsKey(ancenstor)) {
                        children.put(ancenstor, new HashSet<>());
                    }
                    children.get(ancenstor).add(taxa.get(i));
                    if(ancenstor.isRoot()) break;
                    ancenstor = (NetNode) ancenstor.getParents().iterator().next();
                }
            }

            Map<String, Integer> numRetiAbove = new HashMap<>();
            for(String leafname : taxa) {
                int k = 0;
                Set<NetNode> ancestors = getAllAncestors(network.findNode(leafname));
                for(NetNode node : ancestors)
                    if(node.isNetworkNode())
                        k++;
                numRetiAbove.put(leafname, k);
            }
            List<String> underReti = Networks.getTaxaNamesUnderReticulation(network);
            for(String targetName : underReti) {
                if(!dependencies.containsKey(targetName))
                    dependencies.put(targetName, new HashSet<>() );

                for(String taxon : taxa) {
                    if(!taxon.equals(targetName) && numRetiAbove.get(taxon) <= numRetiAbove.get(targetName))
                        dependencies.get(targetName).add(taxon);
                }
            }
            /*for(String targetName : underReti) {
                if(!dependencies.containsKey(targetName))
                    dependencies.put(targetName, new HashSet<>() );

                // Find the first ancenstor reticulation node
                NetNode ancenstor = network.findNode(targetName).getParents().iterator().next();
                while(!ancenstor.isNetworkNode()) {
                    if(ancenstor.isRoot()) break;
                    ancenstor = (NetNode) ancenstor.getParents().iterator().next();
                }

                for(Object parentObj : ancenstor.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    if(children.containsKey(parent))
                        dependencies.get(targetName).addAll(children.get(parent));
                }
            }*/
        }

    }

    private int n_;
    List<NetworkWithInfo> subnetworks_;
    List<NetworkWithInfo> subnetworks_removed_;
    List<NetworkWithInfo> subnetworks_reduceed_;
    private List<String> leafnames_;
    private Set<String> leavesUnderReticulation;
    private Set<String> leavesUnderBubbleReticulation;
    private List<String> buildOrder;
    Map<Tuple<String, String>, List<Double>> ehm; // Extended height matrix
    protected static double eps = 0.01;
    protected static boolean trustReticulationTime = true;
    public static boolean printDetails_ = false;
    protected static boolean reconcileHeights = false;
    public static String outgroup = "Z";
    private double popsize = Double.NaN;

    public static void initNetHeights(Network<NetNodeInfo> network) {
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(network)) {
            if(node.getData() == null) {
                node.setData(new NetNodeInfo(0.0));
            }
            for(NetNode<NetNodeInfo> par : node.getParents()) {
                double dist = node.getParentDistance(par);
                if(par.getData() == null) {
                    par.setData(new NetNodeInfo(node.getData().getHeight() + dist));
                }
            }
        }
    }

    public static Set<NetNode> getAllAncestors(NetNode node) {
        Queue<NetNode> q = new ArrayDeque<>();
        Set<NetNode> result = new HashSet<>();
        q.add(node);
        while(!q.isEmpty()) {
            NetNode curnode = q.poll();
            for(Object parentObj : curnode.getParents()) {
                NetNode parent = (NetNode) parentObj;
                result.add(parent);
                q.add(parent);
            }
        }
        return result;
    }

    private static NetNode<NetNodeInfo> getMCRA(List<NetNode> leaves) {
        Set<NetNode> mcras = getMCRAs(leaves);
        NetNode<NetNodeInfo> lowest = null;
        for(NetNode node : mcras) {
            if(lowest == null || lowest.getData().getHeight() > ((NetNodeInfo)node.getData()).getHeight()) {
                lowest = node;
            }
        }
        return lowest;
    }

    private static Set<NetNode> getMCRAs(List<NetNode> leaves) {
        Set<NetNode> curSet = null;
        for(NetNode leaf : leaves) {
            Set<NetNode> s = getAllAncestors(leaf);
            if(curSet == null) {
                curSet = s;
            } else {
                curSet.retainAll(s);
            }
        }
        return curSet;
    }

    public static Set<NetNode> getMCRAs(Network net, String s1, String s2) {
        NetNode l1 = net.findNode(s1);
        NetNode l2 = net.findNode(s2);
        List<NetNode> leaves = new ArrayList<>();
        leaves.add(l1);
        leaves.add(l2);
        return getMCRAs(leaves);
    }

    public static Tuple<Network, Map<NetNode, NetNode>> getSubNetwork(Network trueNetwork, List<String> selectedLeaves, boolean removeBinaryNodes) {
        int subsize = selectedLeaves.size();
        Network currentNetwork = trueNetwork;

        // bottom-up build subnetwork
        Map<NetNode, BniNetNode> old2new = new HashMap<>();
        Map<NetNode, Integer> hitCount = new HashMap<>();
        for(String leaf : selectedLeaves) {
            Set<NetNode> visited = new HashSet<>();
            Queue<NetNode> queue = new LinkedList<>();
            NetNode node = currentNetwork.findNode(leaf);
            BniNetNode newleaf = new BniNetNode();
            newleaf.setName(node.getName());
            if(node.getData() != null) newleaf.setData(new NetNodeInfo(((NetNodeInfo)node.getData()).getHeight()));
            old2new.put(node, newleaf);
            queue.add(node);
            while(queue.size() > 0) {
                node = queue.poll();
                if(visited.contains(node)) continue;
                if(!hitCount.containsKey(node)) {
                    hitCount.put(node, 0);
                }
                hitCount.put(node, hitCount.get(node) + 1);
                visited.add(node);
                for(Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    queue.add(parent);
                    BniNetNode newParentNode;
                    if(old2new.containsKey(parent)) {
                        newParentNode = old2new.get(parent);
                    } else {
                        newParentNode = new BniNetNode();
                        newParentNode.setName(parent.getName());
                        if(parent.getData() != null) newParentNode.setData(new NetNodeInfo(((NetNodeInfo)parent.getData()).getHeight()));
                        old2new.put(parent, newParentNode);
                    }

                    newParentNode.adoptChild(old2new.get(node), node.getParentDistance(parent));
                    old2new.get(node).setParentProbability(newParentNode, node.getParentProbability(parent));
                    old2new.get(node).setParentSupport(newParentNode, node.getParentSupport(parent));
                }
            }
        }

        Network newsubnetwork = new BniNetwork(old2new.get(currentNetwork.getRoot()));
        newsubnetwork.getRoot().setRootPopSize(currentNetwork.getRoot().getRootPopSize());
        if(removeBinaryNodes) {
            Networks.removeBinaryNodes(newsubnetwork);
        }

        Map<NetNode, NetNode> new2old = new HashMap<>();
        for(NetNode oldNode : old2new.keySet()) {
            new2old.put(old2new.get(oldNode), oldNode);
        }



        return new Tuple<>(newsubnetwork, new2old);
    }

    public static int getCountDistance(NetNode node, NetNode ancestor)  {
        int dist = Integer.MAX_VALUE;
        Queue<Tuple<NetNode, Integer>> q = new LinkedList<>();
        q.offer(new Tuple<>(node, 0));
        while(!q.isEmpty()) {
            Tuple<NetNode, Integer> cur = q.poll();
            if(cur.Item1 == ancestor) return cur.Item2;
            for(Object parentObj : cur.Item1.getParents()) {
                NetNode parent = (NetNode) parentObj;
                q.offer(new Tuple<>(parent, cur.Item2 + 1));
                if(parent == ancestor) return cur.Item2 + 1;
            }
        }

        return dist;
    }

    public static List<String> getTaxaNamesUnderReticulation(NetNode reticulationNode) {
        Set<String> leafset = new HashSet<>();
        Queue<NetNode> queue = new LinkedList<>();
        queue.offer(reticulationNode);
        while(!queue.isEmpty()) {
            NetNode node = queue.poll();
            if(node.isLeaf()) {
                leafset.add(node.getName());
            } else {
                for(Object childObj : node.getChildren()) {
                    NetNode child = (NetNode) childObj;
                    queue.offer(child);
                }
            }
        }
        List<String> results = new ArrayList<>(leafset);
        return results;
    }

    public static int getReticulationNodeDiameter(NetNode reticulationNode) {
        Iterator<NetNode> it = reticulationNode.getParents().iterator();
        NetNode parent1 = it.next();
        NetNode parent2 = it.next();
        Set<NetNode> set1 = getAllAncestors(parent1);
        Set<NetNode> set2 = getAllAncestors(parent2);

        Set<NetNode> set = new HashSet<>(set1);
        set.retainAll(set2);
        int result = Integer.MAX_VALUE;
        for(NetNode n: set) {
            result = Math.min(result, getCountDistance(parent1, n) + getCountDistance(parent2, n));
        }
        return result;
    }

    public static <T> Map<NetNode<T>, NetNode<T>> mapTwoNetworks(Network<T> net1, Network<T> net2) {
        Map<NetNode<T>, NetNode<T>> result = new HashMap<NetNode<T>, NetNode<T>>();

        result.put(net1.getRoot(), net2.getRoot());

        for(NetNode<T> node1 : net1.getLeaves()) {
            NetNode<T> node2 = net2.findNode(node1.getName());
            result.put(node1, node2);
        }

        int nodeCount = 0;
        for(NetNode<T> node1 : net1.dfs()) {
            nodeCount++;
        }

        while(result.size() < nodeCount) {
            NetNode<T> nextNodeToSolve = null;
            NetNode<T> nextNode2ToSolve = null;

            NetNode<T> nextNodeParent = null;
            NetNode<T> nextNodeChild = null;
            for(NetNode<T> node1 : result.keySet()) {
                NetNode<T> node2 = result.get(node1);

                Set<NetNode<T>> neighbors1 = new HashSet<>();
                for(NetNode<T> child : node1.getChildren()) neighbors1.add(child);
                for(NetNode<T> parent : node1.getParents()) neighbors1.add(parent);
                neighbors1.removeAll(result.keySet());

                Set<NetNode<T>> neighbors2 = new HashSet<>();
                for(NetNode<T> child : node2.getChildren()) neighbors2.add(child);
                for(NetNode<T> parent : node2.getParents()) neighbors2.add(parent);
                neighbors2.removeAll(result.values());

                if(neighbors1.size() == 1) {
                    if(neighbors2.size() != 1)
                        return null;

                    nextNodeToSolve = neighbors1.iterator().next();
                    nextNode2ToSolve = neighbors2.iterator().next();
                    break;
                }
            }

            if(nextNodeToSolve == null) {
                double closest = Double.MAX_VALUE;

//                for(NetNode<T> node1 : net1.dfs()) {
//                    if(result.containsKey(node1)) continue;
//                    for(NetNode<T> node2 : net2.dfs()){
//                        if(result.containsKey(node2)) continue;
//
//                        if(node1.getData() instanceof NetNodeInfo && node2.getData() instanceof NetNodeInfo) {
//                            double height1 = ((NetNodeInfo)node1.getData()).getHeight();
//                            double height2 = ((NetNodeInfo)node2.getData()).getHeight();
//
//                            double delta = Math.abs(height1 - height2);
//                            if(delta < closest) {
//                                closest = delta;
//                                nextNodeToSolve = node1;
//                                nextNode2ToSolve = node2;
//                            }
//                        }
//                    }
//                }

                double lowest1 = Double.MAX_VALUE;
                for(NetNode<T> node1 : net1.dfs()) {
                    if (result.containsKey(node1)) continue;
                    if (node1.getData() instanceof NetNodeInfo) {
                        double height1 = ((NetNodeInfo) node1.getData()).getHeight();
                        if (height1 < lowest1) {
                            lowest1 = height1;
                            nextNodeToSolve = node1;
                        }
                    }
                }

                double lowest2 = Double.MAX_VALUE;
                for(NetNode<T> node2 : net2.dfs()) {
                    if (result.containsValue(node2)) continue;
                    if (node2.getData() instanceof NetNodeInfo) {
                        double height2 = ((NetNodeInfo) node2.getData()).getHeight();
                        if (height2 < lowest2) {
                            lowest2 = height2;
                            nextNode2ToSolve = node2;
                        }
                    }
                }

                if(nextNodeToSolve == null)
                    return null;
            }

            result.put(nextNodeToSolve, nextNode2ToSolve);
        }

        return result;
    }

    SuperNetwork3(SNProblem problem) {
        this(problem.subnetworks);
    }

    SuperNetwork3(List<Tuple3<Network, String, Double>> subnetworks) {
        if(printDetails_)
            System.out.println(eps);

        subnetworks_ = new ArrayList<>();
        int index = 0;
        for(Tuple3<Network, String, Double> tuple3 : subnetworks) {
            Network<NetNodeInfo> newnet = Networks.readNetworkWithRootPop(Networks.getFullString(tuple3.Item1));
            if(newnet == null) continue;
            initNetHeights(newnet);
            subnetworks_.add(new NetworkWithInfo(newnet, tuple3.Item2, tuple3.Item3));

        }

    }

    public void ReduceTrinets() {
        Iterator<NetworkWithInfo> it = subnetworks_.iterator();
        while (it.hasNext()) {

            NetworkWithInfo netinfo = it.next();

            Set<String> taxa = new HashSet<>(netinfo.taxa);
            if(!taxa.contains(outgroup)) {
                it.remove();
            }
        }

        if(printDetails_) {
            System.out.println("Reduce number of trinets to " + subnetworks_.size());
        }
    }

    public void ReduceTrinets(Map<String, String> allele2species, String list_path) {
        Set<Set<String>> triset = new HashSet<>();

        try {
            BufferedReader in = new BufferedReader(new FileReader(list_path));
            String s;
            while((s = in.readLine()) != null) {
                List<String> cur_alleles = new ArrayList<>();
                String cur_allele = "";
                Set<String> cur_species = new HashSet<>();
                for(int i = 0 ; i < s.length() ; i++) {
                    if(s.charAt(i) == '[' || s.charAt(i) == ' ') {
                        continue;
                    } else if(s.charAt(i) == ',' || s.charAt(i) == ']') {
                        cur_alleles.add(cur_allele);
                        cur_allele = "";
                    } else {
                        cur_allele += s.charAt(i);
                    }
                }

                for(String allele : cur_alleles) {
                    cur_species.add(allele2species.get(allele));
                }

                if(cur_species.size() == 1) {
                    continue;
                } else if(cur_species.size() == 2) {
                    cur_species.add(outgroup);
                }

                triset.add(cur_species);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }


        subnetworks_reduceed_ = new ArrayList<>();
        Iterator<NetworkWithInfo> it = subnetworks_.iterator();
        while (it.hasNext()) {

            NetworkWithInfo netinfo = it.next();

            Set<String> taxa = new HashSet<>(netinfo.taxa);
            if(!triset.contains(taxa) /*&& !taxa.contains(outgroup)*/) {
                it.remove();
                subnetworks_reduceed_.add(netinfo);
            }
        }

        if(printDetails_) {
            System.out.println("Reduce number of trinets to " + subnetworks_.size());
        }



    }

    void CheckReducedTrinets() {
        Prepare();
        for(int i = 0 ; i < leafnames_.size() ; i++) {
            if(leafnames_.get(i).equals(outgroup)) continue;
            for(int j = i + 1 ; j < leafnames_.size() ; j++) {
                if(leafnames_.get(j).equals(outgroup)) continue;

                Set<String> curLeaves = new HashSet<>();
                curLeaves.add(outgroup);
                curLeaves.add(leafnames_.get(i));
                curLeaves.add(leafnames_.get(j));

//                if(printDetails_) {
//                    System.out.println("Check reduced trinets:");
//                    for (String leaf : curLeaves) {
//                        System.out.print(leaf + " ");
//                    }
//                    System.out.println();
//                }

                for(String targetLeafName : buildOrder) {
                    if(curLeaves.contains(targetLeafName)) continue;
                    //if(printDetails_) {
                    //    System.out.println(targetLeafName);
                    //}
                    boolean ok = false;
                    for(NetworkWithInfo netinfo : subnetworks_) {
                        if(!netinfo.trustTime) continue;
                        if(!netinfo.taxa.contains(targetLeafName)) continue;
                        List<String> leavesIntersection = new ArrayList<>();
                        leavesIntersection.addAll(curLeaves);
                        leavesIntersection.remove(outgroup);
                        leavesIntersection.retainAll(netinfo.taxa);
                        if(leavesIntersection.size() < 2) continue;

//                        if(targetLeafName.equals("F")) {
//                            System.out.println("!!!!!!!");
//
//                            for(String leaf : leavesIntersection) {
//                                System.out.print(leaf + " ");
//                            }
//                            System.out.println();
//                        }

                        ok = true;
                    }
                    if(!ok) {
                        String pickone = null;
                        int index = 0;

                        for(String leaf : curLeaves) {
                            if(!leaf.equals(outgroup)) {
                                pickone = leaf;
                                break;
                            }
                        }
                        List<String> requiredTaxa = new ArrayList<>();
                        requiredTaxa.add(outgroup);
                        requiredTaxa.add(targetLeafName);
                        requiredTaxa.add(pickone);

                        for(NetworkWithInfo netinfo : subnetworks_reduceed_) {
                            if(netinfo.taxa.contains(requiredTaxa.get(0)) && netinfo.taxa.contains(requiredTaxa.get(1)) &&netinfo.taxa.contains(requiredTaxa.get(2)) ) {
                                if(printDetails_) {
                                    System.out.println("Added: " + requiredTaxa.get(0) + " " + requiredTaxa.get(1) + " " + requiredTaxa.get(2));
                                }

                                subnetworks_reduceed_.remove(netinfo);
                                subnetworks_.add(netinfo);
                                break;
                            }
                        }
                    }

                    curLeaves.add(targetLeafName);
                }
            }
        }



        if(printDetails_) {
            System.out.println("After check, number of trinets is " + subnetworks_.size());
        }


    }

    void Prepare() {
        if(reconcileHeights)
            ReconcileSubnetHeights();

        leavesUnderReticulation = new HashSet<>();
        Map<String, Integer> underRetiCount = new HashMap<>();
        Map<String, Map<Integer, Integer>> retiBorderStats = new HashMap<>();
        Set<String> leafnames = new HashSet<>();
        for(NetworkWithInfo netinfo : subnetworks_) {
            //if(netinfo.taxa.contains("B") && netinfo.taxa.contains("D") &&netinfo.taxa.contains("L")) {
            //    System.out.println(netinfo.network);
            //}

            leafnames.addAll(netinfo.taxa);
            leavesUnderReticulation.addAll(Networks.getTaxaNamesUnderReticulation(netinfo.network));
            for(String name : Networks.getTaxaNamesUnderReticulation(netinfo.network)) {
                if(!underRetiCount.containsKey(name)) {
                    underRetiCount.put(name, 0);
                }
                underRetiCount.put(name, underRetiCount.get(name) + 1);

                if(!retiBorderStats.containsKey(name)) {
                    retiBorderStats.put(name, new TreeMap<>());
                }
                Set<NetNode> component1 = getComponent(netinfo, name);
                Set<NetNode> border1 = getBorder(component1);
                int borderSize = border1.size();
                if(!retiBorderStats.get(name).containsKey(borderSize)) {
                    retiBorderStats.get(name).put(borderSize, 0);
                }
                retiBorderStats.get(name).put(borderSize, retiBorderStats.get(name).get(borderSize) + 1);
            }
        }

        //int retiCountThreshold = underRetiCount.size() > 0 ? Collections.max(underRetiCount.values()) / 5 : 0;
        int retiCountThreshold = leafnames.size() - 2;

        for(NetworkWithInfo netinfo : subnetworks_) {
            for(String name : Networks.getTaxaNamesUnderReticulation(netinfo.network)) {
                if(underRetiCount.get(name) < retiCountThreshold / 2) {
                    netinfo.trustTime = false;
                    break;
                }

                Set<NetNode> component1 = getComponent(netinfo, name);
                Set<NetNode> border1 = getBorder(component1);
                int borderSize = border1.size();

                if(retiBorderStats.get(name).get(borderSize) < retiCountThreshold / 2) {
                    netinfo.trustTime = false;
                    break;
                }
            }
        }

        if(printDetails_) {
            System.out.println("Leaves under reticulations:");
            for(String name : underRetiCount.keySet()) {
                System.out.println(name + " " + underRetiCount.get(name));
            }
            System.out.println();
        }

        for(String name : underRetiCount.keySet()) {
            if(underRetiCount.get(name) < retiCountThreshold) {
                leavesUnderReticulation.remove(name);
            }
        }

        leavesUnderBubbleReticulation = new HashSet<>();

        Map<String, Integer> bubbleCount = new HashMap<>();
        for(String leaf : leavesUnderReticulation) {
            bubbleCount.put(leaf, 0);
            for(NetworkWithInfo netinfo : subnetworks_) {
                if(!netinfo.taxa.contains(leaf)) continue;
                Set<NetNode> ancestors = getAllAncestors(netinfo.network.findNode(leaf));
                for(NetNode ancestor : ancestors) {
                    if(ancestor.isNetworkNode()) {
                        if(getReticulationNodeDiameter(ancestor) == 3) {
                            bubbleCount.put(leaf, bubbleCount.get(leaf) + 1);
                        }
                    }
                }
            }

            //if(bubbleCount.get(leaf) >= retiCountThreshold) {
            //    leavesUnderBubbleReticulation.add(leaf);
            //}
            if(underRetiCount.get(leaf) == leafnames.size() - 2) {
                leavesUnderBubbleReticulation.add(leaf);
            }
        }

        leafnames_ = new ArrayList<>();
        leafnames_.addAll(leafnames);
        Collections.sort(leafnames_);

        // Compute build order according to reticulations
        buildOrder = new ArrayList<>();
        Map<String, Integer> indeg = new HashMap<>();
        Map<String, Set<String>> buildBefore = new HashMap<>();

        for(int i = 0 ; i < leafnames_.size() ; i++) {
            indeg.put(leafnames_.get(i), 0);
            buildBefore.put(leafnames_.get(i), new HashSet<>());
        }

        for(NetworkWithInfo netinfo : subnetworks_) {
            if(!netinfo.trustTime) continue;
            for(String r : netinfo.dependencies.keySet()) {
                for(String s : netinfo.dependencies.get(r)) {
                    if(s.equals(r)) continue;
                    if(!buildBefore.containsKey(s))
                        buildBefore.put(s, new HashSet<>());
                    buildBefore.get(s).add(r);
                }
            }
        }

        for(String s : buildBefore.keySet()) {
            for(String r : buildBefore.get(s)) {
                indeg.put(r, indeg.get(r) + 1);
            }
        }

        for(String leaf : leavesUnderBubbleReticulation) {
            indeg.put(leaf, indeg.get(leaf) + 10);
        }

        // Topological sorting
        boolean finished = false;
        while(!finished) {
            String next = null;
            for(String s : indeg.keySet()) {
                /*if(indeg.get(s) == 0) {
                    if(next == null) next = s;
                    else if(leavesUnderReticulation.contains(s))  {
                        next = s;
                        break;
                    }
                } else {
                    if(next == null) next = s;
                    else if(indeg.get(next) > indeg.get(s)) next = s;
                }*/
                if(next == null) next = s;
                else if(indeg.get(next) > indeg.get(s) || (indeg.get(next).equals(indeg.get(s)) && (underRetiCount.getOrDefault(next, 0) < underRetiCount.getOrDefault(s, 0)))) next = s;
            }
            if(next == null) {
                if(indeg.size() == 0) finished = true;
                break;
            }
            buildOrder.add(next);
            for(String r : buildBefore.get(next)) {
                if(indeg.containsKey(r))
                    indeg.put(r, indeg.get(r) - 1);
            }
            buildBefore.remove(next);
            indeg.remove(next);
        }

        if(!finished) {
            for(int i = 0 ; i < leafnames_.size() ; i++) {
                if(!buildOrder.contains(leafnames_.get(i))) {
                    buildOrder.add(leafnames_.get(i));
                }
            }
        }

        if(printDetails_) {
            System.out.println("Build Order:");
            for (String s : buildOrder) {
                System.out.print(s + " ");
            }
            System.out.println();
        }

        // Compute extended height matrix
        ComputeExtendedHeightMatrix();


        // Recompute build order
        Mean popsizemean = new Mean();
        for(NetworkWithInfo netinfo : subnetworks_) {
            if (!netinfo.trustTime) continue;
            popsizemean.increment(netinfo.network.getRoot().getRootPopSize());
        }
        popsize = popsizemean.getResult();

        // Backup subnets
        subnetworks_removed_ = new ArrayList<>();
//        for(NetworkWithInfo netinfo : subnetworks_) {
//            if (!netinfo.trustTime) continue;
//            subnetworks_backup_.add(Networks.getFullString(netinfo.network));
//        }
    }

    void ReconcileSubnetHeights() {
        Map<NetNode, Mean> heights = new HashMap<>();

        for(int ii = 0 ; ii < subnetworks_.size() ; ii++) {
            NetworkWithInfo curNet = subnetworks_.get(ii);
            for(int i = 0 ; i < curNet.taxa.size() ; i++) {
                for(int j = i + 1 ; j < curNet.taxa.size() ; j++) {
                    List<String> selected = new ArrayList<>();
                    selected.add(curNet.taxa.get(i));
                    selected.add(curNet.taxa.get(j));

                    Tuple<Network, Map<NetNode, NetNode>> tuple1 = getSubNetwork(curNet.network, selected, true);

                    for(NetworkWithInfo nextNet : subnetworks_) {
                        if(nextNet == curNet) continue;

                        List<String> interscetion = new ArrayList<>(selected);
                        interscetion.retainAll(nextNet.taxa);
                        if(interscetion.size() != 2) continue;

                        Tuple<Network, Map<NetNode, NetNode>> tuple2 = getSubNetwork(nextNet.network, selected, true);

                        if(Networks.hasTheSameTopology(tuple1.Item1, tuple2.Item1)) {
                            Map<NetNode, NetNode> sc2sn = mapTwoNetworks(tuple1.Item1, tuple2.Item1);
                            Map<NetNode, NetNode> sc2c = tuple1.Item2;
                            Map<NetNode, NetNode> sn2n = tuple2.Item2;

                            if(sc2sn == null) continue;
                            for(NetNode scNode : sc2sn.keySet()) {
                                NetNode snNode = sc2sn.get(scNode);
                                NetNode cNode = sc2c.get(scNode);
                                NetNode nNode = sn2n.get(snNode);

                                if(!heights.containsKey(cNode)) {
                                    heights.put(cNode, new Mean());
                                }
                                heights.get(cNode).increment(((NetNodeInfo)nNode.getData()).getHeight());
                            }
                        }
                    }
                }
            }
        }

        for(NetNode node : heights.keySet()) {
            ((NetNodeInfo) node.getData()).setHeight(heights.get(node).getResult());
        }

        for(NetworkWithInfo net : subnetworks_) {
            resetBranchLengths(net.network);
            net.ehm = getEHMFromNetwork(net.network);
        }


    }

    void ComputeExtendedHeightMatrix(){
        ehm = new HashMap<>();
        Map<Tuple<String, String>, List<Mean>> ehm_mean = new HashMap<>();
        for(NetworkWithInfo netinfo : subnetworks_) {
            if(!netinfo.trustTime) continue;
            for(Tuple<String, String> tuple : netinfo.ehm.keySet()) {
                if(tuple.Item1.compareTo(tuple.Item2) < 0) {
                    if(!ehm.containsKey(tuple)) {
                        ehm.put(tuple, new ArrayList<>());
                    }
                    if(!ehm_mean.containsKey(tuple)) {
                        ehm_mean.put(tuple, new ArrayList<>());
                    }

                    if(netinfo.ehm.get(tuple).size() > ehm.get(tuple).size() || (netinfo.ehm.get(tuple).size() == ehm.get(tuple).size() && netinfo.ehm.get(tuple).get(0) < ehm.get(tuple).get(0))) {
                        ehm.put(tuple, new ArrayList<>(netinfo.ehm.get(tuple)));
                    }

                    if(netinfo.ehm.get(tuple).size() > ehm_mean.get(tuple).size() ) {
                        ehm_mean.put(tuple, new ArrayList<>());
                        for(Double d : netinfo.ehm.get(tuple)) {
                            Mean m = new Mean();
                            m.increment(d);
                            ehm_mean.get(tuple).add(m);
                        }
                    } else if(netinfo.ehm.get(tuple).size() == ehm_mean.get(tuple).size()) {
                        for(int i = 0 ; i < netinfo.ehm.get(tuple).size() ; i++) {
                            ehm_mean.get(tuple).get(i).increment(netinfo.ehm.get(tuple).get(i));
                        }
                    }

                    /*for(double d : netinfo.ehm.get(tuple)) {
                        boolean found = false;
                        for(double d1 : ehm.get(tuple)) {
                            if(Math.abs(d - d1) < eps) {
                                found = true;
                                break;
                            }
                        }
                        if(!found)
                            ehm.get(tuple).add(d);
                    }*/
                }
            }
        }

        for(NetworkWithInfo netinfo : subnetworks_) {
            for (Tuple<String, String> tuple : netinfo.ehm.keySet()) {
                if(!ehm.containsKey(tuple)) {
                    ehm.put(tuple, new ArrayList<>(netinfo.ehm.get(tuple)));
                    break;
                }
            }
        }

//        for(Tuple<String, String> tuple : ehm_mean.keySet()) {
//            ehm.put(tuple, new ArrayList<>());
//            for(int i = 0 ; i < ehm_mean.get(tuple).size() ; i++) {
//                ehm.get(tuple).add(ehm_mean.get(tuple).get(i).getResult());
//            }
//        }

        for(Tuple<String, String> tuple : ehm.keySet()) {
            Collections.sort(ehm.get(tuple));
        }


        // Compute avg pop size

    }

    List<String> getTaxaNames() {
        return leafnames_;
    }

    double computeScore(Network network) {
        double score = 0.0;
        List<String> curleaves = new ArrayList<>();
        for(Object leafObj : network.getLeaves()) {
            curleaves.add(((NetNode)leafObj).getName());
        }
        for (NetworkWithInfo netinfo : subnetworks_) {
            if (netinfo.dirty) {
                netinfo.network = Networks.readNetwork(netinfo.backup);
                netinfo.dirty = false;
            }
            List<String> leavesIntersection = new ArrayList<>();
            leavesIntersection.addAll(curleaves);
            leavesIntersection.retainAll(netinfo.taxa);
            if(leavesIntersection.size() < 2) continue;
            score += compareTwoSubNetwork(network, netinfo.network, leavesIntersection);
        }
        return score;
    }

    Set<NetNode> getComponent(NetworkWithInfo net, String targetLeafName) {

        Set<NetNode> ancestorsOfOthers = new HashSet<>();
        Set<NetNode> ancestorsOfTarget = new HashSet<>();
        List<NetNode> otherLeaves = new ArrayList<>();
        for(NetNode leaf : net.network.getLeaves()) {
            if(!leaf.getName().equals(targetLeafName)) {
                ancestorsOfOthers.addAll(getAllAncestors(leaf));
                otherLeaves.add(leaf);
            } else {
                ancestorsOfTarget.addAll(getAllAncestors(leaf));
            }
        }

        ancestorsOfTarget.removeAll(ancestorsOfOthers);
        ancestorsOfTarget.add(net.network.findNode(targetLeafName));

        return ancestorsOfTarget;
    }

    Set<NetNode> cutBorder(Set<NetNode> component) {
        Set<NetNode> border = getBorder(component);

        // Cut border down
        for(NetNode node : border) {
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                if(!component.contains(parent)) parent.removeChild(node);
                else throw new RuntimeException("Should not come here");
            }
            for(Object childObj : node.getChildren()) {
                NetNode child = (NetNode) childObj;
                if(!component.contains(child)) node.removeChild(child);
            }
        }

        return border;
    }

    Set<NetNode> getBorder(Set<NetNode> component) {
        Set<NetNode> border = new HashSet<>();
        // Find border of component
        for(NetNode node : component) {
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                if(!component.contains(parent)) {
                    border.add(parent);
                }
            }
        }

        return border;
    }

    static Network duplicateNetworkTopDown(Network net) {
        Map<BniNetNode, NetNode> new2old = new HashMap<>();
        Map<NetNode, BniNetNode> old2new = new HashMap<>();
        Queue<NetNode> queue = new LinkedList<>();
        queue.offer(net.getRoot());
        while(!queue.isEmpty()) {
            NetNode oldNode = queue.poll();
            BniNetNode newNode = new BniNetNode();
            if(oldNode.getData() != null) newNode.setData(new NetNodeInfo(((NetNodeInfo)oldNode.getData()).getHeight()));
            if(oldNode.isLeaf()) newNode.setName(oldNode.getName());
            new2old.put(newNode, oldNode);
            old2new.put(oldNode, newNode);


            for(Object childObj : oldNode.getChildren()) {
                NetNode child = (NetNode) childObj;
                if(!old2new.containsKey(child)) queue.offer(child);
            }
        }

        for(NetNode oldNode : old2new.keySet()) {
            BniNetNode newNode = old2new.get(oldNode);
            for(Object childObj : oldNode.getChildren()) {
                NetNode oldChild = (NetNode) childObj;
                BniNetNode newChild = old2new.get(oldChild);

                newNode.adoptChild(newChild, oldChild.getParentDistance(oldNode));
                newChild.setParentSupport(newNode, oldChild.getParentSupport(oldNode));
                newChild.setParentProbability(newNode, oldChild.getParentProbability(oldNode));
                newChild.setParentProbability(newNode, NetNode.NO_PROBABILITY);
            }
        }

        BniNetwork newnet = new BniNetwork(old2new.get(net.getRoot()));
        Networks.removeBinaryNodes(newnet);
        if(!resetBranchLengths(newnet)) return null;
        return newnet;
    }

    public static Iterable<List> Permutate(List list) {
        return new Iterable<List>()
        {
            @Override
            public Iterator<List> iterator()
            {
                return new Iterator<List>()
                {

                    int index[] = null;
                    int size = list.size();

                    @Override
                    public void remove()
                    {
                        throw new UnsupportedOperationException();
                    }


                    @Override
                    public boolean hasNext()
                    {
                        if(index == null) return true;
                        for(int i = 0, j = size - 1 ; i < index.length ; i++, j--) {
                            if(index[j] != i) return true;
                        }
                        return false;
                    }

                    @Override
                    public List next()
                    {
                        if (index == null) {
                            index = new int[size];
                            for(int i = 0 ; i < size; i++)
                                index[i] = i;
                        }
                        else
                        {
                            tryAdvance();
                        }

                        List temp = new ArrayList();
                        for(int i = 0 ; i < size ; i++)
                            temp.add(list.get(index[i]));

                        return temp;
                    }

                    private void tryAdvance()
                    {
                        int j = size - 1;
                        while(j > 0 && index[j - 1] >= index[j]) {
                            j--;
                        }
                        if(j == 0) throw new RuntimeException("Already get to end!");

                        int k = index[j - 1];
                        for(int i = size - 1 ; i >= j ; i--) {
                            if (index[i] > k) {
                                index[j - 1] = index[i];
                                index[i] = k;
                                break;
                            }
                        }
                        Arrays.sort(index, j, size);
                    }
                };
            }
        };

    }

    void enumerateDFS(Network net, Set<NetNode> component, List<NetNode> border, int index, List<Network> results, double prevScore) {
        Network dupnet = duplicateNetworkTopDown(net);

        if(dupnet == null) return;

        double curScore = 0;//index > 0 ? computeScore(dupnet) : Double.MAX_VALUE;
        //if(curScore > prevScore) return;

        if(index  > 0 && dupnet != null)
            results.add(dupnet);

        if(index == border.size()) {
            /*if (Networks.hasCycle(net)) {
                return;
            }
            if (!Networks.isDisconnectedNetwork(net, null)) {
                return;
            }
            results.add(Networks.readNetwork(net.toString()));*/
            return;
        }

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = Networks.getAllEdges(net);
        NetNode toInsert = border.get(index);
        List<Double> scores = new ArrayList<>();
        double bestScore = Double.MAX_VALUE;
        for(Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge : edges) {
            NetNode<NetNodeInfo> child = edge.Item1;
            NetNode<NetNodeInfo> parent = edge.Item2;
            if(component.contains(child) || component.contains(parent)) {
                scores.add(Double.MAX_VALUE);
                continue;
            }
            double distanceBackup = child.getParentDistance(parent);
            parent.removeChild(child);
            parent.adoptChild(toInsert, NetNode.NO_DISTANCE);
            toInsert.adoptChild(child, NetNode.NO_DISTANCE);

            //Network dupnet0 = duplicateNetworkTopDown(net);
            //curScore = computeScore(dupnet0);
            //scores.add(curScore);
            //if(curScore < bestScore) bestScore = curScore;
            enumerateDFS(net, component, border, index + 1, results, curScore);

            parent.removeChild(toInsert);
            toInsert.removeChild(child);
            parent.adoptChild(child, distanceBackup);
        }


        //System.out.println();
    }

    Map<Tuple<String, String>, List<Double>> getEHMFromNetwork(Network<NetNodeInfo> net) {
        Map<Tuple<String, String>, List<Double>> ehm = new HashMap<>();
        List<String> taxa = new ArrayList<>();
        for(NetNode leaf : net.getLeaves()) {
            taxa.add(leaf.getName());
        }
        Collections.sort(taxa);

        for(int i = 0 ; i < taxa.size() ; i++) {
            for(int j = i + 1 ; j < taxa.size() ; j++) {
                List<String> selected = new ArrayList<>();
                selected.add(taxa.get(i));
                selected.add(taxa.get(j));
                Network subnet = getSubNetwork(net, selected, true).Item1;
                initNetHeights(subnet);
                Map<Tuple<String, String>, List<Double>> singleEHM = getEHMFromNetworkHelper(subnet);
                Tuple<String, String> tuple = new Tuple<>(taxa.get(i), taxa.get(j));
                ehm.put(tuple, singleEHM.get(tuple));
            }
        }
        return ehm;
    }

    Map<Tuple<String, String>, List<Double>> getEHMFromNetworkHelper(Network<NetNodeInfo> net) {
        Map<Tuple<String, String>, List<Double>> ehm;
        List<String> taxa = new ArrayList<>();
        for(NetNode leaf : net.getLeaves()) {
            taxa.add(leaf.getName());
        }
        Collections.sort(taxa);

        ehm = new HashMap<>();
        for(int i = 0 ; i < taxa.size() ; i++) {
            for(int j = i + 1 ; j < taxa.size() ; j++) {
                ehm.put(new Tuple<>(taxa.get(i), taxa.get(j)), new ArrayList<>());
                Set<NetNode> MCRAs = getMCRAs(net, taxa.get(i), taxa.get(j));
                for(NetNode node : MCRAs) {
                    if(node.isNetworkNode()) continue;
                    double value = ((NetNodeInfo) node.getData()).getHeight();
                    ehm.get(new Tuple<>(taxa.get(i), taxa.get(j))).add(value);
                }
                Collections.sort(ehm.get(new Tuple<>(taxa.get(i), taxa.get(j))));
            }
        }

        return ehm;
    }

    Map<String, List<Double>> getTargetedHeightsFromEHM(Map<Tuple<String, String>, List<Double>> ehm, String target) {
        Map<String, List<Double>> heights = new HashMap<>();
        for(Tuple<String, String> tuple : ehm.keySet()) {
            String x = null;
            if(tuple.Item1.equals(tuple.Item2)) continue;
            if(tuple.Item1.equals(target)) x = tuple.Item2;
            else if(tuple.Item2.equals(target)) x = tuple.Item1;
            else continue;

            if(x == null) throw new RuntimeException("Should not be here!");
            heights.put(x, new ArrayList<>(ehm.get(tuple)));
        }
        return heights;
    }

    public static double getUpperBoundDfs(NetNode<NetNodeInfo> node) {
        if(node.getData() != null) {
            return node.getData().getHeight();
        }

        double bound = Double.MAX_VALUE;
        for(NetNode<NetNodeInfo> parent : node.getParents()) {
            bound = Math.min(bound, getUpperBoundDfs(parent));
        }
        return bound;
    }

    public static double getLowerBoundDfs(NetNode<NetNodeInfo> node) {
        if(node.getData() != null) {
            return node.getData().getHeight();
        }

        double bound = 0.0;
        for(NetNode<NetNodeInfo> child : node.getChildren()) {
            bound = Math.max(bound, getLowerBoundDfs(child));
        }
        return bound;
    }

    public double[] getLowerAndUpperBound(NetNode<NetNodeInfo> node) {
        double[] bounds = new double[] {Double.MIN_VALUE, Double.MAX_VALUE};
        for(NetNode<NetNodeInfo> child: node.getChildren()) {
            bounds[0] = Math.max(bounds[0], child.getData().getHeight());
        }
        for(NetNode<NetNodeInfo> par: node.getParents()) {
            bounds[1] = Math.min(bounds[1], par.getData().getHeight());
        }
        return bounds;
    }

    static boolean resetBranchLengths(Network<NetNodeInfo> net) {
        for(NetNode node : net.dfs()) {
            if (node.getData() == null) continue;
            for (Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                if (parent.getData() != null) {
                    double newlength = ((NetNodeInfo) parent.getData()).getHeight() - ((NetNodeInfo) node.getData()).getHeight();
                    if(newlength < 0) return false;
                    node.setParentDistance(parent, newlength);
                }
            }
        }
        return true;
    }

    void resetHeightEHM(Network<NetNodeInfo> net) {

        List<NetNode<NetNodeInfo>> leaves = new ArrayList<>();
        for(NetNode<NetNodeInfo> leaf : net.getLeaves()) {
            leaves.add(leaf);
        }

        for(NetNode<NetNodeInfo> node : net.dfs()) {
            if(!node.isLeaf() && !node.isNetworkNode()) {
                for(NetNode<NetNodeInfo> leaf1 : leaves) {
                    for(NetNode<NetNodeInfo> leaf2 : leaves) {
                        if(leaf1 == leaf2) continue;
                        List<NetNode> nodepair = new ArrayList<>();
                        nodepair.add(leaf1);
                        nodepair.add(leaf2);
                        NetNode<NetNodeInfo> mcra = getMCRA(nodepair);
                        if(mcra == node) {
                            List<Double> heightList = leaf1.getName().compareTo(leaf2.getName()) < 0 ? ehm.get(new Tuple<>(leaf1.getName(), leaf2.getName())) : ehm.get(new Tuple<>(leaf2.getName(), leaf1.getName()));
                            if(heightList == null) continue;
                            if(heightList.get(0) < node.getData().getHeight()) {
                                double lowerbound = getLowerBoundDfs(node);
                                if(heightList.get(0) < lowerbound) {
                                    continue;
                                }
                                node.getData().setHeight(heightList.get(0));
                            }
                        }
                    }
                }
            }
        }


    }

    void resetHeightsInComponent(Set<NetNode> component, List<NetNode<NetNodeInfo>> border) {


        Set<NetNode> heightGood = new HashSet<>();
        for(NetNode node : border) {
            if(node.getData() != null)
                heightGood.add(node);
        }

        if(!trustReticulationTime) {
            Set<NetNode> heightToSet = new HashSet<>(component);
            for (NetNode node : component) {
                node.setData(null);
                for (Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    node.setParentDistance(parent, NetNode.NO_DISTANCE);
                }
            }

            while (heightToSet.size() > 0) {
                boolean updated = false;
                for (NetNode node : heightToSet) {
                    node.setData(null);

                    for (Object parentObj : node.getParents()) {
                        NetNode parent = (NetNode) parentObj;
                        if (heightGood.contains(parent)) {
                            if (!node.isLeaf()) {
                                //double desired = ((NetNodeInfo) parent.getData()).getHeight() * 0.95;
                                //if(node.getData() == null || ((NetNodeInfo) node.getData()).getHeight() > desired)
                                //    node.setData(new NetNodeInfo(desired));
                            } else
                                node.setData(new NetNodeInfo(0.0));
                        }
                    }
                    if (node.getData() != null) {
                        for (Object parentObj : node.getParents()) {
                            NetNode parent = (NetNode) parentObj;
                            if (parent.getData() != null)
                                node.setParentDistance(parent, ((NetNodeInfo) parent.getData()).getHeight() - ((NetNodeInfo) node.getData()).getHeight());
                        }
                        heightGood.add(node);
                        heightToSet.remove(node);
                        updated = true;
                        break;
                    }
                }
                if (!updated) break;
            }
        }

        /*
        //Reset branch lengths
        for(NetNode node : component) {
            if(node.getData() == null) continue;
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                if(parent.getData() != null)
                    node.setParentDistance(parent, ((NetNodeInfo)parent.getData()).getHeight() - ((NetNodeInfo)node.getData()).getHeight());
            }
        }*/
    }

    void enumerateDFS_EHM(Network net, Set<NetNode> component, List<NetNode<NetNodeInfo>> border, String targetLeaf, List<Tuple<Double, String>> prevHeightsToResolve, int index, List<Network> results) {
        Network dupnet = duplicateNetworkTopDown(net);
        if(dupnet == null) return;

        if(index  > 0 && dupnet != null)
            results.add(dupnet);

        if(index == border.size()) return;

        List<Tuple<Double, String>> nextHeightsToResolve = new ArrayList<>();
        // Clean up heights to resolve
        if(index > 0) {
            Map<Tuple<String, String>, List<Double>> curehm = getEHMFromNetwork(net);
            Map<String, List<Double>> curHeights = getTargetedHeightsFromEHM(curehm, targetLeaf);

            for(Tuple<Double, String> tuple : prevHeightsToResolve) {
                boolean resolved = false;
                for(Double d : curHeights.get(tuple.Item2)) {
                    if(Math.abs(tuple.Item1 - d) < eps) {
                        resolved = true;
                        break;
                    }
                }
                if(!resolved) nextHeightsToResolve.add(tuple);
            }
        } else {
            nextHeightsToResolve.addAll(prevHeightsToResolve);
        }

        if(nextHeightsToResolve.size() == 0) {
            // Hooray!
            return;
        }
        Tuple<Double, String> currentHeightToResolve = nextHeightsToResolve.get(0);

//        if(targetLeaf.equals("L")) {
//            for(Tuple<Double, String> tuple : nextHeightsToResolve) {
//                System.out.println(tuple.Item1 + " " + tuple.Item2);
//            }
//            System.out.println();
//        }

        Set<NetNode> nodesToCheck = getAllAncestors(net.findNode(currentHeightToResolve.Item2));
        nodesToCheck.add(net.findNode(currentHeightToResolve.Item2));

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edgesToTry = new ArrayList<>();
        for(NetNode node : nodesToCheck) {
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;

                if(parent.getData() != null && node.getData() != null) {
                    double parentHeight = ((NetNodeInfo) parent.getData()).getHeight();
                    double curHeight = ((NetNodeInfo) node.getData()).getHeight();

                    if (parentHeight > currentHeightToResolve.Item1 && curHeight < currentHeightToResolve.Item1) {
                        edgesToTry.add(new Tuple<>(node, parent));
                    }
                } else {
                    if (getUpperBoundDfs(parent) > currentHeightToResolve.Item1 && getLowerBoundDfs(node) < currentHeightToResolve.Item1) {
                        edgesToTry.add(new Tuple<>(node, parent));
                    }
                }
            }
        }

        NetNode<NetNodeInfo> toInsert = border.get(index);
        NetNode<NetNodeInfo> oneChild = toInsert.getChildren().iterator().next();
        List<NetNode<NetNodeInfo>> oneChildChildren = new ArrayList<>();
        for(NetNode<NetNodeInfo> child : oneChild.getChildren()) oneChildChildren.add(child);
        //System.out.println(toInsert.getChildren().iterator().next().getChildCount());

        for(NetNode<NetNodeInfo> childToKeep : oneChildChildren) {
            NetNode<NetNodeInfo> childToRemove = null;
            for(NetNode<NetNodeInfo> child : oneChild.getChildren()) {
                if(child != childToKeep) {
                    childToRemove = child;
                    break;
                }
            }

            double childToRemoveDistBackup = NetNode.NO_DISTANCE;
            if(childToRemove != null) {
                childToRemoveDistBackup = childToRemove.getParentDistance(oneChild);
                oneChild.removeChild(childToRemove);


                for(Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge : edgesToTry) {
                    NetNode<NetNodeInfo> child = edge.Item1;
                    NetNode<NetNodeInfo> parent = edge.Item2;
                    if(component.contains(child) || component.contains(parent)) {
                        //scores.add(Double.MAX_VALUE);
                        continue;
                    }
                    double distanceBackup = child.getParentDistance(parent);
                    parent.removeChild(child);
                    toInsert.setData(new NetNodeInfo(currentHeightToResolve.Item1));
                    parent.adoptChild(toInsert, parent.getData() != null ? parent.getData().getHeight() - toInsert.getData().getHeight() : NetNode.NO_DISTANCE);
                    toInsert.adoptChild(child, child.getData()!= null ? toInsert.getData().getHeight() - child.getData().getHeight() : NetNode.NO_DISTANCE);
                    resetHeightsInComponent(component, border);

                    //Network dupnet0 = duplicateNetworkTopDown(net);
                    //curScore = computeScore(dupnet0);
                    //scores.add(curScore);
                    //if(curScore < bestScore) bestScore = curScore;
                    enumerateDFS_EHM(net, component, border, targetLeaf, nextHeightsToResolve, index + 1, results);

                    parent.removeChild(toInsert);
                    toInsert.removeChild(child);
                    parent.adoptChild(child, distanceBackup);
                    toInsert.setData(null);
                }


                oneChild.adoptChild(childToRemove, childToRemoveDistBackup);
            }
        }

        for(Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge : edgesToTry) {
            NetNode<NetNodeInfo> child = edge.Item1;
            NetNode<NetNodeInfo> parent = edge.Item2;
            if(component.contains(child) || component.contains(parent)) {
                //scores.add(Double.MAX_VALUE);
                continue;
            }
            double distanceBackup = child.getParentDistance(parent);
            parent.removeChild(child);
            toInsert.setData(new NetNodeInfo(currentHeightToResolve.Item1));
            parent.adoptChild(toInsert, parent.getData() != null ? parent.getData().getHeight() - toInsert.getData().getHeight() : NetNode.NO_DISTANCE);
            toInsert.adoptChild(child, child.getData()!= null ? toInsert.getData().getHeight() - child.getData().getHeight() : NetNode.NO_DISTANCE);
            resetHeightsInComponent(component, border);

            //Network dupnet0 = duplicateNetworkTopDown(net);
            //curScore = computeScore(dupnet0);
            //scores.add(curScore);
            //if(curScore < bestScore) bestScore = curScore;
            enumerateDFS_EHM(net, component, border, targetLeaf, nextHeightsToResolve, index + 1, results);

            parent.removeChild(toInsert);
            toInsert.removeChild(child);
            parent.adoptChild(child, distanceBackup);
            toInsert.setData(null);
        }


    }

    List<Network> enumerate_EHM(Network backbone, Set<NetNode> component, Set<NetNode> border, Set<String> curleaves, String targetName) {
        List<NetNode<NetNodeInfo>> borderList = new ArrayList<>();
        for(NetNode b : border) borderList.add((NetNode<NetNodeInfo>)b);

        Collections.sort(borderList, (NetNode<NetNodeInfo> a, NetNode<NetNodeInfo> b)->Double.compare(a.getData().getHeight(), b.getData().getHeight()));

        for(NetNode node : border) {
            node.setData(null);
        }
        List<Network> results = new ArrayList<>();

        for(Object nodeObj : backbone.getLeaves()) {
            NetNode node = (NetNode) nodeObj;
            node.setData(new NetNodeInfo(0.0));
        }

        Map<String, List<Double>> heightToResolve = new HashMap<>();
        List<Tuple<Double, String>> heightVSname = new ArrayList<>();
        for(String leaf : curleaves) {
            List<Double> heightList = leaf.compareTo(targetName) < 0 ? ehm.get(new Tuple<>(leaf, targetName)) : ehm.get(new Tuple<>(targetName, leaf));
            if(heightList == null) {
                //System.out.println("!!! Empty height list !!! " + leaf + " " + targetName);
                continue;
            }
            List<Double> newList = new ArrayList<>(heightList);
            heightToResolve.put(leaf, newList);
            for(double d : newList) {
                heightVSname.add(new Tuple<>(d, leaf));
            }
        }

        Collections.sort(heightVSname, (Tuple<Double, String> a, Tuple<Double, String> b)->Double.compare(a.Item1, b.Item1));

        //List<NetNode<NetNodeInfo>> borderListToTry = borderList;
        for(List<NetNode<NetNodeInfo>> borderListToTry : Permutate(borderList)) {
            enumerateDFS_EHM(backbone, component, borderListToTry, targetName, heightVSname, 0, results);
        }
        return results;
    }

    List<Network> enumerate(Network backbone, Set<NetNode> component, Set<NetNode> border, String targetName) {
        List<NetNode> borderList = new ArrayList<>();
        borderList.addAll(border);

        for(NetNode node : border) {
            node.setData(null);
        }

        for(NetNode node : component) {
            node.setData(null);
        }

        List<Network> results = new ArrayList<>();
        for(List<NetNode> borderListToTry : Permutate(borderList)) {
            enumerateDFS(backbone, component, borderListToTry,0, results, Double.MAX_VALUE);
        }

        return results;
    }

    static double compareTwoNetworkSymmetric(Network net1, Network net2) {
        NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
        return metric.computeDistanceBetweenTwoNetworks(net1, net2);
    }

    static double compareTwoNetworkAsymmetric(Network net1, Network net2) {
        if(Pipeline.isBackboneOf(net1, net2)) {
            return 0.0;
        }

        NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
        return metric.computeDistanceBetweenTwoNetworks(net1, net2);
    }

    static double compareTwoSubNetwork(Network net1, Network net2, List<String> selectedLeaves) {
        Network subnet1 = getSubNetwork(net1, selectedLeaves, true).Item1;
        Network subnet2 = getSubNetwork(net2, selectedLeaves, true).Item1;
        return compareTwoNetworkSymmetric(subnet1, subnet2);
    }

    static double compareTwoSubNetworkA(Network net1, Network net2, List<String> selectedLeaves) {
        Network subnet1 = getSubNetwork(net1, selectedLeaves, true).Item1;
        Network subnet2 = getSubNetwork(net2, selectedLeaves, true).Item1;
        return compareTwoNetworkAsymmetric(subnet1, subnet2);
    }

    Network ElectBackbone(Set<String> curleaves) {
        Network backbone = null;
        int skip = 0;
        double bestScore = Double.MAX_VALUE;

        leavesUnderReticulation.remove(outgroup);
        for(NetworkWithInfo candidate : subnetworks_ ) {
            boolean flag = true;
            for(String taxon : candidate.taxa) {
                if(leavesUnderReticulation.contains(taxon) /*&& netinfo.network.getReticulationCount() == 0*/) {
                    flag = false;
                    break;
                }
            }

            for(Tuple<Network, Map<NetNode, NetNode>> binetT : candidate.binets) {
                List<String> bileaves = new ArrayList<>();
                for(Object leafObj : binetT.Item1.getLeaves()) {
                    NetNode leaf = (NetNode) leafObj;
                    bileaves.add(leaf.getName());
                }

                if(!bileaves.contains(outgroup)) continue;

                double score = 0.0;
                for (NetworkWithInfo netinfo : subnetworks_) {
                    if(!netinfo.trustTime) continue;

                    if (netinfo.dirty) {
                        netinfo.network = Networks.readNetwork(netinfo.backup);
                        initNetHeights(netinfo.network);
                        netinfo.dirty = false;
                    }
                    List<String> leavesIntersection = new ArrayList<>();
                    leavesIntersection.addAll(bileaves);
                    leavesIntersection.retainAll(netinfo.taxa);
                    if(leavesIntersection.size() < 2) continue;
                    score += compareTwoSubNetwork(binetT.Item1, netinfo.network, leavesIntersection) * netinfo.percentage;
                }

                if(score < bestScore) {
                    backbone = binetT.Item1.clone();
                    bestScore = score;
                    curleaves.clear();
                    curleaves.addAll(bileaves);
                }
            }

        }

        return backbone;
    }

    Network ElectBackbone3(Set<String> curleaves) {
        Network backbone = null;
        int skip = 0;
        double bestScore = Double.MAX_VALUE;

        leavesUnderReticulation.remove(outgroup);
        for(NetworkWithInfo candidate : subnetworks_ ) {
            double score = 0.0;

            for(String taxon : candidate.taxa) {
                if(leavesUnderReticulation.contains(taxon) /*&& netinfo.network.getReticulationCount() == 0*/) {
                    score += 1.0;
                    break;
                }
            }

            for (NetworkWithInfo netinfo : subnetworks_) {
                if(!netinfo.trustTime) continue;

                if (netinfo.dirty) {
                    netinfo.network = Networks.readNetwork(netinfo.backup);
                    initNetHeights(netinfo.network);
                    netinfo.dirty = false;
                }
                List<String> leavesIntersection = new ArrayList<>();
                leavesIntersection.addAll(candidate.taxa);
                leavesIntersection.retainAll(netinfo.taxa);
                if(leavesIntersection.size() < 2) continue;
                score += compareTwoSubNetwork(candidate.network, netinfo.network, leavesIntersection) * netinfo.percentage;
            }

            if(score < bestScore) {
                backbone = candidate.network.clone();
                bestScore = score;
                curleaves.clear();
                curleaves.addAll(candidate.taxa);
            }

        }

        if(printDetails_) {
            System.out.println("Backbone score: " + bestScore);
        }
        initNetHeights(backbone);
        return backbone;
    }

    double ComputeBestScore(List<Network> networksToTry, List<Network> bestNetworks) {
        double bestScore = Double.MAX_VALUE;
        for(Network netToTry : networksToTry) {
            double score = 0.0;
            List<Network> backbonesToTry = null;

            //if(netToTry.getReticulationCount() < 5)
            //    backbonesToTry = Pipeline.getAllBackboneNets(netToTry, 5);
            //else
                backbonesToTry = new ArrayList<>();

            //if(netToTry.getReticulationCount() == 5) {
            //    System.out.println("!!!!!!!!! " + backbonesToTry.size());
            //}

            backbonesToTry.add(netToTry.clone());

            for (NetworkWithInfo netinfo : subnetworks_) {
                if(!netinfo.trustTime) continue;

                if (netinfo.dirty) {
                    netinfo.network = Networks.readNetwork(netinfo.backup);
                    initNetHeights(netinfo.network);
                    netinfo.dirty = false;
                }
                List<String> leavesIntersection = new ArrayList<>();
                for(Object leafObj : netToTry.getLeaves()) {
                    NetNode leaf = (NetNode) leafObj;
                    leavesIntersection.add(leaf.getName());
                }
                leavesIntersection.retainAll(netinfo.taxa);
                if(leavesIntersection.size() < 2) continue;

                double minScore = Double.MAX_VALUE;
                for(Network backbone : backbonesToTry) {
                    double curScore = compareTwoSubNetwork(netinfo.network, backbone, leavesIntersection);
                    minScore = Math.min(minScore, curScore);
                    if(minScore == 0) break;
                }
                score += minScore * netinfo.percentage;
            }

            score += netToTry.getReticulationCount() * 8.0;
            if (score < bestScore) {
                bestScore = score;
                bestNetworks.clear();
                bestNetworks.add(netToTry);
            } else if(score == bestScore) {
                bestNetworks.add(netToTry);
            }
        }
        return bestScore;
    }

    double ComputeBestScoreOpt(List<Network> networksToTry, List<Network> bestNetworks) {
        double bestScore = Double.MAX_VALUE;
        for(Network netToTry : networksToTry) {
            double score = 0.0;

            for (NetworkWithInfo netinfo : subnetworks_) {
                if(!netinfo.trustTime) continue;

                if (netinfo.dirty) {
                    netinfo.network = Networks.readNetwork(netinfo.backup);
                    initNetHeights(netinfo.network);
                    netinfo.dirty = false;
                }
                List<String> leavesIntersection = new ArrayList<>();
                for(Object leafObj : netToTry.getLeaves()) {
                    NetNode leaf = (NetNode) leafObj;
                    leavesIntersection.add(leaf.getName());
                }
                leavesIntersection.retainAll(netinfo.taxa);
                if(leavesIntersection.size() < 2) continue;

                Network subnet = getSubNetwork(netinfo.network, leavesIntersection, true).Item1;
                Network subnetNetToTry = getSubNetwork(netToTry.clone(), leavesIntersection, true).Item1;
                List<Network> backbonesToTry = new ArrayList<>();
                if(subnetNetToTry.getReticulationCount() > subnet.getReticulationCount())
                    backbonesToTry = Pipeline.getBackboneNetsWithNumReti(subnetNetToTry, subnet.getReticulationCount());
                else
                    backbonesToTry.add(subnetNetToTry);
                double minScore = Double.MAX_VALUE;
                for(Network backbone : backbonesToTry) {
                    double curScore = compareTwoSubNetwork(netinfo.network, backbone, leavesIntersection);
                    minScore = Math.min(minScore, curScore);
                    if(minScore == 0) break;
                }
                score += minScore * netinfo.percentage;
            }

            score += netToTry.getReticulationCount() * netToTry.getReticulationCount();
            if (score < bestScore) {
                bestScore = score;
                bestNetworks.clear();
                bestNetworks.add(netToTry);
            } else if(score == bestScore) {
                bestNetworks.add(netToTry);
            }
        }
        return bestScore;
    }

    Network compute() {
        Prepare();
        // Ad-hoc test
        /*{
            double bestScore = Double.MAX_VALUE;
            List<Network> bestNetworks = new ArrayList<>();
            List<Network> networksToTry = new ArrayList<>();
            networksToTry.add(Networks.readNetwork("(Z,(((((A,((I,H),((O)#H1,G))),((P)#H2,E)),((N)#H3,(K,J))),((L)#H4,B)),((C,(#H2,M)),(#H3,(#H4,(D,(#H1,F)))))));"));
            bestScore = ComputeBestScore(networksToTry, bestNetworks);
            System.out.println(bestScore);
        }*/

        for (NetworkWithInfo netinfo : subnetworks_) {
            //if(!netinfo.taxa.contains("Lygisaurus_macfarlani")) continue;
            //if(!netinfo.taxa.contains("Carlia_amax")) continue;
            //if(!netinfo.taxa.contains("Carlia_vivax")) continue;

            //if(!netinfo.taxa.contains("F")) continue;
            //if(!netinfo.taxa.contains("G")) continue;

            //System.out.println(Networks.getFullString(netinfo.network));

            //if(netinfo.taxa.contains("F") && netinfo.network.getReticulationCount() > 0)
            //    System.out.println(Networks.getFullString(netinfo.network));

            //if(netinfo.taxa.contains("Z") && netinfo.network.getReticulationCount() > 0)
            //    System.out.println(Networks.getFullString(netinfo.network));

            //if(netinfo.taxa.contains("I") && netinfo.network.getReticulationCount() > 0)
            //    System.out.println(Networks.getFullString(netinfo.network));
        }

        if(printDetails_)
            System.out.println("Start building...");
        Collections.sort(subnetworks_, (NetworkWithInfo a, NetworkWithInfo b)->Double.compare(b.PairwiseDistanceSum, a.PairwiseDistanceSum));
        //Collections.sort(subnetworks_, (NetworkWithInfo a, NetworkWithInfo b)->Double.compare(b.network.getRoot().getData().getHeight(), a.network.getRoot().getData().getHeight()));
        Network backbone = null;
        Set<String> curleaves = new HashSet<>();
        int skip = 0;

//        leavesUnderReticulation.remove(outgroup);
//        for(NetworkWithInfo candidate : subnetworks_ ) {
//            boolean flag = true;
//            for(String taxon : candidate.taxa) {
//                if(leavesUnderReticulation.contains(taxon) /*&& netinfo.network.getReticulationCount() == 0*/) {
//                    flag = false;
//                    break;
//                }
//            }
//            if(flag) {
//                if(skip > 0) {
//                    skip--;
//                    continue;
//                }
//
//                backbone = candidate.network;
//                curleaves.addAll(candidate.taxa);
//                candidate.dirty = true;
//                break;
//            }
//        }
//
//        if(backbone == null) {
//            for(NetworkWithInfo netinfo : subnetworks_ ) {
//                if(netinfo.taxa.contains(outgroup)) {
//                    if(skip > 0) {
//                        skip--;
//                        continue;
//                    }
//
//                    backbone = netinfo.network;
//                    curleaves.addAll(netinfo.taxa);
//                    netinfo.dirty = true;
//                    break;
//                }
//            }
//        }

        backbone = ElectBackbone3(curleaves);

        if(printDetails_) {
            System.out.println("Backbone taxa:");
            for (String s : curleaves) {
                System.out.println(s);
            }
            System.out.println();
        }

        for(String s : curleaves) {
            buildOrder.remove(s);
        }

        if(!trustReticulationTime) {
            for (Object nodeObj : backbone.dfs()) {
                NetNode node = (NetNode) nodeObj;
                if (node.isNetworkNode()) node.setData(null);
            }
        }

        if(printDetails_)
            System.out.println("Backbone: " + backbone);

        Collections.sort(subnetworks_, (NetworkWithInfo a, NetworkWithInfo b)->Integer.compare(a.network.getReticulationCount(), b.network.getReticulationCount()));

        //for(NetworkWithInfo netinfo : subnetworks_) {
        //    initNetHeights(netinfo.network);
        //    System.out.print("");
        //}

        while(curleaves.size() < leafnames_.size()) {
            //NetworkWithInfo target = null;
            String targetLeafName = buildOrder.get(0);
            if(printDetails_) {
                System.out.println("Target leaf name: " + targetLeafName);
            }
            buildOrder.remove(0);

            // Find the leaf want to add
//            for(NetworkWithInfo netinfo : subnetworks_) {
//                if(netinfo.dirty) continue;
//                Set<String> intersection = new HashSet<>();
//                intersection.addAll(curleaves);
//                intersection.retainAll(netinfo.taxa);
//                if(intersection.size() == netinfo.taxa.size() - 1) {
//                    if(netinfo.taxa.contains(targetLeafName) && !intersection.contains(targetLeafName)) {
//                        if(netinfo.network.getRoot().getData() == null) initNetHeights(netinfo.network);
//                        target = netinfo;
//                        break;
//                    }
//                }
//            }

            Map<Integer, NetworkWithInfo> subnetsToTry = new TreeMap<>();

            for(NetworkWithInfo netinfo : subnetworks_) {
                if(netinfo.dirty) continue;
                //if(!netinfo.trustTime) continue;
                if(netinfo.taxa.contains(targetLeafName)) {
                    if(netinfo.network.getRoot().getData() == null) initNetHeights(netinfo.network);

                    Set<NetNode> component1 = getComponent(netinfo, targetLeafName);
                    Set<NetNode> border1 = getBorder(component1);
//                    Set<NetNode> component2 = getComponent(target, targetLeafName);
//                    Set<NetNode> border2 = getBorder(component2);

                    if(subnetsToTry.containsKey(border1.size())) {
                        Set<NetNode> component3 = getComponent(subnetsToTry.get(border1.size()), targetLeafName);
                        if(component1.size() < component3.size()
                            || (component1.size() == component3.size() && netinfo.network.findNode(targetLeafName).getParents().iterator().next().getData().getHeight() < subnetsToTry.get(border1.size()).network.findNode(targetLeafName).getParents().iterator().next().getData().getHeight()))
                            subnetsToTry.put(border1.size(), netinfo);
                    } else {
                        subnetsToTry.put(border1.size(), netinfo);
                    }


                    //if(border1.size() > border2.size()
                    //        || (border1.size() == border2.size() && component1.size() < component2.size())
                    //|| ((border1.size() == border2.size() && component1.size() == component2.size()) && netinfo.network.findNode(targetLeafName).getParents().iterator().next().getData().getHeight() < target.network.findNode(targetLeafName).getParents().iterator().next().getData().getHeight())) target = netinfo;

                }
            }


            // Enumerate how to add
//            if(printDetails_) {
//                System.out.println(targetLeafName);
//                System.out.println("Target: " + target.network);
//                System.out.println(target.filename);
//            }
//            target.dirty = true;

            List<Network> networksToTry = new ArrayList<>();
            for(NetworkWithInfo netinfo : subnetsToTry.values()) {
                netinfo.dirty = true;
                Set<NetNode> component = getComponent(netinfo, targetLeafName);
                if(component.size() > 5) continue;
                Set<NetNode> border = cutBorder(component);
                networksToTry.addAll(enumerate_EHM(backbone, component, border, curleaves, targetLeafName));
                //networksToTry.addAll(enumerate(backbone, component, border, targetLeafName));

                if(printDetails_) {
                    System.out.println("Component size: " + component.size());
                    System.out.println("Border size: " + border.size());
                }
            }
            if(printDetails_) {
                System.out.println("Networks to try: " + networksToTry.size());
                System.out.println("SubNetworks to check: " + subnetworks_.size());
            }
            double bestScore = Double.MAX_VALUE;
            List<Network> bestNetworks = new ArrayList<>();

            bestScore = ComputeBestScoreOpt(networksToTry, bestNetworks);

            // Update backbone
            curleaves.add(targetLeafName);
            if(bestNetworks.size() == 0) return null;
            backbone = duplicateNetworkTopDown(bestNetworks.get(0));

            //initNetHeights(backbone);
            //resetHeightEHM(backbone);
            //resetBranchLengths(backbone);

            // Reassign heights
            Map<NetNode, Mean> heights = new HashMap<>();
            Map<Tuple<NetNode, NetNode>, Mean> probs = new HashMap<>();
            Map<NetNode, Set<NetNode>> full2inferred = new HashMap<>();
            List<NetworkWithInfo> toCheck = new ArrayList<>();
            toCheck.addAll(subnetworks_);
            toCheck.addAll(subnetworks_removed_);
            for (NetworkWithInfo netinfo : toCheck) {
                if (!netinfo.trustTime) continue;
                if (netinfo.dirty) {
                    netinfo.network = Networks.readNetworkWithRootPop(netinfo.backup);
                    netinfo.dirty = false;
                }

                List<String> leavesIntersection = new ArrayList<>();
                leavesIntersection.addAll(curleaves);
                leavesIntersection.retainAll(netinfo.taxa);
                if(leavesIntersection.size() < 3) continue;

                for(int i1 = 0 ; i1 < netinfo.taxa.size() ; i1++) {
                    if(!leavesIntersection.contains(netinfo.taxa.get(i1))) continue;

                    for(int i2 = i1 + 1 ; i2 < netinfo.taxa.size() ; i2++) {
                        if(!leavesIntersection.contains(netinfo.taxa.get(i2))) continue;

                        List<String> selected = new ArrayList<>();
                        selected.add(netinfo.taxa.get(i1));
                        selected.add(netinfo.taxa.get(i2));

                        Tuple<Network, Map<NetNode, NetNode>> tuple = getSubNetwork(backbone, selected, true);
                        Tuple<Network, Map<NetNode, NetNode>> tuple_netinfo = getSubNetwork(netinfo.network, selected, true);

                        if(Networks.hasTheSameTopology(tuple.Item1, tuple_netinfo.Item1)) {
                            Map<NetNode, NetNode> sub2full = tuple.Item2;
                            Map<NetNode, NetNode> sub2inferred = mapTwoNetworks(tuple.Item1, tuple_netinfo.Item1);

                            if(sub2inferred == null)
                                continue;

                            for(NetNode node : sub2inferred.keySet()) {
                                NetNode fnode = sub2full.get(node);
                                NetNode inode = sub2inferred.get(node);
                                if(!full2inferred.containsKey(fnode)) {
                                    full2inferred.put(fnode, new HashSet<>());
                                }
                                full2inferred.get(fnode).add(inode);

                                if(!heights.containsKey(fnode)) {
                                    heights.put(fnode, new Mean());
                                }
                                heights.get(fnode).increment(((NetNodeInfo)inode.getData()).getHeight());

                                if(node.isNetworkNode()) {
                                    for(Object parentObj : node.getParents()) {
                                        NetNode parent = (NetNode) parentObj;

                                        NetNode fparent = sub2full.get(parent);
                                        NetNode iparent = sub2inferred.get(parent);

                                        Tuple<NetNode, NetNode> ftuple = new Tuple<>(fnode, fparent);
                                        if(!probs.containsKey(ftuple)) {
                                            probs.put(ftuple, new Mean());
                                        }
                                        if(!inode.hasParent(iparent)) {
                                            continue;
                                        }

                                        probs.get(ftuple).increment(inode.getParentProbability(iparent));
                                    }
                                }
                            }
                        }
                    }
                }




            }

            if(reconcileHeights) {
                for (NetNode node : heights.keySet()) {
                    if (node.isLeaf()) continue;
                    double oldheight = node.getData() == null ? 0.0 : ((NetNodeInfo) node.getData()).getHeight();
                    node.setData(null);

                    double upperbound = getUpperBoundDfs(node) - 0.0001;
                    double lowerbound = getLowerBoundDfs(node) + 0.0001;
                    double newheight = heights.get(node).getResult();
                    if (newheight > upperbound) newheight = upperbound;
                    else if (newheight < lowerbound) newheight = lowerbound;

                    if (node.getData() == null) node.setData(new NetNodeInfo(newheight));
                    else ((NetNodeInfo) node.getData()).setHeight(newheight);
                    //for(NetNode inode : full2inferred.get(node)) {
                    //    ((NetNodeInfo)inode.getData()).setHeight(heights.get(node).getResult());
                    //}
                    System.out.print("");
                }
                if(printDetails_)
                    System.out.println("Reset branch length: " + resetBranchLengths(backbone));
                else
                    resetBranchLengths(backbone);
                for (NetworkWithInfo netinfo : subnetworks_) {
                    //if (!netinfo.trustTime) continue;
                    //resetBranchLengths(netinfo.network);
                }
                ComputeExtendedHeightMatrix();
            }
            for(Object nodeObj : backbone.dfs()) {
                NetNode node = (NetNode) nodeObj;
                for(Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;

                    Tuple<NetNode, NetNode> tuple = new Tuple<>(node, parent);
                    if(node.isNetworkNode()) {
                        if (probs.containsKey(tuple) && !Double.isNaN(probs.get(tuple).getResult()) && probs.get(tuple).getResult() >= 0 && probs.get(tuple).getResult() <= 1) {
                            node.setParentProbability(parent, probs.get(tuple).getResult());
                            node.setParentProbability(Networks.getOtherParent(node, parent), 1.0 - probs.get(tuple).getResult());
                            break;
                        } else {
                            node.setParentProbability(parent, NetNode.NO_PROBABILITY);
                            node.setParentProbability(Networks.getOtherParent(node, parent), NetNode.NO_PROBABILITY);
                        }
                    } else {
                        node.setParentProbability(parent, NetNode.NO_PROBABILITY);
                    }
                }
            }


            // Remove fulfilled inferred subnetworks
            // Remove subnetworks with problematic reticulations
            List<NetworkWithInfo> toRemove = new ArrayList<>();
            int numRetiAbove = NetworkUtils.GetNumAllRetiAboveLeaf(backbone, targetLeafName);
            for (NetworkWithInfo netinfo : subnetworks_) {
                List<String> leavesIntersection = new ArrayList<>();
                leavesIntersection.addAll(curleaves);
                leavesIntersection.retainAll(netinfo.taxa);
                if(leavesIntersection.size() == netinfo.taxa.size())
                    toRemove.add(netinfo);
                else if(netinfo.taxa.contains(targetLeafName)){
                    int curnum = NetworkUtils.GetNumAllRetiAboveLeaf(Networks.readNetworkWithRootPop(netinfo.backup), targetLeafName);
                    if(curnum > numRetiAbove) {
                        toRemove.add(netinfo);
                    }
                }
            }
            subnetworks_.removeAll(toRemove);
            subnetworks_removed_.addAll(toRemove);


            if(printDetails_) {
                System.out.println("Current number of leaves: " + curleaves.size());
                System.out.println("Best score: " + bestScore);
                System.out.println(Networks.getTopologyString(backbone));
                System.out.println(Networks.getDendroscopeCompatibleString(backbone));
                System.out.println(Networks.getFullString(backbone));
                System.out.println();
            }
        }
        backbone.getRoot().setRootPopSize(popsize);
        return backbone;
    }

    public static void main(String[] args) {
        final SuperNetwork superNetwork0 = new SuperNetwork(new ArrayList<>());
        Random random = new Random(12345);
        for(int i = 0 ; i < 20000 ; ) {
            SuperNetwork superNetwork = new SuperNetwork(new ArrayList<>());
            Network trueNetwork = superNetwork.genRandomNetwork(20, 3);//Networks.readNetwork("(((((((C:1.0,D:1.0)I1:1.0)I2#H1:2.0,B:4.0)I4:3.0,(I2#H1:3.0,E:5.0)I5:2.0)I7:1.0)I8#H2:1.0,A:9.0)I9:2.0,(((F:3.0,G:3.0)I3:3.0,H:6.0)I6:4.0,I8#H2:2.0)I10:1.0)I11;");
            //trueNetwork = Networks.readNetwork("((J:1.809522611909963,E:1.809522611909963):47.443779683538274,((((((B:8.828999653685607,(((Q:3.1360710730314283,(F:0.8037659932417116)#H3:2.332305079789717):3.0386332557943403,A:6.174704328825769):0.7924204274196303)#H1:1.8618748974402077):1.703767690022044,D:10.53276734370765):4.703185062927272,(I:13.893298700567877,M:13.893298700567877):1.3426537060670451):18.315790815909445,(((L:18.915906609467235,O:18.915906609467235):9.640850499772096,(G:26.58996650181392,(S:24.07610666873296,((R:20.61168828635986,C:20.61168828635986):1.6781330084325248,(H:13.946056937106434)#H2:8.343764357685952):1.7862853739405757):2.513859833080957):1.9667906074254127):3.712924589510809,N:32.26968169875014):1.2820615237942263):9.637188428655186,#H3:42.38516565795784):3.405606828619476,(((K:37.07815781017047,#H2:23.132100873064033):6.051880762067171,(T:40.43214528359408,P:40.43214528359408):2.697893288643556):2.5300562958281816,#H1:38.69297011182042):0.9344436117532098):2.6587638156292073);");
            //trueNetwork = Networks.readNetwork("((((((T:1.133394744106908,O:1.133394744106908):14.01105789850394,((L:11.886220106106931,(B:7.922000232953535,(((F:4.478626382697977,J:4.478626382697977):0.5763156066476913)#H2:0.5969268903304865)#H1:2.270131353277381):3.9642198731533957):1.6497798510124717,(S:5.711901499653102)#H3:7.824098457466301):1.6084526854914465):5.817170573438924,(H:17.90964277139056,#H3:12.197741271737456):3.0519804446592147):14.746014637337993,((A:24.28845834008711,I:24.28845834008711):9.863752616521403,(E:30.4898896359965,(M:26.840932210788857,(D:7.124877991145295,#H1:1.4730091114691408):19.71605421964356):3.6489574252076444):3.662321320612012):1.5554268967792524):16.946202450449235,(G:50.348416462918735,(N:46.90669342965074,(C:44.12661457891803,((Q:37.94044455538827,R:37.94044455538827):2.920644877746163,P:40.86108943313443):3.2655251457835988):2.780078850732707):3.4417230332679978):2.305423840918266):1.6642900901934965,(K:39.263036689399094,#H2:34.20809470005342):15.055093704631403);");
            //trueNetwork = Networks.readNetwork("(((B:3.0586500098498957,(((T:0.21737681332974473)#H2:0.6506165917066502)#H1:1.490025662977076,#H2:2.1406422546837263):0.7006309418364247):10.288933052647376,((A:11.237178945612342)#H3:0.86901757748252,(((N:5.97545996630095,F:5.97545996630095):3.480532452392988,I:9.455992418693938):1.0828222348081216,#H1:9.670821248465664):1.5673818695928023):1.2413865394024093):38.244585817584905,((((((Q:15.311501453242885,H:15.311501453242885):7.9553067409298865,((G:17.921294753553088,L:17.921294753553088):1.813935856697416,O:19.735230610250504):3.5315775839222674):6.476454664895083,(C:25.820809702152985,#H3:14.583630756540643):3.9224531569148695):6.735204119920219,((M:33.40394671544683,P:33.40394671544683):1.5301595614594774,E:34.934106276906306):1.544360702081768):4.472187734559647,(K:37.5650879480432,S:37.5650879480432):3.385566765504521):6.6558926194892365,((R:42.45451063917688,D:42.45451063917688):3.8923200366077424,J:46.346830675784624):1.2597166572523335):3.98562154704522);");
            //trueNetwork = Networks.readNetwork("((((A:2.6049934387824605,T:2.6049934387824605):7.992421957251229,(Q:7.801423397786972,(H:5.3797807017101364,I:5.3797807017101364):2.4216426960768356):2.795991998246717):2.116131427512247,B:12.713546823545936):33.01605990339996,((M:14.28905695216133,E:14.28905695216133):27.644943258597912,(((D:17.86665283292104,S:17.86665283292104):1.6702946814656094,(G:9.644962078989135)#H1:9.891985435397515):21.230926011784906,((O:23.829695590850886,((C:21.367153727311322,P:21.367153727311322):1.307934205997249,(N:4.293091461266711)#H2:18.381996472041862):1.1546076575423143):14.45393640344557,((L:27.01958938067096,(#H1:1.7345921399294593,#H2:7.0864627576518835):15.640035161752367):8.194271252791307,((((K:28.67070240578753,(R:10.797087835984318)#H3:17.873614569803213):2.4216393087428294,J:31.09234171453036):1.7709764862459174,F:32.86331820077628):1.3049271171019043,#H3:23.371157481893864):1.0456153155840866):3.0697713608341886):2.4842415318750994):1.166126684587688):3.795606516186652);");
            //trueNetwork = Networks.readNetwork("(((((((K:22.617658768976227,(N:19.05301040791268,D:19.05301040791268):3.564648361063547):2.6730114239527296,(P:9.838367146190123)#H1:15.452303046738834):1.165122333694061,M:26.455792526623018):8.03463275092047,((I:30.140085385664086,(R:11.40529595405491)#H2:18.734789431609176):2.0292316025167096,Q:32.169316988180796):2.321108289362691):5.031153625903066,(J:38.01234738277146,S:38.01234738277146):1.509231520675094):2.995252755803527,((A:41.33618751525881,(#H2:10.056647364919943,(L:13.747346663796062,#H1:3.908979517605939):7.7145966551787915):19.874244196283954):0.4458050684497721,(((E:3.7348407863349498,H:3.7348407863349498):11.940435110049194,((C:4.772570008201653,G:4.772570008201653):9.261526551109947,((T:7.09265639242227,O:7.09265639242227):4.362560651863441,(B:9.177172293411338,F:9.177172293411338):2.278044750874372):2.578879515025889):1.6411793370725434):12.622221907053337)#H3:13.484494780271099):0.7348390755415011):2.242948671491824,#H3:16.462282527304424);");
            //trueNetwork = Networks.readNetwork("((A:6.8191196535208665,(O:5.287988577107271,(S:3.1176958307989935,(N:1.4656508129431318)#H1:1.6520450178558617):2.1702927463082773):1.5311310764135957):42.54198335861826,((((F:40.678920285213806,I:40.678920285213806):1.7951465699313829,Q:42.47406685514519):3.46883920543209,((((J:9.151141986739185,C:9.151141986739185):9.531336338509673,(G:17.549437924672333,(E:13.885915210910277,(((B:10.403285412623116,K:10.403285412623116):1.9302688441885731,M:12.33355425681169):0.6970892583766872)#H3:0.8552716957218998):3.663522713762056):1.133040400576526):18.77570238982458,(H:35.63942148432256,((((D:21.87803272712815,R:21.87803272712815):7.003856507976032,(T:24.984779224678846,P:24.984779224678846):3.8971100104253367):3.8843810326606913,L:32.766270267764874):2.0928272616377157,#H3:21.828454014214213):0.7803239549199716):1.818759230750878):2.2947591321771483)#H2:6.1899662133266915):1.9031873008730642,(#H2:4.273529317766474,#H1:42.56081835207393):3.8196241964332813):1.515009650688782);");
            //trueNetwork = Networks.readNetwork("((((Q:1.8510447476068266,H:1.8510447476068266):5.980526316451876,(O:4.831431720779742,I:4.831431720779742):3.0001393432789607):10.05515862582396,((C:4.6440918069414305)#H1:10.385179963495576,((((A:11.559272884379528,N:11.559272884379528):0.8211261579441924,#H1:7.73630723538229):0.5639181928523431,(G:5.784512729936338)#H3:7.159804505239725):0.5671717890694961)#H2:1.5177827461914468):2.8574579194456575):30.94236725362942,((((S:34.22880828415031,((K:19.471057409620833,J:19.471057409620833):12.474503094814168,(((M:23.06153836548164,L:23.06153836548164):4.2883607751067565,((P:24.160675607479373,R:24.160675607479373):2.0244898708075,#H2:12.673676454041313):1.164733662301522):2.04065897241847,D:29.390558113006865):2.555002391428136):2.2832477797153103):6.2898647033427295,(B:37.5355658042458,T:37.5355658042458):2.983107183247242):4.342096477233369,(E:43.455701248281315,F:43.455701248281315):1.4050682164450947):1.4827480572436826,#H3:40.559004792033754):2.485579421541992);");
            trueNetwork = Networks.readNetwork("(Z:1000.0,((((((G:2.9859839999999997,((I:2.0736,J:2.0736)S20:0.41472,H:2.48832)S19:0.4976639999999999)S11:15.502441889503631,(((((((P:1.0)#H5:0.728::0.36,K:1.728)S17:1.8551807999999996,(M:1.44,L:1.44)S18:2.1431807999999997)S16:1.5765995519999998,D:5.159780351999999)S15:1.0319560703999997,C:6.191736422399999)S12:4.507584115507197)#H4:4.707701036679165::0.63)#H3:3.08140431491727::0.62)S10:8.1349073913816)#H1:11.7142666435895::0.4,(((F:4.299816959999999,E:4.299816959999999)S14:3.1302667468799994,B:7.430083706879999)S13:1.486016741375999)#H2:29.421499476218735::0.39)S6:16.868543966768875,#H3:39.79912231665725::0.38)S3:24.29070331214718,((O:1.2,N:1.2)S5:65.04737266949232,((A:12.839184645488634,#H4:2.139864107581438::0.37)S7:33.16593526388104,(#H1:5.324666656177044::0.6,(#H5:21.186111067404358::0.64,#H2:13.27001061914836::0.61)S9:9.761888869657916)S8:14.057119972307401)S4:20.242252760122646)S2:13.249474533898464)S1:920.5031527966092);");
            System.out.println("True network: " + Networks.getDendroscopeCompatibleString(trueNetwork));

            NetworkUtils.genAllSubNetworks(trueNetwork, 3);
            List<Tuple3<Network,String,Double>> newnetlist = new ArrayList<>();

            for (Network net : superNetwork._subnetworks) {
                newnetlist.add(new Tuple3<>(net, "", 1.0));
                if(net.findNode("I") == null) continue;
                if(net.findNode("T") == null) continue;
                if(net.findNode("B") == null) continue;
                System.out.println(net.toString());
            }


            SuperNetwork3 superNetwork3 = new SuperNetwork3(newnetlist);
            Network result = superNetwork3.compute();

            boolean correctness = Networks.hasTheSameTopology(result, trueNetwork);
            System.out.println("Correctness: " + correctness);
            if(!correctness) {
                break;
            }
        }
    }
}
