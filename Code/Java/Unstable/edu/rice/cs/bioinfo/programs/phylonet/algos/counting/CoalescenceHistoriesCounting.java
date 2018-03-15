package edu.rice.cs.bioinfo.programs.phylonet.algos.counting;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SpeciesNetPriorDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.dimension.AddReticulation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.DataGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.introgression.Coalescence;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 5/24/17
 * Time: 12:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class CoalescenceHistoriesCounting {

    public static int getCoalescenceHistoriesCount(Network net, Tree gt) {
        Network<Double> trueNetwork = Networks.readNetwork(net.toString());
        List<Tree> gts = new ArrayList<>();
        gts.add(gt);
        Coalescence coal = new Coalescence(trueNetwork, gts, null);
        List<Integer> count = coal.getAllCoalescenceHistoriesCount();
        return count.get(0);
    }

    public static int getNumTaxaUnderReticulation(Network net) {
        int count = 0;
        for(Object leafObject : net.getLeaves()) {
            NetNode leaf = (NetNode) leafObject;
            NetNode node = leaf;
            while(!node.isRoot()) {
                if(node.isNetworkNode()) {
                    count++;
                    break;
                }
                node = (NetNode) node.getParents().iterator().next();
            }
        }
        return count;
    }

    public static int getDiscreteDiameterOneNode(NetNode p, Network net) {
        Map<NetNode, Integer> marks = new HashMap<>();
        NetNode p1 = null;
        NetNode p2 = null;
        if(!p.isNetworkNode())
            return 0;
        for(Object parentO : p.getParents()) {
            NetNode parent = (NetNode) parentO;
            if(p1 != null)
                p2 = parent;
            else
                p1 = parent;
        }

        Queue<NetNode> q = new LinkedList<>();
        q.add(p1);
        marks.put(p1, 0);
        while(!q.isEmpty()) {
            NetNode node = q.poll();
            int d = marks.get(node);

            for(Object parentO : node.getParents()) {
                NetNode parent = (NetNode) parentO;
                q.add(parent);
                if(!marks.containsKey(parent))
                    marks.put(parent, Integer.MAX_VALUE);

                marks.put(parent, Math.min(marks.get(parent), d + 1));
            }
        }


        Map<NetNode, Integer> marks2 = new HashMap<>();
        marks2.put(p2, 0);
        int s = Integer.MAX_VALUE;
        NetNode best = null;
        q.clear();
        q.add(p2);
        while(!q.isEmpty()) {
            NetNode node = q.poll();
            if(marks.containsKey(node) && s > marks.get(node)) {
                s = marks.get(node);
                best = node;
            }
            int d = marks2.get(node);

            for(Object parentO : node.getParents()) {
                NetNode parent = (NetNode) parentO;
                q.add(parent);

                if(!marks2.containsKey(parent))
                    marks2.put(parent, Integer.MAX_VALUE);

                marks2.put(parent, Math.min(marks2.get(parent), d + 1));
            }
        }

        s = marks.get(best) + marks2.get(best) + 2;

        return s;


    }

    public static int getDiscreteDiameter(Network net) {
        int result = 0;
        for(Object nodeObject : net.dfs()) {
            NetNode node = (NetNode) nodeObject;
            NetNode networkNode = null;
            if(node.isNetworkNode()) {
                networkNode = node;
            } else
                continue;

            result += getDiscreteDiameterOneNode(networkNode, net);

        }



        return result;
    }

    private static boolean isReticulationDependent(Network net) {
        boolean flag = false;
        for(Object nodeObject : net.dfs()) {
            NetNode node = (NetNode) nodeObject;
            NetNode networkNode = null;
            if(node.isNetworkNode()) {
                networkNode = node;
            } else
                continue;

            Queue<NetNode> q = new LinkedList<>();
            q.add(networkNode);
            while(!q.isEmpty()) {
                NetNode cnode = q.poll();

                for(Object parentO : cnode.getParents()) {
                    NetNode parent = (NetNode) parentO;
                    q.add(parent);
                    if(parent.isNetworkNode())
                        flag = true;
                }
            }

        }
        return flag;
    }

    private static void countPathsHelper(ArrayList<Tuple<NetNode, NetNode>> allEdges, Map<Tuple<NetNode, NetNode>, Integer> edgeIndices, int startEdge, int currentEdge, int f[]) {
        if(currentEdge == startEdge) {
            f[currentEdge] = 1;
            return;
        }

        f[currentEdge] = 0;
        NetNode fromNode = allEdges.get(currentEdge).Item2;
        for(Object toNodeO : fromNode.getChildren()) {
            NetNode toNode = (NetNode) toNodeO;
            int edgeIndex = edgeIndices.get(new Tuple<>(fromNode, toNode));
            countPathsHelper(allEdges, edgeIndices, startEdge, edgeIndex, f);
            f[currentEdge] += f[edgeIndex];
        }
    }

    public static long getCoalescenceHistoriesCount2(Network net, Tree gt) {
        ArrayList<Tuple<NetNode, NetNode>> allEdges = new ArrayList<>();
        Map<Tuple<NetNode, NetNode>, Integer> edgeIndices = new HashMap<>();
        Map<String, Integer> externalEdges = new HashMap<>();
        int rootEdge = -1;
        for(Object nodeO: Networks.postTraversal(net)){
            NetNode node = (NetNode) nodeO;

            if(node.isRoot()) {
                allEdges.add(new Tuple<>(node, node));
                edgeIndices.put(new Tuple<>(node, node), allEdges.size() - 1);
                rootEdge =  allEdges.size() - 1;
                continue;
            }
            for(Object parentO : node.getParents()) {
                NetNode parent = (NetNode) parentO;
                allEdges.add(new Tuple<>(parent, node));
                edgeIndices.put(new Tuple<>(parent, node), allEdges.size() - 1);
                if(node.isLeaf()) {
                    externalEdges.put(node.getName(), allEdges.size() - 1);
                }
            }
        }
        int nEdges = allEdges.size();

        long nPaths[][] = new long[nEdges][nEdges]; // number of paths from edge to edge
        for(int i = 0 ; i < nEdges ; i++) {
            int f[] = new int[nEdges];
            countPathsHelper(allEdges, edgeIndices, i, rootEdge, f);
            for(int j = 0 ; j < nEdges ; j++) {
                if(f[j] > 0)
                    nPaths[i][j] = f[j];
            }
        }

        List<STITreeCluster> clusters = gt.getClusters(null, true);
        clusters.add(new STITreeCluster(clusters.get(0)));
        // add the cluster that contains all taxa
        clusters.get(clusters.size() - 1).getCluster().set(0, clusters.get(clusters.size() - 1).getTaxa().length);
        // store the two subtrees of each cluster
        List<List<STITreeCluster>> clusterContains = new ArrayList<>();
        for(int i = 0 ; i < clusters.size() ; i++) {
            clusterContains.add(new ArrayList<>());
            for(int j = 0 ; j < clusters.size() ; j++) {
                for(int k = 0 ; k < clusters.size() ; k++) {
                    if (i == j || k == i || k == j) continue;

                    BitSet set = (BitSet) clusters.get(j).getCluster().clone();
                    set.or(clusters.get(k).getCluster());

                    // if [i] == [j] & [k] -- find two subtrees of [i]
                    if (set.equals(clusters.get(i).getCluster())) {
                        clusterContains.get(i).add(clusters.get(j));
                        clusterContains.get(i).add(clusters.get(k));
                        break;
                    }
                }
                if(clusterContains.get(i).size() > 0)
                    break;
            }

            // if no subtree found then one subtree of [i] must be a leaf
            if(clusterContains.get(i).size() == 0) {
                for(int j = 0 ; j < clusters.size() ; j++) {
                    if(i == j) continue;

                    // if card([i]) != card([j]) + 1
                    if(clusters.get(j).getCluster().cardinality() + 1 != clusters.get(i).getCluster().cardinality())
                        continue;
                    BitSet set = (BitSet) clusters.get(j).getCluster().clone();
                    set.and(clusters.get(i).getCluster());

                    // [i] is one taxon more than [j]
                    if (set.equals(clusters.get(j).getCluster())) {
                        clusterContains.get(i).add(clusters.get(j));
                        break;
                    }
                }
            }
        }

        // rho(cluster,edge)
        Map<STITreeCluster, long[]> nCoalWays = new HashMap<>();
        for(int i = 0 ; i < clusters.size() ; i++) {
            nCoalWays.put(clusters.get(i), new long[nEdges]);
            if(clusterContains.get(i).size() == 0) { // two subtrees are both leaf ([i] is a lowest cluster)
                int t = clusters.get(i).getCluster().nextSetBit(0);
                int e1 = externalEdges.get(clusters.get(i).getTaxa()[t]);
                t = clusters.get(i).getCluster().nextSetBit(t + 1);
                int e2 = externalEdges.get(clusters.get(i).getTaxa()[t]);

                for(int j = 0 ; j < nEdges ; j++) { // "outer" edges
                    if(j == e2 || j == e1) continue;
                    nCoalWays.get(clusters.get(i))[j] = nPaths[e1][j] * nPaths[e2][j];
                }
            } else if(clusterContains.get(i).size() == 1) { // one subtree is leaf
                BitSet subset = (BitSet) clusterContains.get(i).get(0).getCluster().clone();
                subset.xor(clusters.get(i).getCluster());
                int t = subset.nextSetBit(0);
                int e1 = externalEdges.get(clusters.get(i).getTaxa()[t]);

                for(int j = 0 ; j < nEdges ; j++) { // "outer" edges
                    for(int e2 = 0 ; e2 < nEdges ; e2++) { // "inner" edges
                        //if(j == e1 || j == e2) continue;
                        nCoalWays.get(clusters.get(i))[j] += nPaths[e1][j]
                                * nCoalWays.get(clusterContains.get(i).get(0))[e2] * nPaths[e2][j];
                    }
                }
            } else if(clusterContains.get(i).size() == 2) {
                for(int j = 0 ; j < nEdges ; j++) { // "outer" edges
                    for(int e1 = 0 ; e1 < nEdges ; e1++) { // "inner" edges
                        for(int e2 = 0 ; e2 < nEdges ; e2++) { // "inner" edges
                            //if(j == e1 || j == e2) continue;
                            nCoalWays.get(clusters.get(i))[j] +=
                                    nCoalWays.get(clusterContains.get(i).get(0))[e1] * nPaths[e1][j]
                                    * nCoalWays.get(clusterContains.get(i).get(1))[e2] * nPaths[e2][j];
                        }
                    }
                }
            }
        }

        long sum = 0;
        for(int i = 0 ; i < nEdges ; i++) {
            sum += nCoalWays.get(clusters.get(clusters.size() - 1))[i] * nPaths[i][rootEdge];
        }

        return sum;
    }

    public static void test1(String[] args) {
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

        double gamma = 1 - 0.05;

        double y = 0.1;
        double x = 10000;

        Network<Double> trueNetwork;
        trueNetwork = Networks.readNetwork("((((b:1000.0,c:1000.0)I4:" + y + ")I3#H1:0.0::" + (1 - gamma) + ",a:" + (1000 + y) + ")I1:" + x + ",(I3#H1:0.0::" + gamma + ",d:" + (1000 + y) + ")I2:"+ x +")I0;");


        Coalescence coal = new Coalescence(trueNetwork, gts, null);
        List<Integer> count = coal.getAllCoalescenceHistoriesCount();
        for(int i = 0 ; i < count.size() ; i++) {
            System.out.println(gts.get(i).toString() + " : " + count.get(i) + " " + getCoalescenceHistoriesCount2(trueNetwork, gts.get(i)));
        }

    }

    public static void test2() {
        GeneTreeProbability geneTreeProbability = new GeneTreeProbability();
        List<Double> gtprob = new ArrayList<>();
        List<Tree> gts = new ArrayList<>();
        gts.add(Trees.readTree("(a, (b, c));"));
        gts.add(Trees.readTree("(b, (a, c));"));
        gts.add(Trees.readTree("(c, (a, b));"));

        double gamma = 1 - 0.05;

        double y = 0.1;
        double x = 10000;

        Network<Double> trueNetwork;
        trueNetwork = Networks.readNetwork("(((b:" + (1000.0 + y) + ")I3#H1:0.0::" + (1 - gamma) + ",a:" + (1000 + y) + ")I1:" + x + ",(I3#H1:0.0::" + gamma + ",c:" + (1000 + y) + ")I2:"+ x +")I0;");
        System.out.println("Network = " + trueNetwork.toString());

        Coalescence coal = new Coalescence(trueNetwork, gts, null);
        List<Integer> count = coal.getAllCoalescenceHistoriesCount();
        for(int i = 0 ; i < count.size() ; i++) {
            System.out.println(gts.get(i).toString() + " : " + count.get(i) + " " + getCoalescenceHistoriesCount2(trueNetwork, gts.get(i)));
        }
        System.exit(0);
    }

    public static void processResult() {

        double data[][] = new double[21][40];
        for(int i = 0 ; i < data.length ; i++) {
            for(int j = 0 ; j < data[0].length; j++) {
                data[i][j] = -10;
            }
        }

        String path = "coal_history_log_multigt.txt";

        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            int i = 0;
            while ((line = br.readLine()) != null) {
                if(line.equals("GTs")) {
                    while(!line.equals(""))
                        line = br.readLine();
                    continue;
                }

                String items[] = line.split("\t");
                int a = Integer.parseInt(items[0]);
                int b = Integer.parseInt(items[1]);
                long c = Long.parseLong(items[2]);
                double d = Double.parseDouble(items[3]);
                Network e = Networks.readNetwork(items[4]);
                //if(isReticulationDependent(e))
                    data[a][b] = Math.max(data[a][b], d);
            }
        } catch(Exception e) {
            e.printStackTrace();
        }

        System.out.print('{');
        for(int i = 0 ; i < 21 ; i++) {
            if(i > 0)
                System.out.print(',');
            System.out.print('{');

            for(int j = 0 ; j < 21 ; j++) {

                if(j > 0)
                    System.out.print(',');
                System.out.print(data[i][j]);
            }
            System.out.print('}');
            //System.out.println();
        }
        System.out.print('}');
    }

    public static void robustnessCheck() {
        List<String> treestrings = new ArrayList<>();
        int n = 30;

        String[] leaves;
        int taxa = 7;
        leaves = new String[taxa];
        for(int i = 0 ; i < taxa ; i++) {
            leaves[i] = Integer.toString(i + 1);
        }

        treestrings.clear();
        for(int i = 0 ; i < n ; i++) {
            Tree randomTree = Trees.generateRandomTree(Arrays.copyOfRange(leaves, 0, taxa));
            treestrings.add(randomTree.toNewick());
        }

        int numGTs = 100;
        List<String> gtPool = new ArrayList<>();
        for(int i = 0 ; i < numGTs ; i++) {
            Tree randomTree = Trees.generateRandomTree(Arrays.copyOfRange(leaves, 0, taxa));
            gtPool.add(randomTree.toNewick());
        }

        for(String treestring : treestrings) {
            System.out.println("Tree: " + treestring);

            DataGenerator datagen = new DataGenerator();

            Network net1 = Networks.readNetwork(treestring);
            List<Network> netList = datagen.addReticulations(net1, 10);

            List<Network> netList2 = new ArrayList<>();
            for(int i = 0 ; i < netList.size() ; i++) {
                netList2.addAll(datagen.addReticulations(netList.get(i), 10));
            }
            System.out.println(netList2.size());
            for(Network network : netList2) {
                //network = Networks.readNetwork("((((8,(4,(3,2)I5)I3)I1)I10#H1)I9,I10#H1)I0;");
                //tree1 = Trees.readTree("(((2,3),4),8);");
                for(String gtstring: gtPool) {
                    Tree tree1 = Trees.readTree(gtstring);
                    //System.out.println(network + "\t" + tree1);
                    long k1 = getCoalescenceHistoriesCount(network, tree1);
                    long k2 = getCoalescenceHistoriesCount2(network, tree1);
                    //System.out.println(k1);
                    if (k1 != k2) {
                        System.out.println(k1 + " " + k2 + "!!!!!");
                        System.out.println(network);
                        System.out.println(tree1);
                    }
                }
            }

        }
    }

    public static void main(String[] args) {
        test2();
        if(true)
        {
            //robustnessCheck();
            processResult();
            return;
        }
        //test1(args);

        String path = "/Users/zhujiafan/Documents/SharedFolder/20taxa.trees";
        List<String> treestrings = new ArrayList<>();
        int n = 30;
        /*try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            int i = 0;
            while ((line = br.readLine()) != null) {
                i++;
                if(i > n) break;
                treestrings.add(line);
            }
        } catch(Exception e) {
            e.printStackTrace();
        }*/

        String[] leaves;
        int taxa = 20;
        leaves = new String[taxa];
        for(int i = 0 ; i < taxa ; i++) {
            leaves[i] = Integer.toString(i + 1);
        }

        treestrings.clear();
        for(int i = 0 ; i < n ; i++) {
            Tree randomTree = Trees.generateRandomTree(Arrays.copyOfRange(leaves, 0, taxa));
            treestrings.add(randomTree.toNewick());
        }

        int numGTs = 100;
        List<String> gtPool = new ArrayList<>();
        for(int i = 0 ; i < numGTs ; i++) {
            Tree randomTree = Trees.generateRandomTree(Arrays.copyOfRange(leaves, 0, taxa));
            gtPool.add(randomTree.toNewick());
        }

        try {
            PrintWriter out = new PrintWriter("coal_history_log_multigt.txt");
            //PrintStream out = System.out;

            out.println("GTs");
            for(String gtstring : gtPool) {
                out.println(gtstring);
            }
            out.println();

            for(String treestring : treestrings) {
                System.out.println("Tree: " + treestring);
                Tree tree1 = Trees.readTree(treestring);

                //UltrametricNetwork net = new UltrametricNetwork(treestring);
                //System.out.println(getCoalescenceHistoriesCount(net.getNetwork(), tree1));
                //SpeciesNetPriorDistribution prior = new SpeciesNetPriorDistribution(1);
                //AddReticulation op = new AddReticulation(net);
                DataGenerator datagen = new DataGenerator();
                //prior.logPrior(net.getNetwork());

                Network net1 = Networks.readNetwork(treestring);
                //List<Network> netList = datagen.allPossibleReticulation(net1);
                List<Network> netList = datagen.addReticulations(net1, 10);
                netList.add(0, net1.clone());

                List<Network> netList2 = new ArrayList<>();
                for(int i = 0 ; i < netList.size() ; i++) {
                    //netList2.addAll(datagen.addReticulations(netList.get(i), 10));
                }
                netList2.add(0, net1.clone());
                System.out.println(netList.size());
                long treeCoalHistories = 0;
                long treesum = 0;

                for(Network network : netList) {
                    //network = Networks.readNetwork("((((8,(4,(3,2)I5)I3)I1)I10#H1)I9,I10#H1)I0;");
                    //tree1 = Trees.readTree("(((2,3),4),8);");
                /*long k1 = getCoalescenceHistoriesCount(network, tree1);
                long k2 = getCoalescenceHistoriesCount2(network, tree1);
                if(k1 != k2) {
                    System.out.println(k1 +  " " + k2);
                    System.out.println(network);
                    System.out.println(tree1);
                }*/
                    long sum = 0;
                    for(String gtstring: gtPool) {
                        tree1 = Trees.readTree(gtstring);
                        long coalHistories = getCoalescenceHistoriesCount2(network, tree1);
                        sum += coalHistories;

                    }
                    if (getNumTaxaUnderReticulation(network) == 0) {
                        treesum = sum;
                    }
                    out.println(getNumTaxaUnderReticulation(network) + "\t" + getDiscreteDiameter(network) + "\t" + sum + "\t" + sum * 1.0 / treesum + "\t" + network.toString() );

                }


            }


            out.close();
        } catch (Exception ioe) {
            ioe.printStackTrace();
        }


    }

}
