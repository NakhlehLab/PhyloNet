package edu.rice.cs.bioinfo.programs.phylonet.algos.counting;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 4/15/18
 * Time: 4:55 PM
 * To change this template use File | Settings | File Templates.
 */

// It counts REDUCED networks. When you are running this, change BniNetNode to allow parallel edges!!!
public class NetworkTopologiesCounting {
    public static void enumerateTreeTopologiesHelper(int nTaxa, List<Tree> trees, List<STITreeCluster> remainingClusters, List<STITreeCluster> currentClusters, String taxa[]) {
        if(remainingClusters.size() == 1) {
            Tree newtree = Trees.buildTreeFromClusters(currentClusters);
            for(int i = 0 ; i < trees.size() ; i++) {
                if(Trees.haveSameRootedTopology(newtree, trees.get(i))) {
                    return;
                }
            }
            trees.add(newtree);
            return;
        }

        for(int i = 0 ; i < remainingClusters.size() ; i++) {
            for(int j = i + 1 ; j < remainingClusters.size() ; j++) {
                STITreeCluster tc = new STITreeCluster(taxa);
                tc = tc.merge(remainingClusters.get(i));
                tc = tc.merge(remainingClusters.get(j));
                List<STITreeCluster> nextRemainingClusters = new ArrayList<>();
                nextRemainingClusters.addAll(remainingClusters);
                nextRemainingClusters.remove(remainingClusters.get(i));
                nextRemainingClusters.remove(remainingClusters.get(j));
                nextRemainingClusters.add(tc);

                List<STITreeCluster> nextCurrentClusters = new ArrayList<>();
                nextCurrentClusters.addAll(currentClusters);
                nextCurrentClusters.add(tc);

                enumerateTreeTopologiesHelper(nTaxa, trees, nextRemainingClusters, nextCurrentClusters, taxa);
            }
        }
    }

    public static List<Tree> enumerateTreeTopologies(int nTaxa) {
        List<STITreeCluster> clusters = new ArrayList<>();
        String[] taxa = new String[nTaxa];
        for(int i = 0 ; i < nTaxa ; i++) {
            taxa[i] = String.format("%d", i + 1);
        }
        for(int i = 0 ; i < nTaxa ; i++) {
            BitSet counter = new BitSet(nTaxa);
            counter.set(i, true);
            STITreeCluster tc = new STITreeCluster(taxa);
            tc.setCluster((BitSet) counter.clone());
            clusters.add(tc);
        }
        List<Tree> trees = new ArrayList<>();
        enumerateTreeTopologiesHelper(nTaxa, trees, clusters, new ArrayList<>(), taxa);
        return trees;
    }

    public static void enumerateNetworkTopologiesHelper(int nTaxa, int nReti, List<Network> networkList, Network currentNetwork) {

        if(currentNetwork.getReticulationCount() == nReti) {
            Networks.autoLabelNodes(currentNetwork);
            if(Networks.hasParallelEdges(currentNetwork)) {
                return;
            }

            for(int i = 0 ; i < networkList.size() ; i++) {
                if(Networks.hasTheSameTopology(currentNetwork, networkList.get(i))) {
                    return;
                }
            }
            networkList.add(currentNetwork.clone());
            return;
        }

        Network net = currentNetwork.clone();
        List<Tuple<NetNode, NetNode>> allEdges = Networks.getAllEdges(net);

        allEdges.add(new Tuple<>(net.getRoot(),null));

        for(int i = 0 ; i < allEdges.size() ; i++) {
            for(int j = 0 ; j < allEdges.size() ; j++) {
                //if(i == j) continue;

                Tuple<NetNode, NetNode> edge1 = allEdges.get(i);
                Tuple<NetNode, NetNode> edge2 = allEdges.get(j);

                NetNode v3 = edge1.Item2;
                NetNode v4 = edge1.Item1;
                NetNode v5 = edge2.Item2;
                NetNode v6 = edge2.Item1;
                NetNode v1 = new BniNetNode<>();
                NetNode v2 = new BniNetNode<>();

                if(edge1 != edge2) {

                    if (v3 != null) v3.removeChild(v4);
                    if (v5 != null) v5.removeChild(v6);

                    if (v3 == null) {
                        net.resetRoot(v1);
                    } else {
                        v3.adoptChild(v1, NetNode.NO_DISTANCE);
                    }
                    v1.adoptChild(v4, NetNode.NO_DISTANCE);

                    if (v5 == null) {
                        net.resetRoot(v2);
                    } else {
                        v5.adoptChild(v2, NetNode.NO_DISTANCE);
                    }
                    v2.adoptChild(v6, NetNode.NO_DISTANCE);

                    v1.adoptChild(v2, NetNode.NO_DISTANCE);
                } else {
                    if(v3 != null) {
                        v3.removeChild(v4);
                        v3.adoptChild(v1, NetNode.NO_DISTANCE);
                    } else {
                        net.resetRoot(v1);
                    }
                    v2.adoptChild(v4, NetNode.NO_DISTANCE);
                    v1.adoptChild(v2, NetNode.NO_DISTANCE);
                    v1.adoptChild(v2, NetNode.NO_DISTANCE);

                }

                if (Networks.hasCycle(net) || !Networks.isDisconnectedNetwork(net, null)) {


                } else {
                    enumerateNetworkTopologiesHelper(nTaxa, nReti, networkList, net);
                }

                if(edge1 != edge2) {
                    if (v5 == null) {
                        net.resetRoot(v6);
                    } else {
                        v5.removeChild(v2);
                    }
                    v2.removeChild(v6);

                    if (v3 == null) {
                        net.resetRoot(v4);
                    } else {
                        v3.removeChild(v1);
                    }
                    v1.removeChild(v4);

                    v1.removeChild(v2);

                    if (v5 != null) v5.adoptChild(v6, NetNode.NO_DISTANCE);
                    if (v3 != null) v3.adoptChild(v4, NetNode.NO_DISTANCE);
                } else {
                    v1.removeChild(v2);
                    v1.removeChild(v2);
                    v2.removeChild(v4);
                    if(v3 != null) {
                        v3.removeChild(v1);
                        v3.adoptChild(v4, NetNode.NO_DISTANCE);
                    } else {
                        net.resetRoot(v4);
                    }
                }

                if(!Networks.hasTheSameTopology(net, currentNetwork)) {
                    throw new RuntimeException("Error");
                }
            }
        }
    }

    public static List<Network> enumerateNetworkTopologies(int nTaxa, int nReti) {
        List<Tree> trees = enumerateTreeTopologies(nTaxa);

        List<Network> netReti0 = new ArrayList<>();
        for(Tree tree : trees) {
            netReti0.add(Networks.readNetwork(tree.toNewick()));

        }

        List<Network> results = new ArrayList<>();
        for(Network net : netReti0) {
            enumerateNetworkTopologiesHelper(nTaxa, nReti, results, net);
        }

        return results;
    }

    public static int getNetworkTopologiesCount(int nTaxa, int nReti) {
        return enumerateNetworkTopologies(nTaxa, nReti).size();
    }

    public List<Network> readAllNetwork(String filepath) {
        List<Network> result = new ArrayList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(filepath));
            String s;
            boolean begin = false;

            while((s = br.readLine().trim()) != null) {
                result.add(Networks.readNetwork(s));
            }

            br.close();
        } catch (Exception ex) {
            result = null;
        }
        return result;
    }

    public static void main(String[] args) {
        /*List<Tree> trees = enumerateTreeTopologies(6);
        System.out.println(trees.size());
        for(int i = 0 ; i < trees.size() ; i++) {
            System.out.println(trees.get(i));
        }*/

        List<Network> networks = enumerateNetworkTopologies(5, 1);
        System.out.println(networks.size());
        for(int i = 0 ; i < networks.size() ; i++) {
            System.out.println(networks.get(i));
        }
    }
}
