package edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 3/21/13
 * Time: 11:01 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkBootstrap<T> {

    public void computeSoftwiredNetworkBranchSupport(Network<T> targetNet, List<Network<T>> refNetList){
        clearNetworkBranchSupport(targetNet);

        Map<NetNode<T>, Set<NetworkCluster<T>>> tNode2clusters = generateSoftwiredClusters(targetNet);
        for(Network<T> refNet: refNetList){
            //Networks.removeAllParameters(refNet);
            //System.out.println(refNet);

            Map<NetNode<T>, Set<NetworkCluster<T>>> rNode2clusters = generateSoftwiredClusters(refNet);
            for(Map.Entry<NetNode<T>, Set<NetworkCluster<T>>> tEntry: tNode2clusters.entrySet()){
                NetNode<T> tNode = tEntry.getKey();
                if(tNode.isRoot() || tNode.isLeaf()){
                    continue;
                }
                Set<NetworkCluster<T>> targetNodeClusters = tEntry.getValue();
                int targetSize = targetNodeClusters.size();
                boolean found = false;
                for(Map.Entry<NetNode<T>, Set<NetworkCluster<T>>> rEntry: rNode2clusters.entrySet()){
                    NetNode<T> rNode = rEntry.getKey();
                    if((tNode.isNetworkNode() && !rNode.isNetworkNode()) || (tNode.isTreeNode() && !rNode.isTreeNode())){
                        continue;
                    }
                    Set<NetworkCluster<T>> refNodeClusters = rEntry.getValue();
                    if(targetSize == refNodeClusters.size()){
                        if(targetNodeClusters.equals(refNodeClusters)){
                            found = true;
                            break;
                        }
                    }
                }
                if(found){
                    for(NetNode<T> parent: tNode.getParents()){
                        tNode.setParentSupport(parent, tNode.getParentSupport(parent)+1);
                    }
                }
            }
            //System.out.println(Networks.network2string(targetNet));
        }

        normalizeNetworkBranchSupport(targetNet, refNetList.size());
    }

    private void clearNetworkBranchSupport(Network<T> network){
        for(NetNode<T> node: ((BniNetwork<T>)network).bfs()){
            for(NetNode<T> child: node.getChildren()){
                if(child.isRoot() || child.isLeaf()){
                    continue;
                }
                child.setParentSupport(node, 0);
            }
        }
    }


    private void normalizeNetworkBranchSupport(Network<T> network, int total){
        for(NetNode<T> node: ((BniNetwork<T>)network).bfs()){
            for(NetNode<T> child: node.getChildren()){
                if(child.isRoot() || child.isLeaf()){
                    continue;
                }
                double support = child.getParentSupport(node)/total;
                child.setParentSupport(node, support);
            }
        }
    }


    private Map<NetNode<T>, Set<NetworkCluster<T>>> generateSoftwiredClusters(Network<T> network){
        Map<NetNode<T>, Set<NetworkCluster<T>>> node2clusters = new Hashtable<NetNode<T>, Set<NetworkCluster<T>>>();
        for (NetworkTree<T> nt : Networks.getTrees(network)) {
            for(Map.Entry<NetNode<T>, NetworkCluster<T>> entry: nt.generateNodeClusters().entrySet()){
                NetNode<T> node = entry.getKey();
                NetworkCluster<T> nc = entry.getValue();
                if (!nc.isEmpty()) {
                    Set<NetworkCluster<T>> clusters = node2clusters.get(node);
                    if(clusters == null){
                        clusters = new HashSet<NetworkCluster<T>>();
                        node2clusters.put(node, clusters);
                    }
                    clusters.add(nc);
                }
            }
        }

        return node2clusters;
    }


    public void computeHardwiredNetworkBranchSupport(Network<T> targetNet, List<Network<T>> refNetList){
        clearNetworkBranchSupport(targetNet);

        Map<NetNode<T>, Set<String>> tNode2clusters = generateHardwiredClusters(targetNet);
        for(Network<T> refNet: refNetList){
            Map<NetNode<T>, Set<String>> rNode2clusters = generateHardwiredClusters(refNet);
            for(Map.Entry<NetNode<T>, Set<String>> tEntry: tNode2clusters.entrySet()){
                NetNode<T> tNode = tEntry.getKey();
                if(tNode.isRoot() || tNode.isLeaf()){
                    continue;
                }
                if(rNode2clusters.containsValue(tEntry.getValue())){
                    for(NetNode<T> parent: tNode.getParents()){
                        tNode.setParentSupport(parent, tNode.getParentSupport(parent)+1);
                    }
                }
            }
            //System.out.println(Networks.network2string(targetNet));
        }

        normalizeNetworkBranchSupport(targetNet, refNetList.size());
    }


    private Map<NetNode<T>, Set<String>> generateHardwiredClusters(Network<T> network){
        Map<NetNode<T>, Set<String>> node2cluster = new HashMap<NetNode<T>, Set<String>>();
        for (NetNode<T> node: Networks.postTraversal(network)) {
            if(node.isRoot()){
                break;
            }
            Set<String> cluster = new HashSet<String>();
            if(node.isLeaf()){
                cluster.add(node.getName());
            }
            else{
                for(NetNode<T> child: node.getChildren()){
                    cluster.addAll(node2cluster.get(child));
                }
            }
            node2cluster.put(node, cluster);
        }

        return node2cluster;
    }

    private void clearReticulationBranchSupport(Network<T> network){
        for(NetNode<T> node: ((BniNetwork<T>)network).bfs()){
            for(NetNode<T> child: node.getChildren()){
                if(child.isNetworkNode()){
                    child.setParentSupport(node, 0);
                }
            }
        }
    }


    public void computeHardwiredReticulationBranchSupport(Network<T> targetNet, List<Network<T>> refNetList){
        clearReticulationBranchSupport(targetNet);

        List<Tuple<Tuple<NetNode<T>,NetNode<T>>,Tuple<Set<String>,Set<String>>>> tRetiEdge2clusters = new ArrayList<Tuple<Tuple<NetNode<T>, NetNode<T>>, Tuple<Set<String>, Set<String>>>>();
        generateReticulationBranchAndHardwiredClusters(targetNet, tRetiEdge2clusters);
        for(Network<T> refNet: refNetList){
            List<Tuple<Set<String>,Set<String>>> rClusterSet = new ArrayList<Tuple<Set<String>, Set<String>>>();
            generateHardwiredReticulationBranchClusters(refNet, rClusterSet);

            for(Tuple<Tuple<NetNode<T>,NetNode<T>>,Tuple<Set<String>,Set<String>>> tRetiEdge: tRetiEdge2clusters){
                if(rClusterSet.contains(tRetiEdge.Item2)){
                    NetNode<T> parentNode = tRetiEdge.Item1.Item1;
                    NetNode<T> childNode = tRetiEdge.Item1.Item2;
                    childNode.setParentSupport(parentNode, childNode.getParentSupport(parentNode)+1);
                }
            }
            //System.out.println(Networks.network2string(targetNet));
        }

        normalizeNetworkBranchSupport(targetNet, refNetList.size());
    }


    private void generateReticulationBranchAndHardwiredClusters(Network<T> network, List<Tuple<Tuple<NetNode<T>,NetNode<T>>,Tuple<Set<String>,Set<String>>>> retiEdge2clusters){
        Map<NetNode<T>, Set<String>> node2cluster = new HashMap<NetNode<T>, Set<String>>();
        //List<Tuple<Tuple<NetNode<T>,NetNode<T>>,Tuple<Set<String>,Set<String>>>> retiEdge2clusters = new ArrayList<Tuple<Tuple<NetNode<T>, NetNode<T>>, Tuple<Set<String>, Set<String>>>>();
        for (NetNode<T> node: Networks.postTraversal(network)) {
            if(node.isRoot()){
                break;
            }
            Set<String> cluster = new HashSet<String>();
            if(node.isLeaf()){
                cluster.add(node.getName());
            }
            else{
                for(NetNode<T> child: node.getChildren()){
                    Set<String> childCluster = node2cluster.get(child);
                    cluster.addAll(childCluster);
                    if(child.isNetworkNode()){
                        retiEdge2clusters.add(new Tuple<Tuple<NetNode<T>,NetNode<T>>,Tuple<Set<String>,Set<String>>>(new Tuple<NetNode<T>,NetNode<T>>(node,child), new Tuple<Set<String>, Set<String>>(cluster, childCluster)));
                    }
                }
            }
            node2cluster.put(node, cluster);

        }

        //return retiEdge2clusters;
    }


    private void generateHardwiredReticulationBranchClusters(Network<T> network, List<Tuple<Set<String>,Set<String>>> retiEdgeClusters){
        Map<NetNode<T>, Set<String>> node2cluster = new HashMap<NetNode<T>, Set<String>>();
        for (NetNode<T> node: Networks.postTraversal(network)) {
            if(node.isRoot()){
                break;
            }
            Set<String> cluster = new HashSet<String>();
            if(node.isLeaf()){
                cluster.add(node.getName());
            }
            else{
                for(NetNode<T> child: node.getChildren()){
                    Set<String> childCluster = node2cluster.get(child);
                    cluster.addAll(childCluster);
                    if(child.isNetworkNode()){
                        retiEdgeClusters.add(new Tuple<Set<String>, Set<String>>(cluster, childCluster));
                    }
                }
            }
            node2cluster.put(node, cluster);

        }
    }

}
