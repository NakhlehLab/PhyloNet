package edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization;

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

    public Network computeNetworkBranchSupport(Network<T> targetNet, List<Network<T>> refNetList){
        clearNetworkBranchSupport(targetNet);

        Map<NetNode<T>, Set<NetworkCluster<T>>> tNode2clusters = generateNetworkClusters(targetNet);
        for(Network<T> refNet: refNetList){
            Map<NetNode<T>, Set<NetworkCluster<T>>> rNode2clusters = generateNetworkClusters(refNet);
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
        }

        normalizeNetworkBranchSupport(targetNet, refNetList.size());
        return targetNet;
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


    private Map<NetNode<T>, Set<NetworkCluster<T>>> generateNetworkClusters(Network<T> network){
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


}
