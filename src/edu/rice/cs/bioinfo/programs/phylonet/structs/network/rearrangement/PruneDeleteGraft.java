package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of NetworkRearrangementOperation.
 * It mutates a network by re-grafting the the sub-network rooted at <code>targetNode</code> to a random edge
 */
public class PruneDeleteGraft extends NetworkRearrangementOperation {

    /**
     * This function is to set parameters
     *
     * @param network   the network
     * @param targetNode    a node in <code>network</code>
     */
    public void setParameters(Network network, Tuple<NetNode, NetNode> targetNode) {
        _network = network;
        _targetEdge = targetNode;
    }


    /**
     * This function is to undo the operation
     */
    public void undoOperation() {
    }


    /**
     * This is the main function for mutating the network
     *
     * @return  whether the operation is performed successfully
     *          returns false when the resulting network is invalid
     */
    public boolean performOperation() {
        List<NetNode> parents = new ArrayList<>();
        for (Object parent : _targetEdge.Item1.getParents()) {
            parents.add((NetNode) parent);
        }
        for (NetNode parent : parents) {
            parent.removeChild(_targetEdge.Item1);
        }
        for (Object child : _targetEdge.Item1.getChildren()) {
            if (child != _targetEdge.Item2) {
                _targetEdge.Item1.removeChild((NetNode) child);
            }
        }
        Network tempNet = new BniNetwork((BniNetNode) _targetEdge.Item1);
        for (Object leaf : tempNet.getLeaves()) {
            deleteNode(_network, _network.findNode(((NetNode) leaf).getName()));
        }
        graftToRandomEdge(_network, _targetEdge.Item1);

        return true;
    }


    /**
     * This function is for adjusting the network (i.e., removing binary nodes)
     */
    private void adjustNetwork(Network network) {
        Set<NetNode> leaves = new HashSet<>();
        for (Object nodeO : network.getLeaves()) {
            NetNode node = (NetNode) nodeO;
            if (!node.isNetworkNode()) {
                leaves.add(node);
            }
        }

        boolean update;
        do {
            update = false;
            for (Object nodeO : Networks.postTraversal(network)) {
                NetNode node = (NetNode) nodeO;
                if (node.isLeaf() && !leaves.contains(node)) {
                    node.removeItself();
                    update = true;
                    break;
                }
            }
            Networks.removeBinaryNodes(network);
        } while (update);

        NetNode root = network.getRoot();
        while (root.getChildCount() == 1) {
            NetNode child = (NetNode) root.getChildren().iterator().next();
            if (!child.isLeaf()) {
                root.removeChild(child);
                root = child;
            } else {
                break;
            }
        }

        network.resetRoot(root);
        for (Object nodeO : network.getNetworkNodes()) {
            NetNode node = (NetNode) nodeO;
            if (node.isTreeNode() && !node.isRoot()) {
                NetNode parentNode = (NetNode) (node.getParents().iterator().next());
                node.setParentProbability(parentNode, NetNode.NO_PROBABILITY);
            }
        }

    }


    /**
     * This function is to graft the node <code>toGraft</code> to an random edge of network <code>network</code>
     */
    private void graftToRandomEdge(Network network, NetNode toGraft) {
        List<Tuple<NetNode, NetNode>> allEdges = new ArrayList<>();
        for (Object nodeO : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeO;
            for (Object childO : node.getChildren()) {
                NetNode childNode = (NetNode) childO;
                Tuple<NetNode, NetNode> edge = new Tuple<>(node, childNode);
                allEdges.add(edge);
            }
        }
        if (allEdges.size() == 1) {
            Tuple<NetNode, NetNode> onlyEdge = allEdges.get(0);
            toGraft.adoptChild(onlyEdge.Item2, onlyEdge.Item2.getParentDistance(onlyEdge.Item1));
            onlyEdge.Item1.removeChild(onlyEdge.Item2);
            network.resetRoot(toGraft);
        } else {
            Random random = new Random();
            Tuple<NetNode, NetNode> sourceEdge = allEdges.get(random.nextInt(allEdges.size()));

            double[] brlens = new double[2];
            double[] inheriProbs = new double[2];
            randomlyPartitionAnEdge(sourceEdge, brlens, inheriProbs);
            addNodeToAnEdge(toGraft, sourceEdge, brlens, inheriProbs);
        }
    }


    /**
     * This function is to delete node <code>node</code> in network <code>network</code>
     */
    private void deleteNode(Network network, NetNode node) {
        do {
            NetNode parentNode = (NetNode) (node.getParents().iterator().next());
            if (parentNode.isNetworkNode()) {
                parentNode.removeChild(node);
                break;
            } else {
                Tuple<NetNode, NetNode> temp = findParentAndAnotherChild(parentNode, node);
                if (temp.Item2.isNetworkNode()) {
                    parentNode.removeChild(temp.Item2);
                    parentNode.removeChild(node);
                } else {
                    parentNode.removeChild(node);
                    break;
                }
            }
            node = parentNode;
        } while (true);

        adjustNetwork(network);
    }


}
