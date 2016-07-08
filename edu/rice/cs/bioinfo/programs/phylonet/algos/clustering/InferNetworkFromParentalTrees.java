package edu.rice.cs.bioinfo.programs.phylonet.algos.clustering;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.mast.SteelWarnowMAST;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 7/4/16
 * Time: 4:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferNetworkFromParentalTrees {
    private void removeBinaryNodes(Network<Object> net)
    {
        // Find all binary nodes.
        List<NetNode<Object>> binaryNodes = new LinkedList<NetNode<Object>>();
        for (NetNode<Object> node : net.bfs()) {
            if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                binaryNodes.add(node);
            }
        }

        // Remove them.
        for (NetNode<Object> node : binaryNodes) {
            NetNode<Object> child = node.getChildren().iterator().next();	// Node's only child.
            if(child.getIndeg() != 1){
                continue;
            }
            NetNode<Object> parent = node.getParents().iterator().next();	// Node's only parent.
            double distance = node.getParentDistance(parent) + child.getParentDistance(node);
            double gamma1 = node.getParentProbability(parent)==NetNode.NO_PROBABILITY?1.0:node.getParentProbability(parent);
            double gamma2 = child.getParentProbability(node)==NetNode.NO_PROBABILITY?1.0:child.getParentProbability(node);
            double gamma =  gamma1 * gamma2;
            parent.removeChild(node);
            node.removeChild(child);
            parent.adoptChild(child, distance);
            child.setParentProbability(parent, gamma);
        }
    }

    public int getSubnetworkRetain(Network<Object> net, List<String> leaves) {

        Map<NetNode<Object>, Integer> marks = new HashMap<>();

        int OUT = 1;
        int IN = 2;

        for(NetNode<Object> node : net.bfs()) {
            marks.put(node, IN);
        }

        marks.put(net.getRoot(), OUT);

        for(NetNode<Object> leaf : net.getLeaves()) {
            if(!leaves.contains(leaf.getName())) {
                NetNode<Object> node = leaf;
                while(!node.isRoot()) {
                    marks.put(node, OUT);
                    for(NetNode<Object> parent : node.getParents()) {
                        node = parent;
                        break;
                    }
                }
            }
        }

        int count = 0;
        for(NetNode<Object> node : net.bfs()) {
            if(!node.isRoot()) {
                for (NetNode<Object> parent : node.getParents()) {
                    if(!marks.get(node).equals(marks.get(parent))) {
                        count++;
                    }
                }
            }
        }

        boolean done = false;
        while(!done) {
            done = true;
            for (NetNode<Object> node : net.dfs()) {
                if (!leaves.contains(node.getName()) && node.getChildCount() == 0) {
                    for (NetNode<Object> parent : node.getParents()) {
                        parent.removeChild(node);
                        done = false;
                    }
                }
            }
        }

        return count;
    }

    public List<NetNode> getConnectNodes(Network net, List<String> outLeaves) {
        List<NetNode> connectNodes = new ArrayList<>();
        Map<NetNode, Boolean> marks = new HashMap<>();

        return connectNodes;
    }

    class MountPoint {
        Tree _tree;
        MutableTuple<String, Integer> _position;

        MountPoint(Tree tree, String leaf, int up) {
            _tree = tree;
            _position = new MutableTuple<>(leaf, up);
        }

        @Override
        public boolean equals(Object candidate)
        {
            try
            {
                return equals((MountPoint) candidate);
            }
            catch (ClassCastException e)
            {
                return false;
            }
        }

        public boolean equals(MountPoint candidate)
        {
            if(candidate == null)
                return false;

            if(candidate._tree != _tree)
                return false;

            return candidate._position.equals(_position);
        }


        public boolean isCongruence(MountPoint candidate)
        {
            if(candidate == null)
                return false;

            if(candidate._tree != _tree)
                return false;

            TNode node1 = _tree.getNode(_position.Item1);
            for(int i = 0 ; i < _position.Item2 - 1 ; i++)
                node1 = node1.getParent();

            TNode node2 = _tree.getNode(candidate._position.Item1);
            for(int i = 0 ; i < candidate._position.Item2 - 1 ; i++)
                node2 = node2.getParent();

            return node1 == node2;
        }

        @Override
        public int hashCode()
        {
            return _tree.hashCode() ^ _position.Item1.hashCode() ^ _position.Item2.hashCode();
        }
    }

    public Network<Object> inferNetwork(List<List<MutableTuple<Tree, Double>>> parentalTrees) {
        Network<Object> network = null;
        List<Tree> trees = new ArrayList<>();
        for(int i = 0 ; i < parentalTrees.size() ; i++) {
            trees.add(parentalTrees.get(i).get(0).Item1);
        }

        SteelWarnowMAST steelWarnowMAST = new SteelWarnowMAST();
        Tree mast = steelWarnowMAST.computeRMAST(trees);
        List<String> outLeaves = new ArrayList<>();
        for(String leaf : mast.getLeaves())
            outLeaves.add(leaf);

        System.out.println("MAST: " + mast.toNewick());

        network = Networks.readNetwork(mast.toNewick());

        List<MutableTuple<Set<String>, Set<MountPoint>>> groupedInLeaves = new ArrayList<>();

        for(int i = 0 ; i < parentalTrees.size() ; i++) {
            Network<Object> net = Networks.readNetwork(parentalTrees.get(i).get(0).Item1.toNewick());
            Map<NetNode<Object>, Integer> marks = new HashMap<>();

            int OUT = 1;
            int IN = 2;
            int DONE = 3;

            for(NetNode<Object> node : net.bfs()) {
                marks.put(node, IN);
            }

            marks.put(net.getRoot(), OUT);

            for(NetNode<Object> leaf : net.getLeaves()) {
                if(outLeaves.contains(leaf.getName())) {
                    NetNode<Object> node = leaf;
                    while(!node.isRoot()) {
                        marks.put(node, OUT);
                        for(NetNode<Object> parent : node.getParents()) {
                            node = parent;
                            break;
                        }
                    }
                }
            }

            for(NetNode<Object> leaf : net.getLeaves()) {
                Set<String> currentGroup = new HashSet<>();
                Set<MountPoint> currentMountPoints = new HashSet<>();

                if(!outLeaves.contains(leaf.getName()) && marks.get(leaf) == IN) {
                    NetNode<Object> node = leaf;

                    while(marks.get(node) == IN) {
                        for(NetNode<Object> parent : node.getParents()) {
                            node = parent;
                            break;
                        }
                    }

                    int count = 1;
                    NetNode<Object> outLeaf = node;
                    while(!outLeaf.isLeaf()) {
                        int cc = 0;
                        for(NetNode<Object> son : outLeaf.getChildren()) {
                            if(marks.get(son) == OUT) {
                                outLeaf = son;
                                cc++;
                            }
                        }
                        if(cc > 1)
                            count++;
                    }
                    currentMountPoints.add(new MountPoint(mast, outLeaf.getName(), count));

                    for(NetNode<Object> son : node.getChildren()) {
                        if(marks.get(son) == IN) {
                            node = son;
                            break;
                        }
                    }

                    Queue<NetNode<Object>> queue = new LinkedList<>();
                    queue.offer(node);
                    marks.put(node, DONE);

                    while(queue.size() > 0) {
                        node = queue.poll();
                        if(node.isLeaf()) {
                            currentGroup.add(node.getName());
                        }
                        for(NetNode<Object> son : node.getChildren()) {
                            queue.offer(son);
                            marks.put(son, DONE);
                        }
                    }
                }
                boolean inserted = false;
                for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {
                    Set<String> group = tuple.Item1;
                    Set<String> intersection = new HashSet<String>(group);
                    intersection.retainAll(currentGroup);
                    if(intersection.size() > 0) {
                        inserted = true;
                        group.addAll(currentGroup);
                        tuple.Item2.addAll(currentMountPoints);
                    }
                }

                if(inserted) {
                    boolean done = false;
                    while(!done) {
                        done = true;
                        for(MutableTuple<Set<String>, Set<MountPoint>> tuple1 : groupedInLeaves) {
                            Set<String> group1 = tuple1.Item1;
                            for(MutableTuple<Set<String>, Set<MountPoint>> tuple2 : groupedInLeaves) {
                                Set<String> group2 = tuple2.Item1;
                                if(group1 != group2) {
                                    Set<String> intersection = new HashSet<String>(group1);
                                    intersection.retainAll(group2);
                                    if(intersection.size() > 0) {
                                        group1.addAll(group2);
                                        tuple1.Item2.addAll(tuple2.Item2);
                                        groupedInLeaves.remove(tuple2);
                                        done = false;
                                        break;
                                    }
                                }
                            }
                            if(!done)
                                break;
                        }
                    }
                }

                if(!inserted && currentGroup.size() > 0)
                    groupedInLeaves.add(new MutableTuple<Set<String>, Set<MountPoint>>(currentGroup, currentMountPoints));
            }
        }

        boolean isLevel1 = true;
        for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {
            boolean done = false;
            while(!done) {
                done = true;
                for (MountPoint mountPoint1 : tuple.Item2) {
                    for (MountPoint mountPoint2 : tuple.Item2) {
                        if (mountPoint1 != mountPoint2 && mountPoint2.isCongruence(mountPoint1)) {
                            tuple.Item2.remove(mountPoint2);
                            done = false;
                            break;
                        }
                    }
                    if(!done)
                        break;
                }
            }
            if(tuple.Item2.size() != 2)
                isLevel1 = false;
        }

        if(isLevel1) {
            //combine groups with same mount points
            boolean combineDone = false;
            while(!combineDone) {
                combineDone = true;
                for (MutableTuple<Set<String>, Set<MountPoint>> tuple1 : groupedInLeaves) {
                    for (MutableTuple<Set<String>, Set<MountPoint>> tuple2 : groupedInLeaves) {
                        if(tuple1 != tuple2) {
                            boolean setCongruence = true;
                            for (MountPoint mountPoint1 : tuple1.Item2) {
                                boolean found = false;
                                for (MountPoint mountPoint2 : tuple2.Item2) {
                                    if (mountPoint1.isCongruence(mountPoint2)) {
                                        found = true;
                                        break;
                                    }
                                }
                                if (!found) {
                                    setCongruence = false;
                                    break;
                                }
                            }
                            if (setCongruence) {
                                combineDone = false;
                                tuple1.Item1.addAll(tuple2.Item1);
                                groupedInLeaves.remove(tuple2);
                                break;
                            }
                        }
                    }
                    if(!combineDone)
                        break;
                }
            }

            for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {
                List<String> currentInLeaves = new ArrayList<>(tuple.Item1);

                Network<Object> reticulatePart = null;
                Network<Object> mostFreqPart = null;
                double mostFreq = -1;

                for (int i = 0; i < parentalTrees.size(); i++) {
                    Network<Object> net = Networks.readNetwork(parentalTrees.get(i).get(0).Item1.toNewick());
                    int count = getSubnetworkRetain(net, currentInLeaves);
                    removeBinaryNodes(net);

                    if (count == 1) {
                        if (reticulatePart == null) {
                            reticulatePart = net;
                            removeBinaryNodes(reticulatePart);
                        }
                    } else {
                        if (mostFreq < parentalTrees.get(i).get(0).Item2) {
                            mostFreq = parentalTrees.get(i).get(0).Item2;
                            mostFreqPart = net;
                            removeBinaryNodes(mostFreqPart);
                        }
                    }
                }

                if (reticulatePart == null)
                    reticulatePart = mostFreqPart;

                List<NetNode<Object>> mountNodes = new ArrayList<>();

                for(MountPoint mountPoint : tuple.Item2) {
                    NetNode<Object> node = network.findNode(mountPoint._position.Item1);
                    for(int i = 0 ; i < mountPoint._position.Item2 - 1; i++) {
                        for(NetNode<Object> parent : node.getParents()) {
                            node = parent;
                            break;
                        }
                    }
                    mountNodes.add(node);
                }

                for(NetNode<Object> mountNode : mountNodes) {
                    NetNode<Object> newnode = (NetNode<Object>) new BniNetNode<Object>();
                    NetNode<Object> parent = mountNode.getParents().iterator().next();
                    parent.removeChild(mountNode);
                    parent.adoptChild(newnode, NetNode.NO_DISTANCE);
                    newnode.adoptChild(mountNode, NetNode.NO_DISTANCE);
                    newnode.adoptChild(reticulatePart.getRoot(), NetNode.NO_DISTANCE);
                }


            }
        } else {
            System.out.println("I think it is not a level-1 network.");
        }

        removeBinaryNodes(network);
        return network;
    }
}
