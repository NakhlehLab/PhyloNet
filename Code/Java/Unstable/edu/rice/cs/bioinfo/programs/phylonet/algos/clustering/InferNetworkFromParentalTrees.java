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

    public NetNode<Object> getMRCA(Network<Object> net, List<String> leaves) {

        Map<NetNode<Object>, Integer> marks = new HashMap<>();

        int OUT = 1;
        int IN = 2;

        for(NetNode<Object> node : net.bfs()) {
            marks.put(node, OUT);
        }

        marks.put(net.getRoot(), IN);

        for(NetNode<Object> leaf : net.getLeaves()) {
            if(leaves.contains(leaf.getName())) {
                NetNode<Object> node = leaf;
                while(!node.isRoot()) {
                    marks.put(node, IN);
                    for(NetNode<Object> parent : node.getParents()) {
                        node = parent;
                        break;
                    }
                }
            }
        }

        NetNode<Object> mrcaRoot = net.getRoot();
        while(!mrcaRoot.isLeaf()) {
            int count = 0;
            NetNode<Object> inChild = null;
            for(NetNode<Object> child : mrcaRoot.getChildren()) {
                if(marks.get(child) == IN) {
                    count++;
                    inChild = child;
                }
            }
            if(count == 1)
                mrcaRoot = inChild;
            else
                break;
        }

        return mrcaRoot;
    }

    public int getDeep(Network<Object> net) {

        Map<NetNode<Object>, Integer> marks = new HashMap<>();
        for(NetNode<Object> node : net.bfs()) {
            marks.put(node, -1);
        }

        NetNode<Object> mrca = net.getRoot();
        while(mrca.getChildCount() == 1) {
            mrca = mrca.getChildren().iterator().next();
        }

        marks.put(mrca, 0);



        int count = Integer.MAX_VALUE;
        for(NetNode<Object> node : net.bfs()) {
            if(marks.get(node) >= 0) {
                for(NetNode<Object> child : node.getChildren()) {
                    marks.put(child, marks.get(node) + 1);
                }
            }
        }

        for(NetNode<Object> leaf : net.getLeaves()) {
            if(marks.get(leaf) > count)
                count = marks.get(leaf);
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

        Map<Set<String>, Integer> numberMountPointsPerGroup = new HashMap<>();
        Map<MountPoint, List<List<Set<String>>>> mountPointsWithLeaves = new HashMap<>();
        Map<MountPoint, List<Set<String>>> orderOnMountPoint = new HashMap<>();
        List<MutableTuple<Set<String>, Set<MountPoint>>> groupedInLeaves = new ArrayList<>();

        for(int i = 0 ; i < parentalTrees.size() ; i++) {
            Map<MountPoint, List<Set<String>>> currentMountPointsWithLeaves = new HashMap<>();
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

            //Find all group in current parental tree
            for(NetNode<Object> outNode : net.bfs()) {
                Set<String> currentGroup = new HashSet<>();
                Set<MountPoint> currentMountPoints = new HashSet<>();

                if(marks.get(outNode) == OUT) {
                    NetNode<Object> inNode = null;
                    for(NetNode<Object> child : outNode.getChildren()) {
                        if(marks.get(child) == IN) {
                            inNode = child;
                            break;
                        }
                    }
                    if(inNode == null)
                        continue;

                    //Find all leaves in current group in current parental tree
                    Queue<NetNode<Object>> queue = new LinkedList<>();
                    queue.offer(inNode);
                    marks.put(inNode, DONE);

                    while(queue.size() > 0) {
                        NetNode<Object> node = queue.poll();
                        if(node.isLeaf()) {
                            currentGroup.add(node.getName());
                        }
                        for(NetNode<Object> child : node.getChildren()) {
                            queue.offer(child);
                            marks.put(child, DONE);
                        }
                    }

                    //Find Mount Point of current group
                    int count = 1;
                    NetNode<Object> outLeaf = inNode.getParents().iterator().next();
                    while(!outLeaf.isLeaf()) {
                        int cc = 0;
                        for(NetNode<Object> son : outLeaf.getChildren()) {
                            if(marks.get(son) == OUT) {
                                outLeaf = son;
                                cc++;
                            }
                        }
                        //if out-degree > 1 in MAST
                        if(cc > 1)
                            count++;
                    }
                    MountPoint currentMountPoint = new MountPoint(mast, outLeaf.getName(), count);
                    currentMountPoints.add(currentMountPoint);

                    //Find corresponding Order of current Mount Point
                    List<Set<String>> order = null;
                    for(MountPoint mountPoint : currentMountPointsWithLeaves.keySet()) {
                        if(mountPoint.isCongruence(currentMountPoint)) {
                            order = currentMountPointsWithLeaves.get(mountPoint);
                            break;
                        }
                    }

                    if(order == null) {
                        currentMountPointsWithLeaves.put(currentMountPoint, new ArrayList<>());
                        order = currentMountPointsWithLeaves.get(currentMountPoint);
                    }

                    order.add(new HashSet<String>(currentGroup));

                    //register current Group
                    boolean inserted = false;
                    //combine sets of leaves with intersection
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

            //combine current Order with Mount Point in current parental tree with global one
            for(MountPoint mountPoint1 : currentMountPointsWithLeaves.keySet()) {
                boolean found = false;
                for(MountPoint mountPoint2 : mountPointsWithLeaves.keySet()) {
                    if(mountPoint1.isCongruence(mountPoint2)) {
                        found = true;
                        mountPointsWithLeaves.get(mountPoint1).add(currentMountPointsWithLeaves.get(mountPoint1));
                        break;
                    }
                }

                if(!found) {
                    mountPointsWithLeaves.put(mountPoint1, new ArrayList<>());
                    mountPointsWithLeaves.get(mountPoint1).add(currentMountPointsWithLeaves.get(mountPoint1));
                }
            }


        }

        Map<String, Set<String>> groupFinder = new HashMap<>();
        boolean isLevel1 = true;

        //build Group finder
        for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {

            for(String leaf : tuple.Item1) {
                groupFinder.put(leaf, tuple.Item1);
            }

        }

        Map<Set<String>, Set<Set<String>>> subgroups = new IdentityHashMap<>();
        for(MountPoint mountPoint : mountPointsWithLeaves.keySet()) {
            //get number of SubGroup of one Group shown in current Mount Point
            //amount all Orders, keep the largest number
            Map<Set<String>, Integer> numberShown = new IdentityHashMap<>();
            for(List<Set<String>> order : mountPointsWithLeaves.get(mountPoint)) {
                Map<Set<String>, Integer> currentNumberShown = new IdentityHashMap<>();
                for(Set<String> set : order) {
                    Set<String> groupBelong = groupFinder.get(set.iterator().next());
                    if(currentNumberShown.get(groupBelong) == null)
                        currentNumberShown.put(groupBelong, 0);
                    currentNumberShown.put(groupBelong, currentNumberShown.get(groupBelong) + 1);
                }
                for(Set<String> group : currentNumberShown.keySet()) {
                    if(numberShown.get(group) == null)
                        numberShown.put(group, 0);
                    if(numberShown.get(group) < currentNumberShown.get(group))
                        numberShown.put(group, currentNumberShown.get(group));
                }

            }

            //initialize the SubGroups of one Group shown in current Mount Point
            //paritalOrder is the graph of Orders
            Map<Set<String>, Map<Set<String>, Integer>> partialOrder = new IdentityHashMap<>();
            Map<Set<String>, List<Set<String>>> subGroupOnMountPoint = new IdentityHashMap<>();
            for(Set<String> group : numberShown.keySet()) {
                subGroupOnMountPoint.put(group, new ArrayList<>());
                for(int i = 0 ; i < numberShown.get(group) ; i++) {
                    subGroupOnMountPoint.get(group).add(new HashSet<>());
                    partialOrder.put(subGroupOnMountPoint.get(group).get(i), new IdentityHashMap<>());
                }

            }

            //get all SubGroups of one Group shown in current Mount Point
            for(List<Set<String>> order : mountPointsWithLeaves.get(mountPoint)) {
                Map<Set<String>, Integer> currentNumberShown = new IdentityHashMap<>();
                Map<Set<String>, List<Set<String>>> curSubGroupOnMountPoint = new IdentityHashMap<>();

                //get all SubGroups of one Group shown in current Order in current Mount Point
                for(Set<String> set : order) {
                    Set<String> groupBelong = groupFinder.get(set.iterator().next());
                    if(currentNumberShown.get(groupBelong) == null)
                        currentNumberShown.put(groupBelong, 0);
                    currentNumberShown.put(groupBelong, currentNumberShown.get(groupBelong) + 1);
                    if(curSubGroupOnMountPoint.get(groupBelong) == null)
                        curSubGroupOnMountPoint.put(groupBelong, new ArrayList<>());
                    curSubGroupOnMountPoint.get(groupBelong).add(set);
                }

                //when number of SubGroup of one Group is correct in this Order -> update
                for(Set<String> group : numberShown.keySet()) {
                    if(numberShown.get(group).equals(currentNumberShown.get(group))) {
                        for(int i = 0 ; i < numberShown.get(group) ; i++)
                            subGroupOnMountPoint.get(group).get(i).addAll(curSubGroupOnMountPoint.get(group).get(i));
                    }

                }

                //when number of SubGroup of one Group is correct in this Order -> put into graph of Orders
                Map<Set<String>, Integer> currentScanNumberShown = new IdentityHashMap<>();
                Set<String> prevset = null;
                for(Set<String> set : order) {
                    Set<String> groupBelong = groupFinder.get(set.iterator().next());
                    if(currentScanNumberShown.get(groupBelong) == null)
                        currentScanNumberShown.put(groupBelong, 0);
                    Set<String> curset = subGroupOnMountPoint.get(groupBelong).get(currentScanNumberShown.get(groupBelong));
                    currentScanNumberShown.put(groupBelong, currentScanNumberShown.get(groupBelong) + 1);

                    if (numberShown.get(groupBelong).equals(currentNumberShown.get(groupBelong))) {

                        if(prevset != null) {
                            partialOrder.get(prevset).put(curset, 1);
                        }
                        prevset = curset;
                    }
                }




            }

            for(Set<String> group : subGroupOnMountPoint.keySet()) {
                if(subgroups.get(group) == null)
                    subgroups.put(group, new HashSet<>());
                subgroups.get(group).addAll(subGroupOnMountPoint.get(group));
            }

            //do topological sorting to get the correct Order
            List<Set<String>> currentOrder = new ArrayList<>();
            int totalSubMountPoints = partialOrder.size();
            Map<Set<String>, Integer> indegree = new IdentityHashMap<>();

            for(Set<String> set : partialOrder.keySet()) {
                indegree.put(set, 0);
            }

            for(Set<String> set : partialOrder.keySet()) {
                for(Set<String> nextset : partialOrder.get(set).keySet()) {
                    indegree.put(nextset, indegree.get(nextset) + 1);
                }
            }

            int subMountPointsCount = totalSubMountPoints;
            while(subMountPointsCount > 0) {
                boolean found = false;

                for(Set<String> set : partialOrder.keySet()) {
                    if(indegree.get(set) == 0) {
                        found = true;
                        currentOrder.add(set);
                        for(Set<String> nextset : partialOrder.get(set).keySet()) {
                            indegree.put(nextset, indegree.get(nextset) - 1);
                        }
                        indegree.put(set, -1);
                        subMountPointsCount--;
                        break;
                    }
                }

                if(subMountPointsCount == 0) break;
                //when there is a conflict, do breaking
                if(!found) {
                    for(Set<String> set : partialOrder.keySet()) {
                        if(indegree.get(set) > 0) {
                            found = true;
                            Set<String> newset = new HashSet<>(set);
                            indegree.put(newset, 0);
                            partialOrder.put(newset, partialOrder.get(set));
                            partialOrder.get(set).clear();
                            Set<String> groupBelong = groupFinder.get(set.iterator().next());
                            numberShown.put(groupBelong, numberShown.get(groupBelong) + 1);
                            subMountPointsCount++;
                            totalSubMountPoints++;
                            break;
                        }
                    }
                }
            }

            for(Set<String> subgroup : currentOrder) {
                Set<String> groupBelong = groupFinder.get(subgroup.iterator().next());
                if(numberMountPointsPerGroup.get(groupBelong) == null)
                    numberMountPointsPerGroup.put(groupBelong, 0);
                numberMountPointsPerGroup.put(groupBelong, numberMountPointsPerGroup.get(groupBelong) + 1);
            }
            orderOnMountPoint.put(mountPoint, currentOrder);


        }

        //check if every Group is associated with two Mount Point
        for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {

            if(numberMountPointsPerGroup.get(tuple.Item1) > 2)
                isLevel1 = false;
        }

        if(!isLevel1) {
            System.out.println("I don't think it is a level-1 network.");
        }

        if(false) {
            //remove duplicated mount points associated with one set of leaves
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

            }

            //combine groups with identical set of mount points
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

            Map<NetNode<Object>, Integer> marks = new HashMap<>();
            int ORIGINAL = 0;
            int NEW = 1;
            for(NetNode<Object> node : network.bfs()) {
                marks.put(node, ORIGINAL);
            }

            for(MountPoint mountPoint : orderOnMountPoint.keySet()) {
                NetNode<Object> node = network.findNode(mountPoint._position.Item1);
                for(int i = 0 ; i < mountPoint._position.Item2 - 1; i++) {
                    for(NetNode<Object> parent : node.getParents()) {
                        node = parent;
                        break;
                    }
                }

                for(int i = 0 ; i < orderOnMountPoint.get(mountPoint).size() ; i++) {
                    NetNode<Object> newnode = (NetNode<Object>) new BniNetNode<Object>();
                    NetNode<Object> parent = node.getParents().iterator().next();
                    parent.removeChild(node);
                    parent.adoptChild(newnode, NetNode.NO_DISTANCE);
                    newnode.adoptChild(node, NetNode.NO_DISTANCE);
                    marks.put(newnode, NEW);
                }
            }


            for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {
                List<String> currentInLeaves = new ArrayList<>(tuple.Item1);

                Network<Object> reticulatePart = null;
                Network<Object> bestFairPart = null;
                double mostFreq = -1;
                int minDeep = Integer.MAX_VALUE;

                //find the best reticulated part of one Group
                for (int i = 0; i < parentalTrees.size(); i++) {
                    Network<Object> net = Networks.readNetwork(parentalTrees.get(i).get(0).Item1.toNewick());
                    int count = getSubnetworkRetain(net, currentInLeaves);
                    removeBinaryNodes(net);

                    if (count == 1) {
                        if (reticulatePart == null) {
                            reticulatePart = net;
                        }
                    } else {
                        int deep = 0;
                        for(Set<String> subgroup : subgroups.get(tuple.Item1)) {
                            Network<Object> subnet = Networks.readNetwork(net.toString());
                            getSubnetworkRetain(subnet, new ArrayList<>(subgroup));
                            deep += getDeep(subnet);
                        }
                        if (deep < minDeep) {
                            minDeep = deep;
                            bestFairPart = net;
                        }
                    }
                }

                if (reticulatePart == null)
                    reticulatePart = bestFairPart;

                //mount best reticulated part to MAST
                for(MountPoint mountPoint : tuple.Item2) {
                    NetNode<Object> node = network.findNode(mountPoint._position.Item1);
                    for(int i = 0 ; i < mountPoint._position.Item2 - 1; i++) {
                        do {
                            node = node.getParents().iterator().next();
                        } while(marks.get(node) == NEW);
                    }

                    for(Set<String> subgroup : orderOnMountPoint.get(mountPoint)) {
                        if(tuple.Item1.containsAll(subgroup)) {
                            int subMountPoint = orderOnMountPoint.get(mountPoint).indexOf(subgroup);
                            subMountPoint = orderOnMountPoint.get(mountPoint).size() - subMountPoint;
                            NetNode<Object> connectNode = node;
                            for(int i = 0 ; i < subMountPoint ; i++)
                                connectNode = connectNode.getParents().iterator().next();

                            //TODO: Maybe bug when encounter partial SubGroup first

                            if(tuple.Item1.equals(subgroup) && reticulatePart.getRoot() != network.getRoot()) {
                                connectNode.adoptChild(reticulatePart.getRoot(), NetNode.NO_DISTANCE);
                            } else {

                                NetNode<Object> newConnectNode = new BniNetNode<>();
                                NetNode<Object> reticulateNode = getMRCA(reticulatePart, new ArrayList<>(subgroup));
                                if (reticulateNode.isRoot()) {
                                    newConnectNode.adoptChild(reticulateNode, NetNode.NO_DISTANCE);
                                    reticulatePart.resetRoot(newConnectNode);
                                } else {
                                    NetNode<Object> parent = reticulateNode.getParents().iterator().next();
                                    parent.removeChild(reticulateNode);
                                    parent.adoptChild(newConnectNode, NetNode.NO_DISTANCE);
                                    newConnectNode.adoptChild(reticulateNode, NetNode.NO_DISTANCE);
                                }
                                connectNode.adoptChild(newConnectNode, NetNode.NO_DISTANCE);
                            }
                        }
                    }
                }
            }
        }

        removeBinaryNodes(network);
        return network;
    }
}
