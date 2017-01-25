package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by dw20 on 5/31/15.
 */
public class NetworkProposal {

    private double[] _weights;
    private Random _rand;
    private int _maxReti;

    private String _operationName;
    private NetworkGTTOperation _op;

    public NetworkProposal(double[] weights, int maxNumReticulations, long seed) {
        _maxReti = maxNumReticulations;
        _rand = new Random(seed);
        _weights = weights;
    }

    public double propose(Network network) {
        return propose(network, null);
    }

        /**
         * Propose a new network
         * @return ln(hastings_ratio)
         */
    public double propose(Network network, Map<String, List<String>> net2treeMap) {
        List<Tuple<NetNode, NetNode>> allEdges = new ArrayList<Tuple<NetNode, NetNode>>();
        List<Tuple<NetNode, NetNode>> edgesNeedBrlens = new ArrayList<Tuple<NetNode, NetNode>>();
        List<Tuple<NetNode, NetNode>> allReticulationEdges = new ArrayList<Tuple<NetNode, NetNode>>();
        List<Tuple<NetNode, NetNode>> removableReticulationEdges = new ArrayList<>();
        Set<NetNode> leafNodes = new HashSet<NetNode>();
        int reticulations = getNetworkInfo(network, net2treeMap, allEdges, edgesNeedBrlens,
                allReticulationEdges, removableReticulationEdges, leafNodes);


        double kappa = _rand.nextDouble();
        /* dimension changing */
        if(kappa < _weights[0]) {
            if(reticulations == 0 || _rand.nextDouble() < _weights[1]) {
                double d = (reticulations == 0) ? (1.0 - _weights[1]) : (1.0 - _weights[1]) / _weights[1];
                return Math.log(d) + proposeAddReticulation(network, allEdges); // done
            } else {
                double d = (reticulations == 1) ? (1.0 - _weights[1]) : (1.0 - _weights[1]) / _weights[1];
                return -Math.log(d) + proposeDeleteReticulation(network, removableReticulationEdges); // done
            }
        } else {
            /* non-topology changing */
            double omega = _rand.nextDouble();
            if(omega < _weights[2]) {
                if(reticulations == 0 || _rand.nextDouble() < _weights[3]) {
                    return proposeChangeLength(network, allEdges); //edgesNeedBrlens);
                } else {
                    return proposeChangeProbability(network, allReticulationEdges); // done
                }
            } else {
                /* topology changing */
                double epsilon = _rand.nextDouble();
                if(removableReticulationEdges.size() == 0 || epsilon < _weights[4]) {
                    return proposeMoveTail(network, allEdges);
                } else if(epsilon < _weights[5]) {
                    return proposeMoveHead(network, removableReticulationEdges, allEdges);
                } else {
                    return proposeFlipReticulation(network, removableReticulationEdges);
                }
            }
        }
    }

    /**
     * Undo the operation
     */
    public void undo() {
        _op.undoOperation();
    }

    /**
     * Gets the name of operation
     * @return name
     */
    public String getOperationName() {
        return _operationName;
    }

    /**
     * targetEdge
     */
    private double proposeFlipReticulation(Network network,
                                           List<Tuple<NetNode, NetNode>> removableReticulationEdges) {
        _operationName = "Flip-Reticulation";
        int tar = _rand.nextInt(removableReticulationEdges.size());
        Tuple<NetNode, NetNode> tarEdge = removableReticulationEdges.get(tar);

        _op = new FlipReticulation();
        _op.setParameters(network, tarEdge, -1, -1, null, null, null, null, null, null);
        if(_op.performOperation()) return _op.getLogHR();
        return Double.MIN_VALUE;
    }

    private double proposeMoveHead(Network network, List<Tuple<NetNode, NetNode>> removableReticulationEdges,
                                   List<Tuple<NetNode, NetNode>> allEdges ) {
        _operationName = "Move-Head";
        int tar = _rand.nextInt(removableReticulationEdges.size());
        Tuple<NetNode, NetNode> tarEdge = removableReticulationEdges.get(tar);
        Tuple<NetNode, NetNode> destEdge = allEdges.get(_rand.nextInt(allEdges.size()));
        while(destEdge.Item1 == tarEdge.Item1 && destEdge.Item2 == tarEdge.Item2) {
            destEdge = allEdges.get(_rand.nextInt(allEdges.size()));
        }
        if(tarEdge == null
                || destEdge.Item2.equals(tarEdge.Item1) || destEdge.Item1.equals(tarEdge.Item1)
                || destEdge.Item2.equals(tarEdge.Item2) || destEdge.Item1.equals(tarEdge.Item2)){
            return Double.MIN_VALUE;
        }
        _op = new MoveHead();
        _op.setParameters(network, tarEdge, null, destEdge);
        if(_op.performOperation()) return _op.getLogHR();
        return Double.MIN_VALUE;
    }

    private double proposeMoveTail(Network network, List<Tuple<NetNode, NetNode>> allEdges) {
        _operationName = "Move-Tail";
        int tar = _rand.nextInt(allEdges.size());
        Tuple<NetNode, NetNode> tarEdge = allEdges.get(tar);
        while(tarEdge.Item1.isNetworkNode()) {
            tar = _rand.nextInt(allEdges.size());
            tarEdge = allEdges.get(tar);
        }

        int dest = _rand.nextInt(allEdges.size());
        Tuple<NetNode, NetNode> destEdge = dest == tar ? null : allEdges.get(dest);

        if(tarEdge.Item1.isRoot()){
            for(Object anotherChild: tarEdge.Item1.getChildren()){
                if(anotherChild!=tarEdge.Item2 && ((NetNode)anotherChild).isNetworkNode()){
                    return Double.MIN_VALUE;
                }
            }
        }
        if(tarEdge != null && destEdge != null &&
                (tarEdge.Item1.equals(destEdge.Item1) || tarEdge.Item1.equals(destEdge.Item2) ||
                        tarEdge.Item2.equals(destEdge.Item1) || tarEdge.Item2.equals(destEdge.Item2))) {
            return Double.MIN_VALUE;
        }
        _op = new MoveTail();
        _op.setParameters(network, tarEdge, null, destEdge);
        if(_op.performOperation()) return _op.getLogHR();
        return Double.MIN_VALUE;
    }

    private double proposeAddReticulation(Network network, List<Tuple<NetNode, NetNode>> allEdges) {
        _operationName = "Add-Reticulation";
        if(network.getReticulationCount() == _maxReti) return Double.MIN_VALUE;
        Tuple<NetNode, NetNode> srcEdge, destEdge;
        int n = allEdges.size();
        if(n < 2) return Double.MIN_VALUE;
        int r1 = _rand.nextInt(n), r2 = _rand.nextInt(n);
        while(r2 == r1) r2 = _rand.nextInt(n);
        srcEdge = allEdges.get(r1);
        destEdge = allEdges.get(r2);
        _op = new AddReticulation();
        _op.setParameters(network, null, srcEdge, destEdge);
        if(_op.performOperation()) {
            double hr = _op.getLogHR();
            return hr;
        }
        return Double.MIN_VALUE;
    }

    private double proposeDeleteReticulation(Network network,
                                             List<Tuple<NetNode, NetNode>> removableReticulationEdges) {
        _operationName = "Delete-Reticulation";

        if(removableReticulationEdges.size() == 0) return Double.MIN_VALUE;

        int tar = _rand.nextInt(removableReticulationEdges.size());
        Tuple<NetNode, NetNode> tarEdge = removableReticulationEdges.get(tar);
        if(!removableReticulationEdges.contains(tarEdge)) return Double.MIN_VALUE;
        _op = new DeleteReticulation();
        _op.setParameters(network, tarEdge, null, null);
        if(_op.performOperation()) {
            double hr = _op.getLogHR();
            return hr;
        }
        return Double.MIN_VALUE;
    }

    private double proposeChangeLength(Network network, List<Tuple<NetNode, NetNode>> edgesNeedBrlens) {
        _operationName = "Change-Length";
        _op = new ChangeLength();
        Tuple<NetNode, NetNode> tarEdge = edgesNeedBrlens.get(_rand.nextInt(edgesNeedBrlens.size()));
        _op.setParameters(network, tarEdge, -1, -1, null, null, null, null, null, null);
        if(_op.performOperation()) return _op.getLogHR();
        return Double.MIN_VALUE;
    }

    private double proposeChangeProbability(Network network, List<Tuple<NetNode, NetNode>> reticulationEdges) {
        _operationName = "Change-Probability";
        _op = new ChangeProbability();
        Tuple<NetNode, NetNode> tarEdge = reticulationEdges.get(_rand.nextInt(reticulationEdges.size()));
        _op.setParameters(network, tarEdge, -1, -1, null, null, null, null, null, null);
        if(_op.performOperation()) return _op.getLogHR();
        return Double.MIN_VALUE;
    }

    private int getNetworkInfo(Network network, Map<String, List<String>> net2treeMap,
                               List<Tuple<NetNode, NetNode>> allEdges,
                                List<Tuple<NetNode, NetNode>> edgesNeedBrlens,
                                List<Tuple<NetNode, NetNode>> allReticulationEdges,
                                List<Tuple<NetNode, NetNode>> removableReticulationEdges,
                                Set<NetNode> leafNodes){
        int numReticulations = 0;
        Map<NetNode,Set<String>> node2leaves = new HashMap<>();
        for(Object nodeO: Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeO;
            Set<String> leaves = new HashSet<>();
            if(node.isLeaf()){
                leafNodes.add(node);
                if(net2treeMap == null) {
                    leaves.add(node.getName());
                } else {
                    leaves.addAll(net2treeMap.get(node.getName()));
                }
            } else if(node.isNetworkNode()){
                numReticulations++;
            }
            for(Object childO: node.getChildren()){
                NetNode childNode = (NetNode)childO;
                Set<String> childLeaves = node2leaves.get(childNode);
                leaves.addAll(childLeaves);
                Tuple<NetNode,NetNode> edge = new Tuple<>(node, childNode);
                allEdges.add(edge);
                if(childLeaves.size()!=1){
                    edgesNeedBrlens.add(edge);
                }
                if(childNode.isNetworkNode()) {
                    if (node.isTreeNode()) {
                        removableReticulationEdges.add(edge);
                    }
                    allReticulationEdges.add(edge);
                }

            }
            node2leaves.put(node,leaves);
        }
        return numReticulations;
    }
}
