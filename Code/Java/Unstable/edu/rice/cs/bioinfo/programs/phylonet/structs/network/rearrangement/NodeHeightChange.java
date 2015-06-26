package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.UnorderedPair;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by yunyu on 6/25/15.
 */
public class NodeHeightChange extends EdgeParameterChange{
    private Network _network;
    private Map<UnorderedPair, Double> _pairwiseTimeLimit;

    public NodeHeightChange(Map<UnorderedPair, Double> pairwiseTimeLimit){
        _pairwiseTimeLimit = pairwiseTimeLimit;
    }

    public void setInfo(Network network){
        _network = network;
    }

    private Tuple<Double,Double> computeBounds(NetNode targetNode){
        Map<NetNode, Double> node2height = new HashMap<>();
        Map<NetNode, Set<String>> node2taxa = new HashMap<>();
        double upperBound = Double.POSITIVE_INFINITY;
        double lowerBound = -1;
        int parentCount = targetNode.getParentCount();
        if(parentCount==0) parentCount = -1;
        int parentVisited = 0;
        for(Object o: Networks.postTraversal(_network)){
            NetNode node = (NetNode)o;
            double height = 0;
            Set<String> taxa = new HashSet<>();
            if(node.isLeaf()){
                taxa.add(node.getName());
            }
            Set<String> intersection = null;
            List<NetNode> childNodes = null;
            boolean isParentOfTargetNode = false;
            for(Object childO: node.getChildren()){
                NetNode childNode = (NetNode)childO;
                if(childNode == targetNode){
                    parentVisited++;
                    isParentOfTargetNode = true;
                }
                if(targetNode == node){
                    lowerBound = Math.max(lowerBound, node2height.get(childNode));
                }
                if(height==0){
                    height = node2height.get(childNode) + childNode.getParentDistance(node);
                }
                if(node == targetNode && !node.isNetworkNode()) {
                    if(childNodes == null){
                        childNodes = new ArrayList<>();
                    }
                    childNodes.add(childNode);
                    if (intersection == null) {
                        intersection = new HashSet<>();
                        intersection.addAll(node2taxa.get(childNode));
                    } else {
                        intersection.retainAll(node2taxa.get(childNode));
                    }
                }
                taxa.addAll(node2taxa.get(childNode));
            }
            if(isParentOfTargetNode){
                upperBound = Math.min(upperBound, height);
            }
            if(node == targetNode && !node.isNetworkNode()){
                for(int i=0; i<childNodes.size(); i++){
                    Set<String> taxa1 = node2taxa.get(childNodes.get(i));
                    for(int j=i+1; j<childNodes.size(); j++){
                        Set<String> taxa2 = node2taxa.get(childNodes.get(j));
                        for(String taxon1: taxa1){
                            if(intersection.contains(taxon1))
                                continue;
                            for(String taxon2: taxa2){
                                if(intersection.contains(taxon2))
                                    continue;
                                upperBound = Math.min(upperBound, _pairwiseTimeLimit.get(new UnorderedPair(taxon1, taxon2)));
                            }
                        }
                    }
                }
            }
            node2height.put(node, height);
            node2taxa.put(node, taxa);
            if(parentVisited == parentCount){
                break;
            }
        }
        return new Tuple<>(lowerBound,upperBound);
    }

    public boolean performOperation(){
        double newNodeHeight = _targetEdgeBrlen;
        if(newNodeHeight==-1){
            //_targetEdgeBrlen = _targetEdge.Item2.getParentDistance(_targetEdge.Item1);
            Tuple<Double,Double> bounds = computeBounds(_targetEdge.Item2);
            newNodeHeight = drawRandomParameter(_targetEdgeBrlen, bounds.Item1, bounds.Item2);
        }
        _targetEdge.Item2.setParentDistance(_targetEdge.Item1, newNodeHeight);
        return true;
    }


    public void undoOperation(){
        performOperation();
    }
}
