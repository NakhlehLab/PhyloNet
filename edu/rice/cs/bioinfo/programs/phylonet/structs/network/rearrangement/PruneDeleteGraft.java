package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by yunyu on 11/7/14.
 */
public class PruneDeleteGraft extends NetworkRearrangementOperation{

    public void setParameters(Network network, Tuple<NetNode,NetNode> targetNode){
        _network = network;
        _targetEdge = targetNode;
    }

    public void undoOperation(){
    }

    public boolean performOperation(){
        List<NetNode> parents = new ArrayList<>();
        for(Object parent: _targetEdge.Item1.getParents()){
            parents.add((NetNode)parent);
        }
        for(NetNode parent: parents){
            parent.removeChild(_targetEdge.Item1);
        }
        for(Object child: _targetEdge.Item1.getChildren()){
            if(child!=_targetEdge.Item2) {
                _targetEdge.Item1.removeChild((NetNode)child);
            }
        }
        Network tempNet = new BniNetwork((BniNetNode)_targetEdge.Item1);
        //System.out.println(tempNet.toString());
        //Set<String> leafSet = new HashSet<>();
        for(Object leaf: tempNet.getLeaves()){
            //System.out.println("Removing leaf " + ((NetNode)leaf).getName());
            deleteNode(_network, _network.findNode(((NetNode) leaf).getName()));
            //System.out.println(_network);
        }
        graftToRandomEdge(_network, _targetEdge.Item1);

        return true;
    }


    private void adjustNetwork(Network network){

        //System.out.print(network);


        Set<NetNode> leaves = new HashSet<>();
        for(Object nodeO: network.getLeaves()){
            NetNode node = (NetNode)nodeO;
            if(!node.isNetworkNode()){
                leaves.add(node);
            }
        }

        boolean update;
        do{
            update = false;
            for(Object nodeO: Networks.postTraversal(network)){
                NetNode node = (NetNode)nodeO;
                if(node.isLeaf() && !leaves.contains(node)){
                    node.removeItself();
                    update = true;
                    break;
                }
            }
            //System.out.println("Before: " + network);
            Networks.removeBinaryNodes(network);
            //System.out.println("Here:" + network);
        }while(update);

        NetNode root = network.getRoot();
        while(root.getChildCount()==1){
            NetNode child = (NetNode)root.getChildren().iterator().next();
            if(!child.isLeaf()) {
                root.removeChild(child);
                root = child;
            }
            else{
                break;
            }
        }

        network.resetRoot(root);

        for(Object nodeO: network.getNetworkNodes()){
            NetNode node = (NetNode)nodeO;
            if(node.isTreeNode() && !node.isRoot()){
                NetNode parentNode = (NetNode)(node.getParents().iterator().next());
                node.setParentProbability(parentNode, NetNode.NO_PROBABILITY);
            }
        }

        //System.out.println("done");

    }

    private void graftToRandomEdge(Network network, NetNode toGraft){
        List<Tuple<NetNode,NetNode>> allEdges = new ArrayList<>();
        for(Object nodeO: Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeO;
            for(Object childO: node.getChildren()){
                NetNode childNode = (NetNode)childO;
                Tuple<NetNode,NetNode> edge = new Tuple<>(node, childNode);
                allEdges.add(edge);
            }
        }
        if(allEdges.size()==1){
            Tuple<NetNode,NetNode> onlyEdge = allEdges.get(0);
            toGraft.adoptChild(onlyEdge.Item2, onlyEdge.Item2.getParentDistance(onlyEdge.Item1));
            onlyEdge.Item1.removeChild(onlyEdge.Item2);
            network.resetRoot(toGraft);
        }
        else {
            Random random = new Random();
            //System.out.println("After delete: " + network.toString());
            Tuple<NetNode, NetNode> sourceEdge = allEdges.get(random.nextInt(allEdges.size()));
            //System.out.println("Regraft to edge: (" + sourceEdge.Item1.getName() + "," + sourceEdge.Item2.getName() + ")\n");

            double[] brlens = new double[2];
            double[] inheriProbs = new double[2];
            randomlyPartitionAnEdge(sourceEdge, brlens, inheriProbs);
            addNodeToAnEdge(toGraft, sourceEdge, brlens, inheriProbs);
        }
    }

/*
    private void deleteNode(Network network, NetNode node){
        System.out.println(node.getName());
        List<Object> parents = new ArrayList<>();
        for(Object parent: node.getParents()){
            parents.add(parent);
        }
        for(Object parent: parents){
            NetNode parentNode = (NetNode)parent;
            if(parentNode.isNetworkNode()){
                deleteNode(network, parentNode);
            }
            else{
                Tuple<NetNode,NetNode> source = findParentAndAnotherChild(parentNode, node);
                if(source.Item1!=null){
                    double[] temp1 = new double[2];
                    double[] temp2 = new double[2];
                    removeNodeFromAnEdge(parentNode,source,temp1,temp2);
                    if(source.Item2.isNetworkNode()){
                        Tuple<NetNode,NetNode> further = findAnotherParentAndChild(source.Item2, source.Item1);
                        removeNodeFromAnEdge(source.Item2,further,temp1,temp2);
                        source.Item1.removeChild(source.Item2);
                    }
                }
                else{
                    parentNode.removeChild(node);
                }

            }
        }
        System.out.println(network);
    }
*/
private void deleteNode(Network network, NetNode node){
    do{
        NetNode parentNode = (NetNode)(node.getParents().iterator().next());
        if(parentNode.isNetworkNode()){
            parentNode.removeChild(node);
            break;
        }
        else{
            Tuple<NetNode,NetNode> temp = findParentAndAnotherChild(parentNode, node);
            if(temp.Item2.isNetworkNode()){
                parentNode.removeChild(temp.Item2);
                parentNode.removeChild(node);
            }
            else{
                parentNode.removeChild(node);
                break;
            }
        }
        node = parentNode;
    }while(true);

    adjustNetwork(network);
}




    public static void main(String[] args){
        Network network = Networks.readNetwork(" ((((e:1.0,d:1.0)I5:1.1236846065108836,f:1.0)I4:1.486119306006081)I2#H1:1.0981674846410259::0.26353934643510035,(a:3.748252512132769,((b:1.1252926387906674,c:2.61086859954818)I1:0.5945071891725606,I2#H1:0.036486372127869995::0.7364606535648996)I8:0.9079903913505394)I3:1.203317552596218)I0;");
        PruneDeleteGraft pdg = new PruneDeleteGraft();
        NetNode node = network.findNode("d");
        pdg.deleteNode(network,node);
        node = network.findNode("f");
        pdg.deleteNode(network,node);
        node = network.findNode("e");
        pdg.deleteNode(network,node);
        System.out.println(network);
    }

}
