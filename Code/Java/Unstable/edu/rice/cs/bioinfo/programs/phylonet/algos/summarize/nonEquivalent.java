package edu.rice.cs.bioinfo.programs.phylonet.algos.summarize;
/*
 *@ClassName: nonEquivalent
 *@Description: This class implements the labeling method in
 *              "Bayesian Inference of Species Networks from Multilocus Sequence Data"
 *              https://academic.oup.com/mbe/article/35/2/504/4705834
 *@Author: Zhen Cao
 *@Date:  2019-07-09 19:57
 *@Version: 1.0
 */

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.HashSet;
import java.util.Map;

//todo find method for summarize common structure after labeling

public class nonEquivalent {


    public static void summarize(Map<Network, Double> networkMap){

        double size = 0.0;
        for (Network<NetNodeInfo> network: networkMap.keySet()){
            System.out.println(network.toString());
            Networks.autoLabelNodes(network);
            labelNodes(network, networkMap.get(network));
            size += networkMap.get(network);

            System.out.println("--------------------------------------------");
        }

        for (Network<NetNodeInfo> network: networkMap.keySet()){
            computeSupport(network, size);
        }


    }

    public static void labelNodes(Network<NetNodeInfo> network, double num){
        if(Utils._leafLabel.isEmpty()){
            for(NetNode node : network.getLeaves()){
                Utils._leafLabel.put(node.getName(), Utils._counter);
                Utils._counter++;
            }
        }

        labelNode(network.getRoot(), num);


    }

    public static int labelNode(NetNode<NetNodeInfo> node, double num){
        int label;
        if(node.getData() != null){
            return node.getData().getLabel();
        }
        else if(node.isLeaf()){
            label = Utils._leafLabel.get(node.getName());
            node.setData(new NetNodeInfo(label));

        }
        else{
            HashSet<Integer> childSet = new HashSet<>();
            for (NetNode child: node.getChildren()){
                childSet.add(labelNode(child, num));
            }
            if (Utils._setLabel.containsKey(childSet)){
                label = Utils._setLabel.get(childSet);
                double count = Utils._labelCounter.get(label);
                Utils._labelCounter.remove(label);
                Utils._labelCounter.put(label, count + num);
                node.setData(new NetNodeInfo(label));

            }
            else{
                label = Utils._counter;
                Utils._setLabel.put(childSet, label);
                Utils._labelCounter.put(label, num);
                node.setData(new NetNodeInfo(label));
                Utils._counter++;

            }
        }

        return label;
    }

    public static void computeSupport(Network<NetNodeInfo> network, double size){
        for(NetNode<NetNodeInfo>  node: network.dfs()){
            try{
                if(!node.isLeaf()){
                    node.getData().setSupport(Utils._labelCounter.get(node.getData().getLabel())*1.0/size);
                }
            }catch (Exception e){
                e.printStackTrace();
            }

        }
    }


    public static void main(String[] args) {

    }


}
