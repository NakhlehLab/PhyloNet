package edu.rice.cs.bioinfo.programs.phylonet.algos.summarize;
/*
 *@ClassName:   majorTree
 *@Description  This class is to summarize the major tree by removing the edge with lower
 *              inheritance probabilities for each reticulation node
 *@Author: Zhen Cao
 *@Date:  2019-07-09 18:08
 *@Version: 1.0
 */

import com.google.common.collect.Lists;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

public class majorTree {
    public static boolean _debug = false;
    public static class VComparator implements Comparator<Map.Entry<Tree, Double>>
    {
        public int compare(Map.Entry<Tree, Double> mp1, Map.Entry<Tree, Double> mp2)
        {
            if (mp2.getValue() - mp1.getValue() > 0.00001){
                return 1;
            }
            else if (mp2.getValue() - mp1.getValue() < -0.00001){
                return -1;

            }
            else{
                return 0;
            }
        }


    }

    public static void removeReti(Network network, List<NetNode> retiList, int retiIndex){
        if (retiList == null || retiList.size() == 0) return;
        if (network.getReticulationCount() == 0) return;
        if (retiIndex >= retiList.size()) return;

        NetNode curReti = retiList.get(retiIndex);
        List<NetNode<NetNodeInfo>> parents = IterableHelp.toList(curReti.getParents());
        if (_debug) {
            System.out.println("retiname"+curReti.getName());
        }
        if(parents.size() != 2) {
            System.err.println("parents size() = " + parents.size() + "\n" + network.toString());
        }
        if( Math.abs(curReti.getParentProbability(parents.get(0)) + curReti.getParentProbability(parents.get(1)) - 1.0) > 0.000001) {
            System.err.println("parents probability != 1.0 " + curReti.getName() + "\n" + curReti.toString());
        }

        NetNode parent = curReti.getParentProbability(parents.get(0))<curReti.getParentProbability(parents.get(1))?parents.get(0):parents.get(1);
        if (parent.getChildCount() == 1){
            removeReti(network, retiList, retiIndex+1);
        }
        parents = IterableHelp.toList(curReti.getParents());
        parent = curReti.getParentProbability(parents.get(0))<curReti.getParentProbability(parents.get(1))?parents.get(0):parents.get(1);
        if (parent.getChildCount() == 2){
            parent.removeChild(curReti);
            Networks.removeBinaryNodes(network);
            removeReti(network, retiList, retiIndex+1);
        }
    }

    public static Map<Tree, Double> summarize(Map<Network, Double> networkMap){
        Map<Tree, Double> treeCount = new HashMap<>();

        for (Network<String> keynet: networkMap.keySet()){
            Network<String> network = keynet.clone();
            Networks.autoLabelNodes(network);
            double count = networkMap.get(keynet);
            List<NetNode> netnodeList = Lists.newArrayList(network.getNetworkNodes());
            if (_debug) {
                System.out.println(network.toString());
            }
            removeReti(network, netnodeList, 0);
            boolean inMap = false;
            for(Tree key: treeCount.keySet()){
                if(Trees.haveSameRootedTopology(key, Trees.readTree(network.toString()))){
                    double num = treeCount.get(key);
                    treeCount.remove(key);
                    treeCount.put(key, num+count);
                    inMap = true;
                    break;
                }
            }

            if(!inMap){
                treeCount.put(Trees.readTree(network.toString()), count);

            }
        }
        if (_debug){
            VComparator vc=new VComparator();
            List<Map.Entry<Tree, Double>> MapList = new ArrayList<>();
            MapList.addAll(treeCount.entrySet());
            Collections.sort(MapList, vc);
        }
        return treeCount;
    }

    public static void test(){
        Network<String> n3 = Networks.readNetwork("((L:1.0)#H1:1.0::0.697176914245199,(((((K:1.0,(P:1.0)#H3:1.0::0.6779314464962011):1.0884181724604034)#H2:5.934653324530301::0.38618142938381267,(C:1.0,#H1:1.0::0.302823085754801):5.9084237144406755):5.9397822589449465,((#H3:1.0::0.3220685535037989)#H4:1.0::0.017142685450424855,(#H4:1.0::0.9828573145495751,O:1.0):0.00808293691880342):2.7714903993551276):5.904173024132873,(#H2:2.012640927894127::0.6138185706161874,F:1.0):5.904657745927051):5.937801402829414);");
        Map<Network, Double> networkMap = new HashMap<>();
        networkMap.put(n3, 1.0);
        Map<Tree, Double>  treeCount = summarize(networkMap);
        for (Tree tree : treeCount.keySet()){
            System.out.println(tree.toString());
        }
    }

    public static void main(String[] args) {
        test();
    }
}
