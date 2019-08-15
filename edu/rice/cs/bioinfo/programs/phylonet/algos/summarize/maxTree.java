package edu.rice.cs.bioinfo.programs.phylonet.algos.summarize;
/*
 *@ClassName: maxTree
 *@Description: This class is to summarize the trees displayed by a set of networks
 *@Author: Zhen Cao
 *@Date:  2019-07-09 17:48
 *@Version: 1.0
 */
import com.google.common.collect.Lists;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;
//tested
public class maxTree {
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
    
    /**
     * @Description:        This function is to get the displayed trees and their support
     * @Param: networkMap   Networks and their support
     * @Author: Zhen Cao
     * @Date: 2019-08-15 
     */
    public static Map<Tree, Double> summarize(Map<Network, Double> networkMap){
        Map<Tree, Double> treeCount = new HashMap<>();
        for (Network<String> network: networkMap.keySet()){
            Networks.autoLabelNodes(network);
            double count = networkMap.get(network);
            if (_debug){
                System.out.println(network);
                System.out.println(count);
                System.out.println(IterableHelp.toList(Networks.getTrees(network)).size());

            }
            List<NetworkTree> ntlist = Lists.newArrayList(Networks.getTrees(network));

            List<Tree> treeList = new ArrayList<>();
            for (NetworkTree nettree: ntlist){
                Tree t1 = Trees.readTree(nettree.makeTree().toString());
                boolean dup = false;
                for (Tree t2 : treeList){
                    if (Trees.haveSameRootedTopology(t1, t2)){
                        dup = true;
                        break;
                    }
                }
                if (!dup){
                    treeList.add(t1);
                }
            }
            for (Tree nt: treeList) {
                boolean inMap = false;
                if (_debug){
                    System.out.println(nt.toString());
                }
                for(Tree key: treeCount.keySet()){
                    if (Trees.haveSameRootedTopology(key, nt)){
                        double num = treeCount.get(key);
                        treeCount.remove(key);
                        treeCount.put(key, num+count);
                        inMap = true;
                        break;

                    }
                }
                if(!inMap){
                    treeCount.put(Trees.readTree(nt.toString()), count);
                }
            }
            if (_debug){
                System.out.println("-----------------------");
            }
        }
        if (_debug) {
            VComparator vc=new VComparator();
            List<Map.Entry<Tree, Double>> MapList = new ArrayList<>();
            MapList.addAll(treeCount.entrySet());
            Collections.sort(MapList, vc);
        }

        return treeCount;
    }

}
