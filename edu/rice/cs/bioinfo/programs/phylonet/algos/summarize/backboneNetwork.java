package edu.rice.cs.bioinfo.programs.phylonet.algos.summarize;

/*
 *@ClassName: backboneNetwork
 *@Description: This class is to summarize the backbone networks for a set of candidate networks
 *              and supports.
 *
 *@Author: Zhen Cao
 *@Date:  2019-07-09 18:08
 *@Version: 1.0
 */

import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.Pipeline;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;


import java.util.*;


public class backboneNetwork {
    private static boolean _debug = false;

    public static class VComparator implements Comparator<Map.Entry<Network, Double>>
    {
        public int compare(Map.Entry<Network, Double> mp1, Map.Entry<Network, Double> mp2)
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
     * @Description: This function is to get all the backbone networks and their support
     * @Param: networkMap   input network and its support
     * @Author: Zhen Cao
     * @Date: 2019-08-15
     */
    public static Map<Network, Double> summarize(Map<Network, Double> networkMap){

        Map<Network, Double> backboneCount = new HashMap<>();
        for (Network<String> network: networkMap.keySet()){
            Networks.removeInternalNodeNames(network);
            List<Network> backboneList2 = Pipeline.getAllBackboneNets(network, Integer.MAX_VALUE);
            backboneList2.add(0, network.clone());
            List<Network> backboneList = new ArrayList<>();
            List<Integer> sameList = new ArrayList<>();
            for(Network b2: backboneList2){
                boolean dup = false;
                for (Network b1: backboneList){
                    if(Networks.hasTheSameTopology(b2, b1)){
                        dup = true;
                        break;
                    }
                }
                if(!dup){
                    backboneList.add(b2);
                }
            }

            double count = networkMap.get(network);
            if(_debug){
                System.out.println(network);
                System.out.println(backboneList.size());

            }
            for(Network bb: backboneList){
                if(_debug){
                    System.out.println(bb.toString());
                }

                boolean inMap = false;
                for(Network key: backboneCount.keySet()){
                    if(Networks.hasTheSameTopology(bb, key)){
                        double num = backboneCount.get(key);
                        backboneCount.put(key, count+num);
                        inMap = true;
                        break;
                    }

                }
                if(!inMap){
                    backboneCount.put(bb, count);
                }
            }
            if (_debug){
                System.out.println("------------------------");

            }

        }
        if(_debug){
            System.out.println("\n\n");
            VComparator vc=new VComparator();
            List<Map.Entry<Network, Double>> backboneMapList = new ArrayList<>();
            backboneMapList.addAll(backboneCount.entrySet());
            Collections.sort(backboneMapList, vc);
            for (Map.Entry<Network, Double> entry: backboneMapList){
                System.out.println(entry.getKey()+"\t"+entry.getValue());
            }
        }

        return backboneCount;
    }
}
