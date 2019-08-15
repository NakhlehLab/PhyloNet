package edu.rice.cs.bioinfo.programs.phylonet.algos.summarize;

/*
 *@ClassName: decomposedTrees
 *@Description: This class is to summarize a map of decomposed trees and their support for a set of candidate networks
 *@Author: Zhen Cao
 *@Date:  2019-07-09 16:32
 *@Version: 1.0
 */

import com.google.common.collect.Lists;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

//checked tree has two same name leaf or network, otherwise rename as x1 x2
//when the reticulation node has different names under different networks, add a label for reticulation node
public class decomposedTrees {
    private static boolean _debug = false;



    /**
     * @Description:
     * @Param: networkMap
     * @Author: Zhen Cao
     * @Date: 2019-08-15
     */
    public static Map<Network, Double> summarize(Map<Network, Double> networkMap){
        Map<Network, Double> subNetworkMap = new HashMap<>();
        for(Network key: networkMap.keySet()){
            Network network = key.clone();
            double count = networkMap.get(key);
            Set<Network>  networkSet = NetworkPartition(network);



            for(Network childnetwork: networkSet){
                boolean inMap = false;
                for(Network subnet: subNetworkMap.keySet()){
                    if(Networks.hasTheSameTopology(childnetwork, subnet)){
                        double num = subNetworkMap.get(subnet);
                        subNetworkMap.put(subnet, num + count);
                        inMap = true;
                        break;
                    }
                }
                if(!inMap){
                    subNetworkMap.put(childnetwork, count);
                }
                if(_debug){
                    System.out.println(childnetwork.toString());

                }
            }
            if (_debug) System.out.println("---------------");
        }
        return subNetworkMap;
    }


    /**
     * @Description: this function is to get a set of k+1 trees decomposed trees for a network with #reti=k.
     * @Param:       the network to get the set of decomposed trees.
     * @Author: Zhen Cao
     * @Date: 2019-07-29
     */
    public static Set<Network> NetworkPartition(Network<Integer> network){
        Set<Network> ChildNetworkSet = new HashSet<>();
        Networks.autoLabelNodes(network);
        if (_debug) System.out.println("original:"+ network.toString());

//        System.out.println(network);
        ChildNetworkSet.add(network);
        for(NetNode<Integer> netnode : network.getNetworkNodes()){
            int i = 1;
            if (_debug) System.out.println(netnode.getName());
            Iterable<NetNode<Integer>> parents = netnode.getParents();
            List<NetNode> parentList = Lists.newArrayList(parents);//note: has to be changed to list

            for(NetNode parent: parentList){

                double dist = netnode.getParentDistance(parent);
                BniNetNode newleaf = new BniNetNode();
                newleaf.setName("#"+netnode.getName());
                parent.removeChild(netnode);
                parent.adoptChild(newleaf, dist);
                i++;
            }
            BniNetwork newNetwork = new BniNetwork((BniNetNode)netnode);
            ChildNetworkSet.add(newNetwork);
        }
        renameNetwork(ChildNetworkSet);
        return ChildNetworkSet;
    }


    /**
     * @Description:
     * @Param:
     * @Author: Zhen Cao
     * @Date: 2019-08-15
     */
    public static String removeRepeatChar(String s) {
        if (s == null) {
            return "";
        }

        String underline = "";
        for(int i = 0; i < s.length(); i++){
            if(s.charAt(i)=='_'){
                underline += "_";
            }
        }
        String[] chs = s.replace("_","").split("");

        List stringlist = Lists.newArrayList(chs);
        Collections.sort(stringlist);
        s = String.join("", stringlist);
        StringBuffer sb = new StringBuffer();
        int i = 0;
        int len = s.length();
        while (i < len) {
            char c = s.charAt(i);
            sb.append(c);
            i++;
            while (i < len && s.charAt(i) == c) {
                i++;
            }
        }
        return sb.toString()+underline;
    }


    /**
     * @Description:        This function is to rename the leaves of decomposed trees
     * @Param: networkSet   the set of decomposed trees to rename
     * @Author: Zhen Cao
     * @Date: 2019-08-15
     */

    public static void renameNetwork(Set<Network> networkSet){
        Map<NetNode, List<NetNode>> leafrootmap = new HashMap<>();
        for (Network retinet: networkSet){
            Iterable<NetNode<Integer>> leaves = retinet.getLeaves();
            List<NetNode> leafList = Lists.newArrayList(leaves);
            leafrootmap.put(retinet.getRoot(), leafList);
        }
        boolean isHybrid = true;

        Map<String, String> nameMap = new HashMap<>();
        List<NetNode> keyList = new ArrayList<>(leafrootmap.keySet());
        int i = 0;
        while(!keyList.isEmpty()){
            i %= keyList.size();
            NetNode node = keyList.get(i);
            isHybrid = false;
            String name = "";
            for(NetNode leaf: leafrootmap.get(node)){
                String leafName = leaf.getName();
                if(leafName.contains("#")){

                    String newName = nameMap.get(leafName.substring(1));
                    if (newName!=null){
                        leaf.setName(newName);
                        name += newName;
                    }
                    else {
                        isHybrid = true;
                    }
                }
                else{
                    name += leafName;
                }
            }
            if(isHybrid){
                i ++;
            }
            else{

                name = removeRepeatChar(name);
                name =  name + "_";
                nameMap.put(node.getName(), name);
                node.setName(name);
                keyList.remove(i);
            }
        }


    }

    public static void Test(){
//        String n1 = "((A,(((F,(E,Y#H2:::0.1)))X#H1:::0.2,(B,(C,(D)Y#H2:::0.9)))),((X#H1:::0.8,G),H));";
        String n1 = "((A,(((F,(E,Y#H2)))X#H1,(B,(C,(D)Y#H2)))),((X#H1,G),H));";
        String n2 = "((A,(((F,(E,Y#H2:::0.2)))X#H1:::0.1,(B,(C,(D)Y#H2:::0.8)))),((X#H1:::0.9,H),G));";
        String n3 = "((B,(((F,(E,Y#H2:::0.8)))X#H1:::0.9,(A,(C,(D)Y#H2:::0.2)))),((X#H1:::0.1,H),G));";
        String n4 = "((A,((B)Y#H1)X#H2),(X#H2,(Y#H1,C)));";
        String n5 = "((A,((B)Y#H1)X#H2),(X#H2,(Y#H1,C)));";
//        String n6 = "(((((((E:6.191736422399999,(F:4.299816959999999)#H1:1.8919194623999998::0.7500000000000001)S13:2.7243640258559987,(N:1.2,M:1.2)S14:7.716100448255998)S12:1.7832200896511985,(D:7.430083706879999)#H2:3.2692368310271975::0.8000000000000002)S4:35.30579937146248,(((((((((J:2.0736,I:2.0736)S15:3.0861803519999995,#H1:0.859963392::0.2499999999999999)S11:7.679404293488635,C:12.839184645488634)S10:2.5678369290977265,B:15.407021574586361)S9:3.08140431491727,A:18.48842588950363)S8:3.697685177900727,#H2:14.75602736052436::0.19999999999999984)S7:4.437222213480872,(((H:2.48832,G:2.48832)S18:0.4976639999999999,(((O:1.0,P:1.0)S20:0.43999999999999995,L:1.44)S19:0.28800000000000003,K:1.728)S17:1.2579839999999998)S16:0.5971967999999999)#H4:23.04015248088523::0.7000000000000001)S6:5.324666656177044,#H4:28.364819137062273::0.29999999999999993)S5:6.389599987412456)#H3:7.667519984894945::0.65)S3:20.242252760122646,(#H3:16.868543966768875::0.35)#H5:11.041228778248716::0.6)S2:13.249474533898464,#H5:24.29070331214718::0.4)S1:20.503152796609214,Z:100.0);";
        Map<Network, Double> networkMap = new HashMap<>();
//        Network n = Networks.readNetwork(n1);
        networkMap.put(Networks.readNetwork(n1),1.0);
        networkMap.put(Networks.readNetwork(n2),2.0);
        networkMap.put(Networks.readNetwork(n3),1.0);
        networkMap.put(Networks.readNetwork(n4),1.0);
//        networkMap.put(Networks.readNetwork(n6),1.0);
        Map<Network, Double> subNetworkMap = summarize(networkMap);
        System.out.println("size="+subNetworkMap.size());
        for(Network sub:subNetworkMap.keySet()){
            System.out.println(sub.toString());
            System.out.println(subNetworkMap.get(sub));
            System.out.println("----------");
            System.out.println();
        }
//        String x = "(Y)X;";
//        Network<String> n = Networks.readNetwork(x);
//        for(NetNode<String> node: n.getLeaves()){
//            System.out.println(node.getName());
//        }
    }

    public static void main(String[] args) {
        Test();

    }

}
