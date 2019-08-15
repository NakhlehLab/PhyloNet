package edu.rice.cs.bioinfo.programs.phylonet.algos.summarize;
/*
 *@ClassName: TestMethods
 *@Description
 *@Author: Zhen Cao
 *@Date:  2019-07-09 19:17
 *@Version: 1.0
 */

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkInfo;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.commands.NetworkTransformer;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;



public class TestMethods {


    public static void Test1(int testmethod){
        String n1 = "((A,(((F,(E,Y#H2:::0.1)))X#H1:::0.2,(B,(C,(D)Y#H2:::0.9)))),((X#H1:::0.8,G),H));";
        String n2 = "((A,(((F,(E,Y#H2:::0.2)))X#H1:::0.1,(B,(C,(D)Y#H2:::0.8)))),((X#H1:::0.9,H),G));";
        String n3 = "((B,(((F,(E,Y#H2:::0.8)))X#H1:::0.9,(A,(C,(D)Y#H2:::0.2)))),((X#H1:::0.1,H),G));";

//        List<Network> networkList = new ArrayList<>();
//        networkList.add(Networks.readNetwork(n1));
//        networkList.add(Networks.readNetwork(n2));
//        networkList.add(Networks.readNetwork(n3));
//        Network<NetNodeInfo> n = Networks.readNetwork(n1);
//        Networks.autoLabelNodes(n);


        Map<Network, Double> networkMap = new HashMap<>();
        networkMap.put(Networks.readNetwork(n1), 1.0);
        networkMap.put(Networks.readNetwork(n2), 2.0);
        networkMap.put(Networks.readNetwork(n3), 1.0);

//        String file = "/Users/zhen/Desktop/Zhen/research/phylogenetics/summarize/nex/wheat.nex";
//        String file = "/Users/zhen/Desktop/Zhen/research/phylogenetics/summarize/nex/mouse.nex";
//        Map<Network, Double> networkMap = getInputNetworkMap(file);


//        int testmethod = 0;
        switch (testmethod){
            case 0:
                Map<Tree, Double> maxtreemap = maxTree.summarize(networkMap);
                for(Tree t: maxtreemap.keySet()){
                    System.out.println(t.toString()+"\tcount:"+maxtreemap.get(t).toString());
                }
                break;
            case 1:
                Map<Network, Double> backboneMap = backboneNetwork.summarize(networkMap);
                for(Network backnet: backboneMap.keySet()){
                    System.out.println(backnet.toString()+"\t\t"+backboneMap.get(backnet).toString());
                }
                break;
            case 2:
                Map<Tree, Double> majortreemap = majorTree.summarize(networkMap);
                for(Tree t: majortreemap.keySet()){
                    System.out.println(t.toString()+"\t\tcount:"+majortreemap.get(t).toString());
                }
                break;
            case 3:
                nonEquivalent.summarize(networkMap);
//                nonEquivalent.labelNodes(n);
                break;
            case 4:
                Map<Network, Double> decomposedtreemap = decomposedTrees.summarize(networkMap);
                for(Network t: decomposedtreemap.keySet()){
                    System.out.println(t.toString()+"\t\tcount:"+decomposedtreemap.get(t).toString());
                }
                break;
            case 5:
                Map<netNodeTuple, Double> tripart = tripartition2.summarize(networkMap);
                break;
            default:
                break;

        }

//        for (NetNode node: n.dfs()){
//            System.out.println(node.getName());
//        }
//        System.out.println("-----------");
//        for (Object o: Networks.postTraversal(n)){
//            NetNode netnode = (NetNode) o;
//            System.out.println(netnode.getName());
//        }

//        NetworkTransformer.toENewick(new NetworkNonEmpty(n1));


    }
    public static void Test2(int testmethod){
//        String file = "/Users/zhen/Desktop/Zhen/research/phylogenetics/summarize/nex/wheat.nex";
        String file = "/Users/zhen/Desktop/Zhen/research/phylogenetics/bookChapter/data/Reti3_C/out4/MCMC_GT_pl8_2_true.nex.out";
//        String file = "/Users/zhen/Desktop/Zhen/research/phylogenetics/summarize/nex/mouse.nex";
        Map<Network, Double> networkMap = getInputNetworkMap(file);


        switch (testmethod){
            case 0:
                Map<Tree, Double> maxtreemap = maxTree.summarize(networkMap);
                for(Tree t: maxtreemap.keySet()){
                    System.out.println(t.toString()+"\tcount:"+maxtreemap.get(t).toString());
                }
                break;
            case 1:
                Map<Network, Double> backboneMap = backboneNetwork.summarize(networkMap);
                for(Network backnet: backboneMap.keySet()){
                    System.out.println(backnet.toString()+"\t\t"+backboneMap.get(backnet).toString());
                }
                break;
            case 2:
                Map<Tree, Double> majortreemap = majorTree.summarize(networkMap);
                for(Tree t: majortreemap.keySet()){
                    System.out.println(t.toString()+"\t\tcount:"+majortreemap.get(t).toString());
                }
                break;
            case 3:
                nonEquivalent.summarize(networkMap);
//                nonEquivalent.labelNodes(n);
                break;
            case 4:
                Map<Network, Double> decomposedtreemap = decomposedTrees.summarize(networkMap);
                for(Network t: decomposedtreemap.keySet()){
                    System.out.println(t.toString()+"\t\tcount:"+decomposedtreemap.get(t).toString());
                }
                break;
            case 5:
                Map<netNodeTuple, Double> tripart = tripartition2.summarize(networkMap);
                for(netNodeTuple t: tripart.keySet()){
                    System.out.println(t.toString()+"\t\tcount:"+tripart.get(t).toString());
                }
                break;
            default:
                break;

        }
    }

    public static Map<Network, Double> getInputNetworkMap(String file){
        Map<Network, Double> networkMap = new HashMap<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            String s;
            String[] ss;
            while((s = br.readLine()) != null){
                s = s.trim();
                ss = s.split("\\s");
                if(ss.length == 2){
                    networkMap.put(Networks.readNetwork(ss[0]), Double.parseDouble(ss[1]));
                }
            }

        }catch (Exception e){
            e.printStackTrace();
        }
        return networkMap;
    }

    public static void main(String[] args) {
//        Test1(5);
        Test2(5);

    }


}
