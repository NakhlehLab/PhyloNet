package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.simulation;
/*
 * @ClassName:   OutgroupAdder
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        12/14/23 5:31 PM
 */
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;

import java.util.*;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.Pipeline.initNetHeights;

public class OutgroupAdder {
    private String original_net = null;
    private Network _network = null;
    /* Constructor */
    public OutgroupAdder(String net) {
        original_net = net;
    }

    public boolean isUltrametric() {
//        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network)) {
//            double height = node.getData().getHeight();
//            if(node.isLeaf()) {
//                if(height != Utils.DEFAULT_NET_LEAF_HEIGHT) {
//                    System.err.println(height + " vs " + Utils.DEFAULT_NET_LEAF_HEIGHT);
//                    return false;
//                }
//            } else {
//                for(NetNode<NetNodeInfo> child : node.getChildren()) {
//                    double temp = child.getData().getHeight() + child.getParentDistance(node);
//                    if(Math.abs(temp - height) > 0.000001) {
//                        System.err.println(node.getName() + " - " + height + " vs " + temp);
//                        return false;
//                    }
//                }
//            }
//        }
        return true;
    }

    // check the length of each tip to the root, if the length is not the same, change the length to the same
    public Network makeUltrametric(Network network){
        // for each leaf, get all path from the root to it, and get the length of the path
        Map<NetNode, Double> nodeDist2Root = new HashMap<>();
//        for(Object o : Networks.postTraversal(network)) {
//            NetNode node = (NetNode) o;
//            nodeDist2Root.put(node, new ArrayList<>());
//
//        }
        makeUltrametricHelper(network.getRoot(), 0, nodeDist2Root);
        double maxTipDistance = 0;
        for(Object o : Networks.postTraversal(network)) {
            NetNode node = (NetNode) o;
//            if (node.isNetworkNode()) {
//
//            }
            if (node.isLeaf()){
                maxTipDistance = Math.max(maxTipDistance, nodeDist2Root.get(node));
            }

        }

        System.out.println(maxTipDistance);
        for(Object o : Networks.postTraversal(network)) {
            NetNode node = (NetNode) o;
//            if (node.isNetworkNode()) {
//
//            }

            if (node.isLeaf()){
                if (maxTipDistance != nodeDist2Root.get(node)) {
                    NetNode parent = (NetNode) node.getParents().iterator().next();
                    System.out.println(nodeDist2Root.get(node));
                    node.setParentDistance(parent, maxTipDistance - nodeDist2Root.get(parent));
                }

            }
        }
        return network;
    }

    // check the length of each node to the root, if the length is not the same, change the branch length to be consistent with the max length
    private void makeUltrametricHelper(NetNode node, double currentDistance, Map<NetNode, Double> nodeDist2Root) {
//        double maxDistance = 0;
        if (node.isLeaf()){
            return;
        }
        for (Object o : node.getChildren()) {
            NetNode child = (NetNode) o;
            if (child.isNetworkNode()){
                double childDistance = currentDistance + child.getParentDistance(node);

                if (nodeDist2Root.containsKey(child) && (nodeDist2Root.get(child) > childDistance)) {
                    child.setParentDistance(node, nodeDist2Root.get(child) - currentDistance);
                }
                else if(nodeDist2Root.containsKey(child) && (childDistance > nodeDist2Root.get(child))){
                    for (Object op : child.getParents()) {
                        NetNode otherParent = (NetNode) op;
                        if (otherParent != node){
                            child.setParentDistance(otherParent, childDistance - nodeDist2Root.get(otherParent));
                            break;
                        }
                    }
                    nodeDist2Root.put(child, childDistance);
                }
                else {
                    nodeDist2Root.put(child, childDistance);
                }
            }
            else{
                nodeDist2Root.put(child, currentDistance + child.getParentDistance(node));
            }

            makeUltrametricHelper(child, nodeDist2Root.get(child), nodeDist2Root);
//            maxDistance = Math.max(maxDistance, childDistance);
        }
//        return maxDistance;
    }

    public String addOutgroup(String outgroup, double branchLength) {
        Network<NetNodeInfo> outgroup_net = Networks.readNetwork(original_net);
        Networks.removeBinaryNodes(outgroup_net);
        initNetHeights(outgroup_net);

        NetNode<NetNodeInfo> outgroup_node = new BniNetNode();
        outgroup_node.setName(outgroup);
        NetNode<NetNodeInfo> root = new BniNetNode();

        root.adoptChild(outgroup_net.getRoot(), branchLength - outgroup_net.getRoot().getData().getHeight());
        root.adoptChild(outgroup_node, branchLength);
        outgroup_net.resetRoot(root);


        System.out.println(outgroup_net.toString());
        return outgroup_net.toString();
    }

    public static void toyTest(){
//        String net = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);";
        String net = "(((((t31:5.921989392,((t25:0.08603025499,t18:0.08603025499):2.966215911,(t14:0.3723561354,(t36:0.2720867863,t24:0.2720867863):0.1002693491):2.67989003):2.869743226):0.1357092011,((t8:0.6792094518,(t29:0.4836171251)#H43:0.1955923267::0.5341611349):4.253411065,(((((t35:1.971968257,#H43:1.488351132::0.4658388651):0.7573741672,t28:2.729342424):0.2081455566,(((t38:0.2401680339,t32:0.2401680339):1.345023173,#H41:0.5722789213::0.4464888709):0.8927140673,(t10:0.1843753954,t4:0.1843753954):2.293529879):0.4595827061):0.5598953526,((t39:1.012912286)#H41:0.07691345254::0.5535111291,t6:1.089825738):2.407557595):0.7011607985,t22:4.198544132):0.734076385):1.125078076):0.3023215241,t3:6.360020117):0.8773627279,(t23:2.693809225,t13:2.693809225):4.54357362):3.607386478);";
        OutgroupAdder outgroupAdder = new OutgroupAdder(net);
        String network = outgroupAdder.addOutgroup("E", 100);
        Network n = Networks.readNetwork(network);
        Networks.autoLabelNodes(n);
        Network ultrametricNet = outgroupAdder.makeUltrametric(n);
        System.out.println(ultrametricNet.toString());
    }

    // read the strings in the given path
    public static String[] read_model_phylogeny_add_outgroup(String path){
        String netnewicks = null;
        try {
            netnewicks = new String(java.nio.file.Files.readAllBytes(java.nio.file.Paths.get(path)));
        } catch (Exception e) {
            e.printStackTrace();
        }
        return netnewicks.split("\n");
    }

    //write the strings in the given path
    public static void write_model_phylogeny_add_outgroup(String path, String[] netnewicklist){
        try {
            java.nio.file.Files.write(java.nio.file.Paths.get(path), Arrays.asList(netnewicklist));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void processData(){
        String dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/merge/data/multinet/";
//        String[] tips = {"20", "40", "80", "160"};
        String[] tips = {"80"};
        for (String net_tip : tips){
            String path = dir_path+net_tip+"/network_newicks1500.txt";

            String data_level_path = dir_path+net_tip+"/network_newicks_outgroup1500.txt";
            String[] netnewicklist = read_model_phylogeny_add_outgroup(path);
            String[] outgroupnewick = new String[netnewicklist.length];
            int i = 0;
            for (String newick : netnewicklist){
                OutgroupAdder outgroupAdder = new OutgroupAdder(newick);
                System.out.println(i);
                outgroupnewick[i] = outgroupAdder.addOutgroup("Z", 100);
                i += 1;
            }
            write_model_phylogeny_add_outgroup(data_level_path, outgroupnewick);
            System.out.println("Done with "+net_tip);
        }
    }



    public static void main(String[] args) {
        processData();
//        toyTest();
    }
}
