package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 2019-04-21
 * Time: 18:01
 * To change this template use File | Settings | File Templates.
 */
public class MULTreeUtils {

    static <T> NetNode<T> duplicateSubnetworkTopdown(NetNode<T> root, int[] dupIndex) {
        Map<BniNetNode, NetNode> new2old = new HashMap<>();
        Map<NetNode, BniNetNode> old2new = new HashMap<>();
        Queue<NetNode> queue = new LinkedList<>();
        queue.offer(root);
        while(!queue.isEmpty()) {
            NetNode oldNode = queue.poll();
            BniNetNode newNode = new BniNetNode();
            if(oldNode.isLeaf()) newNode.setName(oldNode.getName() + "_Dup" + dupIndex[0]);
            dupIndex[0]++;
            new2old.put(newNode, oldNode);
            old2new.put(oldNode, newNode);


            for(Object childObj : oldNode.getChildren()) {
                NetNode child = (NetNode) childObj;
                if(!old2new.containsKey(child)) queue.offer(child);
            }
        }

        for(NetNode oldNode : old2new.keySet()) {
            BniNetNode newNode = old2new.get(oldNode);
            for(Object childObj : oldNode.getChildren()) {
                NetNode oldChild = (NetNode) childObj;
                BniNetNode newChild = old2new.get(oldChild);

                newNode.adoptChild(newChild, oldChild.getParentDistance(oldNode));
                newChild.setParentSupport(newNode, oldChild.getParentSupport(oldNode));
                newChild.setParentProbability(newNode, oldChild.getParentProbability(oldNode));
                newChild.setParentProbability(newNode, NetNode.NO_PROBABILITY);
            }
        }

        BniNetwork newnet = new BniNetwork(old2new.get(root));
        Networks.removeBinaryNodes(newnet);
        return old2new.get(root);
    }

    public static <T> List<String> getTaxaUnderNode(NetNode<T> node) {
        Queue<NetNode<T>> q = new LinkedList<>();
        Set<NetNode<T>> visited = new HashSet<>();
        q.offer(node);
        visited.add(node);
        List<String> leaves = new ArrayList<>();
        if(node.isLeaf()) {
            leaves.add(node.getName());
            return leaves;
        }

        while(!q.isEmpty()) {
            NetNode<T> cur = q.poll();
            for(NetNode<T> child : cur.getChildren()) {
                if(visited.contains(child)) continue;
                visited.add(child);
                if(child.isLeaf()) {
                    leaves.add(child.getName());
                } else {
                    q.offer(child);
                }
            }
        }

        return leaves;
    }

    public static Tuple<Network, Map<String, Double>> GetMULTree(Network network) {
        Network multree = network.clone();
        Queue<NetNode<SNAPPData[]>> queue = new LinkedList<>();
        Map<String, Double> probs = new HashMap<>();

        for (Object node : Networks.postTraversal(multree)) {
            ((NetNode<SNAPPData[]>) node).setData(new SNAPPData[1]);
            queue.offer((NetNode<SNAPPData[]>) node);
        }

        int []dupIndex = new int[1];
        dupIndex[0] = 1;
        while(!queue.isEmpty()) {
            NetNode<SNAPPData[]> h = queue.poll();
            if(h.isNetworkNode()) {
                Iterator it = h.getParents().iterator();
                NetNode<SNAPPData[]> u = (NetNode<SNAPPData[]> )it.next();
                NetNode<SNAPPData[]> v = (NetNode<SNAPPData[]> )it.next();
                NetNode<SNAPPData[]> w = (NetNode<SNAPPData[]> ) h.getChildren().iterator().next();

                NetNode<SNAPPData[]> wp = duplicateSubnetworkTopdown(w, dupIndex);
                double distance1 = w.getParentDistance(h) + h.getParentDistance(v);
                double inheritanceProb1 = h.getParentProbability(v);
                double popsize1 = h.getParentSupport(v);
                v.removeChild(h);
                v.adoptChild(wp, distance1);
                wp.setParentSupport(v, popsize1);

                double distance0 = w.getParentDistance(h) + h.getParentDistance(u);
                double inheritanceProb0 = h.getParentProbability(u);
                double popsize0 = h.getParentSupport(u);
                h.removeChild(w);
                u.removeChild(h);
                u.adoptChild(w, distance0);
                w.setParentSupport(u, popsize0);

                List<String> leaves0 = getTaxaUnderNode(w);
                List<String> leaves1 = getTaxaUnderNode(wp);

                for(String leaf1 : leaves1) {
                    if(!probs.containsKey(leaf1)) probs.put(leaf1, 1.0);

                    String oldName = leaf1.substring(0, leaf1.lastIndexOf("_Dup"));

                    probs.put(leaf1, probs.getOrDefault(oldName, 1.0) * inheritanceProb1);
                }

                for(String leaf0 : leaves0) {
                    if(!probs.containsKey(leaf0)) probs.put(leaf0, 1.0);
                    probs.put(leaf0, probs.get(leaf0) * inheritanceProb0);
                }


            }
        }

        return new Tuple<>(multree, probs) ;
    }

    public static void main(String args[]) {
        double alpha1 = 0.5;
        double alpha2 = 0.4;
        double alpha3 = 0.6;
        double alpha4 = 0.8;
        double a = 0;
        double b = 0;
        double c = 0;
        double theta = 0.1;
        double gamma = 0.55;

        String netstring = String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1-alpha2, gamma, alpha1, alpha4-alpha1, alpha3, alpha3-alpha2, 1.0-gamma, alpha4-alpha3);
        netstring = "[0.01]((((B:0.01:0.01)#H1:0.05:0.01:0.5)#H2:0.01:0.01:0.5,A:0.07:0.01):0.03:0.01,((C:0.02:0.01,#H1:0.01:0.01:0.5):0.06:0.01,#H2:0.02:0.01:0.5):0.02:0.01);";

        Network net = Networks.readNetwork(netstring);
        Network multree = GetMULTree(net).Item1;

        System.out.println(multree);
    }
}
