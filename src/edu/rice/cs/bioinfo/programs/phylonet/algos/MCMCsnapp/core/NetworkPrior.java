package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Given a species network, return log prior(net)
 *
 * Created by wendingqiao on 10/17/14.
 */
public class NetworkPrior {

    private Network _net;

    private double k; // k = 1.0 -> prior on diameters; k = 0.0 -> no prior on diameters.

    private List<NetNode> startLeaves;

    public NetworkPrior(Network net, double k, List<NetNode> startLeaves) {
        this._net = net;
        this.k = k;
        this.startLeaves = startLeaves;
    }

    /**
     * Calculate prior given network.
     * @return    Log of prior
     */
    public double calculate() {
        if(!isValid(_net)) return Double.MIN_VALUE;
        return priorOnReticulations() + priorOnBranchLengths() + priorOnDiameters();
    }

    /*
     * prior on the number of reticulations
     */
    private double priorOnReticulations() {
        double lambda = k; // parameter of poisson distribution
        double priorReti = 0; // e^(-lambda) * (lambda^k / k!);  normalize e^(-lambda)
        int numReti = _net.getReticulationCount();
        int numEdges = _net.getLeafCount() * 2 - 2;
        for(int i = 1; i <= numReti; i++) {
            priorReti += Math.log(lambda / i) - Math.log((double)(numEdges + i*3) * (numEdges + i*3 - 1));
        }
        return priorReti;
    }

    /*
     * prior on branch lengths
     */
    private double priorOnBranchLengths() {
        double priorBranches = 0.0;
        for(Object node : _net.dfs()) {
            NetNode n = (NetNode) node;
            List<NetNode> pars = IterableHelp.toList(n.getParents());
            for(NetNode par : pars) {
                priorBranches -= n.getParentDistance(par);
            }
        }
        return priorBranches;
    }

    /**
     * prior on diameters
     */
    private double priorOnDiameters() {
        Network net = Networks.readNetwork(_net.toString());
        Map<NetNode, Double> distMap = new HashMap<NetNode, Double>();
        Map<NetNode, Tuple<NetNode, NetNode>> parentsMap = new HashMap<>();
        for(Object o : net.getNetworkNodes()) {
            NetNode node = (NetNode) o;
            List<NetNode> parents = IterableHelp.toList(node.getParents());
            if(parents.size() != 2) throw new IllegalArgumentException("Invalid network " + net.toString());
            distMap.put(node, 0.0);
            parentsMap.put(node, new Tuple<NetNode, NetNode>(parents.get(0), parents.get(1)));
            for(NetNode par : parents) {
                distMap.put(node, distMap.get(node) + node.getParentDistance(par));
                if(node.getParentProbability(par) < 0.50) {
                    par.removeChild(node);
                }
            }
        }
        // exponential distribution where lambda = 1.0;
        // ln(lambda * e^(-lambda * x) ) = ln(lambda) - lambda * x = -x;
        double priorDiameter = 0.0;
        for(NetNode key: distMap.keySet()) {
            distMap.put(key, distMap.get(key) + getReticulationNodeDiameter(net, parentsMap.get(key)));
            priorDiameter -= distMap.get(key);
        }
        return priorDiameter;
    }

    public static Map<NetNode, Double> getDiameterMap(String netStr) {
        Network net = Networks.readNetwork(netStr);
        Map<NetNode, Double> distMap = new HashMap<NetNode, Double>();
        Map<NetNode, Tuple<NetNode, NetNode>> parentsMap = new HashMap<>();
        for(Object o : net.getNetworkNodes()) {
            NetNode node = (NetNode) o;
            List<NetNode> parents = IterableHelp.toList(node.getParents());
            if(parents.size() != 2) throw new IllegalArgumentException("Invalid network " + net.toString());
            distMap.put(node, 0.0);
            parentsMap.put(node, new Tuple<>(parents.get(0), parents.get(1)));
            for(NetNode par : parents) {
                distMap.put(node, distMap.get(node) + node.getParentDistance(par));
                if(node.getParentProbability(par) < 0.50) {
                    par.removeChild(node);
                }
            }
        }
        for(NetNode key: distMap.keySet()) {
            distMap.put(key, distMap.get(key) + getReticulationNodeDiameter(net, parentsMap.get(key)));
        }
        return distMap;
    }

    public static double getReticulationNodeDiameter(Network net, Tuple<NetNode, NetNode> parents) {
        List<NetNode> list1 = getParentsList(net, parents.Item1);
        List<NetNode> list2 = getParentsList(net, parents.Item2);
        NetNode mrca = net.getRoot();
        for(NetNode n1: list1) {
            for(NetNode n2 : list2) {
                if(n1.equals(n2)) {
                    mrca = n1;
                    break;
                }
            }
        }
        return getDistance(parents.Item1, mrca) + getDistance(parents.Item2, mrca);
    }

    public static double getDistance(NetNode node, NetNode mrca)  {
        double dist = 0.0;
        NetNode tmp = node;
        while(!tmp.equals(mrca)) {
            NetNode par = (NetNode) tmp.getParents().iterator().next();
            dist += tmp.getParentDistance(par);
            tmp = par;
        }
        return dist;
    }

    public static List<NetNode> getParentsList(Network net, NetNode node) {
        List<NetNode> list = new ArrayList<>();
        NetNode tmp = node;
        while(tmp != net.getRoot()) {
            list.add(tmp);
            tmp = (NetNode) tmp.getParents().iterator().next();
        }
        return list;
    }

    public boolean isValid(Network network) {
        if(Networks.hasCycle(network) || !Networks.isDisconnectedNetwork(network, null))return false;
        List<Object> networkLeafNode, networkDFSNode;
        try{
            networkLeafNode = IterableHelp.toList(network.getLeaves());
            networkDFSNode = IterableHelp.toList(network.dfs());
        } catch (Exception ex) {
            return false;
        }
        List<String> leafNodes = new ArrayList();
        for(NetNode leaf : startLeaves) {
            leafNodes.add(leaf.getName());
        }
        int count = 0;
        for(Object leaf: networkLeafNode){
            NetNode l = (NetNode)leaf;
            if(leafNodes.contains(l.getName())){
                count++;
            }else{
                return false;
            }
        }
        if(count!=leafNodes.size()) {
            return false;
        }
        for(Object nodeobj: networkDFSNode){

            NetNode node = (NetNode)nodeobj;

            double totalProb = 0;
            List<NetNode> pars = IterableHelp.toList(node.getParents());
            for (NetNode parent : pars) {
                totalProb += node.getParentProbability(parent);
            }
            if(node.isNetworkNode()){
                double p = node.getParentProbability(pars.get(0));

                if(Math.abs(totalProb - 1) > 0.00001) {
                    if(pars.size() == 2) {

                        if(p < 1.0 && p > 0.0) {
                            node.setParentProbability(pars.get(1), 1 - p);
                        } else {
                            node.setParentProbability(pars.get(0), 1 - node.getParentProbability(pars.get(1)));
                        }
                    } else {
                        return false;
                    }
                }
            }
            else if(!node.isRoot()){
                if(!Double.isNaN(totalProb)){
                    if(pars.size() == 1) {
                        NetNode par = (NetNode) node.getParents().iterator().next();
                        node.setParentProbability(par, node.NO_PROBABILITY);
                    } else {
                        System.out.println("check if at gamma == 1.0");
                        return false;
                    }
                }
            }

        }
        return true;
    }

}
