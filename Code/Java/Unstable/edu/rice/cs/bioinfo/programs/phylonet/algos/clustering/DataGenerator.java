package edu.rice.cs.bioinfo.programs.phylonet.algos.clustering;


import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.ReticulationEdgeAddition;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 7/1/16
 * Time: 11:52 AM
 * To change this template use File | Settings | File Templates.
 */
public class DataGenerator {
    private String _msdir = "/Users/zhujiafan/Documents/Luay/msdir/ms";

    private Set<String> _singleAlleleSpecies;
    private Tuple<NetNode,NetNode> _targetEdge;
    private Tuple<NetNode,NetNode> _sourceEdge;
    private Tuple<NetNode,NetNode> _destinationEdge;
    private double _targetEdgeBrlen = 1;
    private double _targetEdgeInheriProb = 0.5;
    private boolean _printDetails;


    private int _setsPerNetwork;
    private int _numNetwork;
    private Random _random;
    private int [] _sizes;
    Func2<Network, Integer, List> _simulator;

    private void getNetworkInfo(Network network, ArrayList<Tuple<NetNode, NetNode>> allEdges, ArrayList<Tuple<NetNode, NetNode>> edgesNeedBrlens, ArrayList<Tuple<NetNode, NetNode>> allReticulationEdges, ArrayList<Tuple<NetNode, NetNode>> removableReticulationEdges, Set<String> leafSet){

        _singleAlleleSpecies = new HashSet<String>();
        _singleAlleleSpecies.clear();
        for(Object leaf: network.getLeaves()){
            String name = ((NetNode)leaf).getName();
            _singleAlleleSpecies.add(name);
        }

        int numReticulations = 0;
        Map<NetNode,Set<String>> node2leaves = new HashMap<>();
        for(Object nodeO: Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeO;
            Set<String> leaves = new HashSet<>();
            if(node.isLeaf()){
                leafSet.add(node.getName());
                leaves.add(node.getName());
            }
            else if(node.isNetworkNode()){
                numReticulations++;
            }
            for(Object childO: node.getChildren()){
                NetNode childNode = (NetNode)childO;
                Set<String> childLeaves = node2leaves.get(childNode);
                leaves.addAll(childLeaves);
                Tuple<NetNode,NetNode> edge = new Tuple<>(node, childNode);
                allEdges.add(edge);
                if(childLeaves.size()!=1 || !_singleAlleleSpecies.containsAll(childLeaves)){
                    edgesNeedBrlens.add(edge);
                }
                if(childNode.isNetworkNode()) {
                    if (node.isTreeNode()) {
                        removableReticulationEdges.add(edge);
                    }
                    allReticulationEdges.add(edge);
                }

            }
            node2leaves.put(node,leaves);
        }

    }

    private boolean isNetworkValid(Network network, Set<String> leaves){
        int count = 0;
        for(Object leaf: network.getLeaves()){
            if(leaves.contains(((NetNode)leaf).getName())){
                count++;
            }
            else{
                return false;
            }
        }
        if(count!=leaves.size())return false;
        if(!Networks.isDisconnectedNetwork(network,null))return false;
        for(Object node: Networks.postTraversal(network)){
            double totalProb = 0;
            for (Object parent : ((NetNode) node).getParents()) {
                totalProb += ((NetNode) node).getParentProbability((NetNode) parent);
            }
            if(((NetNode)node).getChildCount()==1 && ((NetNode)node).getParentCount()<2){
                return false;
            }
            if(totalProb!=NetNode.NO_PROBABILITY && ((NetNode)node).isNetworkNode()){
                if(Math.abs(totalProb - 1) > 0.00001) {
                    throw new RuntimeException(network.toString());
                }
            }
            else if(!((NetNode)node).isRoot()){
                if(totalProb != NetNode.NO_PROBABILITY){
                    throw new RuntimeException(network.toString());
                }
            }
        }
        return true;
    }

    protected String printEdge(Tuple<NetNode, NetNode> edge) {
        return "(" + edge.Item1.getName() + "," + edge.Item2.getName() + ")";
    }

    private void setParametersForReticulationEdgeAddition(ArrayList<Tuple<NetNode,NetNode>> allEdges){
        _targetEdge = null;
        int size = allEdges.size();

        int sourceID = _random.nextInt(size);
        int destinationID = sourceID;
        while(sourceID==destinationID){
            destinationID = _random.nextInt(size);
        }

        _sourceEdge = allEdges.get(sourceID);
        _destinationEdge = allEdges.get(destinationID);

        if(_printDetails){
            System.out.print("Add " + printEdge(_sourceEdge) + " to " + printEdge(_destinationEdge));
        }

    }

    private boolean isFullReticulation(Network<Object> network){
        //if(true)
            //return true;
        int number = 0;
        for (NetNode<Object> node : network.dfs()) {
            if(node.isNetworkNode())
                number++;
        }
        int count = 0;
        for (NetworkTree<Object> nt : Networks.getTrees(network)) {
            count++;
        }
        if((1 << number) == count)
            return true;
        else
            return false;
    }

    public List<Network> allPossibleReticulation(Network<Object> originalNetwork) {
        List<Network> results = new ArrayList<>();
        ArrayList<Tuple<NetNode, NetNode>> allEdges = new ArrayList<>();
        ArrayList<Tuple<NetNode, NetNode>> allEdgesNeedBrlens = new ArrayList<>();
        ArrayList<Tuple<NetNode, NetNode>> allReticulationEdges = new ArrayList<>();
        ArrayList<Tuple<NetNode, NetNode>> allRemovableReticulationEdges = new ArrayList<>();
        Set<String> taxa = new HashSet<>();

        Network network = originalNetwork.clone();
        //results.add(originalNetwork.clone());
        ReticulationEdgeAddition reticulationEdgeAddition = new ReticulationEdgeAddition();
        getNetworkInfo(network, allEdges, allEdgesNeedBrlens, allReticulationEdges, allRemovableReticulationEdges, taxa);
        int size = allEdges.size();

        for(int sourceID = 0 ; sourceID < size ; sourceID++) {
            for(int destinationID = 0 ; destinationID < size ; destinationID++) {
                if(sourceID == destinationID)
                    continue;

                _sourceEdge = allEdges.get(sourceID);
                _destinationEdge = allEdges.get(destinationID);
                _targetEdge = null;

                reticulationEdgeAddition.setParameters(network, _targetEdge, _targetEdgeBrlen, _targetEdgeInheriProb, _sourceEdge, null, null, _destinationEdge, null, null);
                if (reticulationEdgeAddition.performOperation()) {
                    if (Networks.hasCycle(network) || !isNetworkValid(network, taxa) /*|| !isFullReticulation(network)*/) {

                    } else {
                        results.add(network.clone());
                    }
                    reticulationEdgeAddition.undoOperation();
                }
            }
        }



        return results;


    }

    public List<Network> addReticulations(Network<Object> originalNetwork, int n) {
        List<Network> results = new ArrayList<>();
        for(int i = 0 ; i < n ; i++) {
            Network net = originalNetwork.clone();
            addReticulation(net);
            results.add(net);
        }
        return results;
    }

    public void addReticulation(Network<Object> network){
        ArrayList<Tuple<NetNode, NetNode>> allEdges = new ArrayList<>();
        ArrayList<Tuple<NetNode, NetNode>> allEdgesNeedBrlens = new ArrayList<>();
        ArrayList<Tuple<NetNode, NetNode>> allReticulationEdges = new ArrayList<>();
        ArrayList<Tuple<NetNode, NetNode>> allRemovableReticulationEdges = new ArrayList<>();
        Set<String> taxa = new HashSet<>();

        getNetworkInfo(network, allEdges, allEdgesNeedBrlens, allReticulationEdges, allRemovableReticulationEdges, taxa);

        /*for(Tuple<NetNode, NetNode> edge : allEdges) {
            if(edge.Item1 == network.getRoot())
                allEdges.remove(edge);
        }*/

        ReticulationEdgeAddition reticulationEdgeAddition = new ReticulationEdgeAddition();

        boolean successRearrangment;
        do{
            Set<Integer> previousTriedEdges = new HashSet<>();
            successRearrangment = false;
            while(!successRearrangment){
                setParametersForReticulationEdgeAddition(allEdges);
                reticulationEdgeAddition.setParameters(network, _targetEdge, _targetEdgeBrlen, _targetEdgeInheriProb, _sourceEdge, null, null, _destinationEdge, null, null);

                if (reticulationEdgeAddition.performOperation()) {

                    if (Networks.hasCycle(network) || !isNetworkValid(network, taxa) /*|| !isFullReticulation(network)*/) {
                        successRearrangment = false;
                        reticulationEdgeAddition.undoOperation();
                    }
                    else{
                        successRearrangment = true;
                        if(_printDetails){
                            System.out.println();
                        }
                    }
                }
            }
        }while(!successRearrangment);
    }

    public Network<Object> getRandomNetwork(int taxa) {
        String[] leaves;
        leaves = new String[taxa];
        for(int i = 0 ; i < taxa ; i++) {
            leaves[i] = Character.toString((char)('A' + i));
        }
        int maxReticulations = 1;

        Tree randomTree = Trees.generateRandomTree(Arrays.copyOfRange(leaves, 0, taxa));
        Network<Object> trueNetwork = Networks.readNetwork(randomTree.toNewick());

        _printDetails = false;
        _targetEdgeBrlen = 1;
        _targetEdgeInheriProb = 0.5;

        _singleAlleleSpecies = new HashSet<String>();
        _singleAlleleSpecies.clear();
        for(Object leaf: trueNetwork.getLeaves()){
            String name = ((NetNode)leaf).getName();
            _singleAlleleSpecies.add(name);
        }

        for(int j = 1 ; j <= maxReticulations ; j++){
            _targetEdgeBrlen = 0;
            addReticulation(trueNetwork);
        }

        return trueNetwork;
    }

    public Network<Object> getNewNetwork1() {
        List<Double> blPool = new ArrayList<>();
        for (int i = 0; i < 12; i++)
            blPool.add(0.7 + _random.nextDouble() * (1.3 - 0.7));

        blPool.set(6, blPool.get(0) + blPool.get(1) + blPool.get(2) - blPool.get(4) - blPool.get(5));
        double bl2 = blPool.get(0) + blPool.get(1) + blPool.get(2);
        double bl6 = blPool.get(4) + blPool.get(5) + blPool.get(6);
        blPool.set(4, blPool.get(4) / bl6 * bl2);
        blPool.set(5, blPool.get(5) / bl6 * bl2);
        blPool.set(6, blPool.get(6) / bl6 * bl2);

        double blTotal = 8.0;
        double blK = 0.2; //t
        double blAB = blTotal - blPool.get(0) - blPool.get(1) - blPool.get(2) - blPool.get(3);
        double blG = blTotal - blPool.get(0) - blPool.get(1);
        double blH = blTotal - blPool.get(0);
        double blI = blTotal - blPool.get(4);
        double blJ = blTotal - blPool.get(4) - blPool.get(5);
        double blEF = blTotal - blPool.get(4) - blPool.get(5) - blPool.get(6) - blPool.get(7);
        double blCD = blTotal - blPool.get(0) - blPool.get(1) - blPool.get(2) - blK;

        double prob1 = 0.35;
        double prob2 = 1 - prob1;
        String newick = "(((((A:" + blAB + ",B:" + blAB + "):" + blPool.get(3) + ",((C:" + blCD + ",D:" + blCD + ")K:" + blK + ")R#H1:0.0::" + prob1 + "):" + blPool.get(2) + ",G:" + blG + "):" + blPool.get(1) + ",H:" + blH + "):" + blPool.get(0) + ",((((E:" + blEF + ",F:" + blEF + "):" + blPool.get(7) + ",R#H1:0.0::" + prob2 + "):" + blPool.get(6) + ",J:" + blJ + "):" + blPool.get(5) + ",I:" + blI + "):" + blPool.get(4) + ");";

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork(newick);

        return trueNetwork;
    }

    public Network<Object> getNewNetworkWithXY() {
        int numXY = 6;
        List<Double> blPool = new ArrayList<>();
        List<String> leafPool = new ArrayList<>();

        for(int i = 0 ; i < 12 ; i++)
            leafPool.add("" + (char)('W' - i));

        for (int i = 0; i < 12; i++)
            blPool.add(0.7 + _random.nextDouble() * (1.3 - 0.7));

        double sumXY = 0;

        for(int i = 0 ; i < numXY + 2 ; i++) {
            sumXY += blPool.get(i);
        }

        sumXY /= 2;
        double sumX = 0;
        for(int i = 0 ; i < numXY / 2 + 1 ; i++) {
            sumX += blPool.get(i);
        }
        for(int i = 0 ; i < numXY / 2 + 1 ; i++) {
            blPool.set(i, blPool.get(i) / sumX * sumXY);
        }

        double sumY = 0;
        for(int i = numXY / 2 + 1 ; i < numXY + 2 ; i++) {
            sumY += blPool.get(i);
        }
        for(int i = numXY / 2 + 1 ; i < numXY + 2 ; i++) {
            blPool.set(i, blPool.get(i) / sumY * sumXY);
        }

        double blTotal = 8.0;
        double blK = 0.2;
        double blJ = blPool.get(numXY + 3);
        double blL = blPool.get(numXY + 4);
        double prob1 = 0.35;
        double prob2 = 1 - prob1;


        Network<Object> trueNetwork = Networks.readNetwork("Root;");
        NetNode<Object> node = trueNetwork.getRoot();

        sumX = 0;
        for(int i = 0 ; i < numXY / 2 + 1 ; i++) {
            NetNode<Object> newnode = new BniNetNode<Object>();
            node.adoptChild(newnode, blPool.get(i));
            sumX += blPool.get(i);
            NetNode<Object> leafnode = new BniNetNode<Object>();
            if(i + 1 == numXY / 2 + 1) {
                newnode.setName("X");
            } else {
                leafnode.setName(leafPool.get(i));
                newnode.adoptChild(leafnode, blTotal - sumX);
            }
            node = newnode;
        }
        NetNode<Object> nodeX = node;

        sumY = 0;
        node = trueNetwork.getRoot();
        for(int i = numXY / 2 + 1 ; i < numXY + 2 ; i++) {
            NetNode<Object> newnode = new BniNetNode<Object>();
            node.adoptChild(newnode, blPool.get(i));
            sumY += blPool.get(i);
            NetNode<Object> leafnode = new BniNetNode<Object>();
            if(i + 1 == numXY + 2) {
                newnode.setName("Y");
            } else {
                leafnode.setName(leafPool.get(i));
                newnode.adoptChild(leafnode, blTotal - sumY);
            }
            node = newnode;
        }
        NetNode<Object> nodeY = node;

        node = new BniNetNode<Object>();
        node.setName("I");

        nodeX.adoptChild(node, 0.0);
        node.setParentProbability(nodeX, prob1);
        nodeY.adoptChild(node, 0.0);
        node.setParentProbability(nodeY, prob2);

        NetNode<Object> nodeK = new BniNetNode<Object>();
        node.adoptChild(nodeK, blK);
        NetNode<Object> nodeC = new BniNetNode<Object>();
        NetNode<Object> nodeD = new BniNetNode<Object>();
        nodeK.adoptChild(nodeC, blTotal - blK - sumXY);
        nodeK.adoptChild(nodeD, blTotal - blK - sumXY);
        nodeK.setName("K");
        nodeC.setName("C");
        nodeD.setName("D");

        NetNode<Object> nodeJ = new BniNetNode<Object>();
        nodeX.adoptChild(nodeJ, blJ);
        NetNode<Object> nodeA = new BniNetNode<Object>();
        NetNode<Object> nodeB = new BniNetNode<Object>();
        nodeJ.adoptChild(nodeA, blTotal - blJ - sumXY);
        nodeJ.adoptChild(nodeB, blTotal - blJ - sumXY);
        nodeJ.setName("J");
        nodeA.setName("A");
        nodeB.setName("B");

        NetNode<Object> nodeL = new BniNetNode<Object>();
        nodeY.adoptChild(nodeL, blL);
        NetNode<Object> nodeE = new BniNetNode<Object>();
        NetNode<Object> nodeF = new BniNetNode<Object>();
        nodeL.adoptChild(nodeE, blTotal - blL - sumXY);
        nodeL.adoptChild(nodeF, blTotal - blL - sumXY);
        nodeL.setName("L");
        nodeE.setName("E");
        nodeF.setName("F");

        return trueNetwork;
    }

    public Network<Object> getNewNetworkWithAlpha() {
        double alpha = 7.0;
        int numXY = 0;
        List<Double> blPool = new ArrayList<>();
        List<String> leafPool = new ArrayList<>();

        for(int i = 0 ; i < 12 ; i++)
            leafPool.add("" + (char)('W' - i));

        for (int i = 0; i < 12; i++)
            blPool.add(0.7 + _random.nextDouble() * (1.3 - 0.7));

        double sumXY = 0;

        for(int i = 0 ; i < numXY + 2 ; i++) {
            sumXY += blPool.get(i);
        }
        sumXY = alpha;

        sumXY /= 2;
        double sumX = 0;
        for(int i = 0 ; i < numXY / 2 + 1 ; i++) {
            sumX += blPool.get(i);
        }
        for(int i = 0 ; i < numXY / 2 + 1 ; i++) {
            blPool.set(i, blPool.get(i) / sumX * sumXY);
        }

        double sumY = 0;
        for(int i = numXY / 2 + 1 ; i < numXY + 2 ; i++) {
            sumY += blPool.get(i);
        }
        for(int i = numXY / 2 + 1 ; i < numXY + 2 ; i++) {
            blPool.set(i, blPool.get(i) / sumY * sumXY);
        }

        double blTotal = 8.0;
        double blK = 0.2;
        double blJ = blPool.get(numXY + 3);
        double blL = blPool.get(numXY + 4);
        double prob1 = 0.35;
        double prob2 = 1 - prob1;


        Network<Object> trueNetwork = Networks.readNetwork("Root;");
        NetNode<Object> node = trueNetwork.getRoot();

        sumX = 0;
        for(int i = 0 ; i < numXY / 2 + 1 ; i++) {
            NetNode<Object> newnode = new BniNetNode<Object>();
            node.adoptChild(newnode, blPool.get(i));
            sumX += blPool.get(i);
            NetNode<Object> leafnode = new BniNetNode<Object>();
            if(i + 1 == numXY / 2 + 1) {
                newnode.setName("X");
            } else {
                leafnode.setName(leafPool.get(i));
                newnode.adoptChild(leafnode, blTotal - sumX);
            }
            node = newnode;
        }
        NetNode<Object> nodeX = node;

        sumY = 0;
        node = trueNetwork.getRoot();
        for(int i = numXY / 2 + 1 ; i < numXY + 2 ; i++) {
            NetNode<Object> newnode = new BniNetNode<Object>();
            node.adoptChild(newnode, blPool.get(i));
            sumY += blPool.get(i);
            NetNode<Object> leafnode = new BniNetNode<Object>();
            if(i + 1 == numXY + 2) {
                newnode.setName("Y");
            } else {
                leafnode.setName(leafPool.get(i));
                newnode.adoptChild(leafnode, blTotal - sumY);
            }
            node = newnode;
        }
        NetNode<Object> nodeY = node;

        node = new BniNetNode<Object>();
        node.setName("I");

        nodeX.adoptChild(node, 0.0);
        node.setParentProbability(nodeX, prob1);
        nodeY.adoptChild(node, 0.0);
        node.setParentProbability(nodeY, prob2);

        NetNode<Object> nodeK = new BniNetNode<Object>();
        node.adoptChild(nodeK, blK);
        NetNode<Object> nodeC = new BniNetNode<Object>();
        NetNode<Object> nodeD = new BniNetNode<Object>();
        nodeK.adoptChild(nodeC, blTotal - blK - sumXY);
        nodeK.adoptChild(nodeD, blTotal - blK - sumXY);
        nodeK.setName("K");
        nodeC.setName("C");
        nodeD.setName("D");

        NetNode<Object> nodeJ = new BniNetNode<Object>();
        nodeX.adoptChild(nodeJ, blJ);
        NetNode<Object> nodeA = new BniNetNode<Object>();
        NetNode<Object> nodeB = new BniNetNode<Object>();
        nodeJ.adoptChild(nodeA, blTotal - blJ - sumXY);
        nodeJ.adoptChild(nodeB, blTotal - blJ - sumXY);
        nodeJ.setName("J");
        nodeA.setName("A");
        nodeB.setName("B");

        NetNode<Object> nodeL = new BniNetNode<Object>();
        nodeY.adoptChild(nodeL, blL);
        NetNode<Object> nodeE = new BniNetNode<Object>();
        NetNode<Object> nodeF = new BniNetNode<Object>();
        nodeL.adoptChild(nodeE, blTotal - blL - sumXY);
        nodeL.adoptChild(nodeF, blTotal - blL - sumXY);
        nodeL.setName("L");
        nodeE.setName("E");
        nodeF.setName("F");

        return trueNetwork;
    }

    public Network<Object> getNewNetworkR2_1() {
        List<Double> blPool = new ArrayList<>();

        for (int i = 0; i < 12; i++)
            blPool.add(0.7 + _random.nextDouble() * (1.3 - 0.7));

        double blTotal = 8.0;
        double blI1 = 2.0;
        double blJ1 = 2.0;
        double blX = blPool.get(0);
        double blY = blPool.get(1);
        double blS = blPool.get(2);
        double blT = blPool.get(3);
        double blV = blPool.get(4);
        double blW = blPool.get(5);
        double blO = blPool.get(6);
        double blP = blPool.get(7);
        double blQ = blPool.get(8);
        double blR = blPool.get(9);
        double prob1 = 0.5;
        double prob2 = 1 - prob1;
        double prob3 = 0.5;
        double prob4 = 1 - prob3;


        Network<Object> trueNetwork = Networks.readNetwork("Root;");

        NetNode<Object> nodeX = new BniNetNode<Object>();
        trueNetwork.getRoot().adoptChild(nodeX, blX);

        NetNode<Object> nodeY = new BniNetNode<Object>();
        trueNetwork.getRoot().adoptChild(nodeY, blY);

        double sumOP = (blS + blT + blO + blP) / 2;
        double sumO = blS + blO;
        double sumP = blT + blP;
        blS = blS / sumO * sumOP;
        blO = blO / sumO * sumOP;
        blT = blT / sumP * sumOP;
        blP = blP / sumP * sumOP;

        NetNode<Object> nodeS = new BniNetNode<Object>();
        nodeX.adoptChild(nodeS, blS);
        NetNode<Object> nodeK = new BniNetNode<Object>();
        nodeK.setName("K");
        nodeS.adoptChild(nodeK, blTotal - blS - blX);
        NetNode<Object> nodeO = new BniNetNode<Object>();
        nodeS.adoptChild(nodeO, blO);


        NetNode<Object> nodeT = new BniNetNode<Object>();
        nodeX.adoptChild(nodeT, blT);
        NetNode<Object> nodeL = new BniNetNode<Object>();
        nodeL.setName("L");
        nodeT.adoptChild(nodeL, blTotal - blT - blX);
        NetNode<Object> nodeP = new BniNetNode<Object>();
        nodeT.adoptChild(nodeP, blP);

        NetNode<Object> nodeI = new BniNetNode<Object>();
        nodeO.adoptChild(nodeI, 0.0);
        nodeP.adoptChild(nodeI, 0.0);
        nodeI.setParentProbability(nodeO, prob1);
        nodeI.setParentProbability(nodeP, prob2);
        NetNode<Object> nodeI1 = new BniNetNode<Object>();
        nodeI.adoptChild(nodeI1, blI1);

        NetNode<Object> nodeA = new BniNetNode<Object>();
        nodeA.setName("A");
        NetNode<Object> nodeB = new BniNetNode<Object>();
        nodeB.setName("B");
        NetNode<Object> nodeC = new BniNetNode<Object>();
        nodeC.setName("C");
        NetNode<Object> nodeD = new BniNetNode<Object>();
        nodeD.setName("D");

        nodeO.adoptChild(nodeA, blTotal - sumOP - blX);
        nodeI1.adoptChild(nodeB, blTotal - sumOP - blX - blI1);
        nodeI1.adoptChild(nodeC, blTotal - sumOP - blX - blI1);
        nodeP.adoptChild(nodeD, blTotal - sumOP - blX);




        double sumQR = (blV + blQ + blW + blR) / 2;
        double sumQ = blV + blQ;
        double sumR = blW + blR;
        blV = blV / sumQ * sumQR;
        blQ = blQ / sumQ * sumQR;
        blW = blW / sumR * sumQR;
        blR = blR / sumR * sumQR;

        NetNode<Object> nodeV = new BniNetNode<Object>();
        nodeY.adoptChild(nodeV, blV);
        NetNode<Object> nodeM = new BniNetNode<Object>();
        nodeM.setName("M");
        nodeV.adoptChild(nodeM, blTotal - blV - blY);
        NetNode<Object> nodeQ = new BniNetNode<Object>();
        nodeV.adoptChild(nodeQ, blQ);


        NetNode<Object> nodeW = new BniNetNode<Object>();
        nodeY.adoptChild(nodeW, blW);
        NetNode<Object> nodeN = new BniNetNode<Object>();
        nodeN.setName("N");
        nodeW.adoptChild(nodeN, blTotal - blW - blY);
        NetNode<Object> nodeR = new BniNetNode<Object>();
        nodeW.adoptChild(nodeR, blR);

        NetNode<Object> nodeJ = new BniNetNode<Object>();
        nodeQ.adoptChild(nodeJ, 0.0);
        nodeR.adoptChild(nodeJ, 0.0);
        nodeJ.setParentProbability(nodeQ, prob3);
        nodeJ.setParentProbability(nodeR, prob4);
        NetNode<Object> nodeJ1 = new BniNetNode<Object>();
        nodeJ.adoptChild(nodeJ1, blJ1);

        NetNode<Object> nodeE = new BniNetNode<Object>();
        nodeE.setName("E");
        NetNode<Object> nodeF = new BniNetNode<Object>();
        nodeF.setName("F");
        NetNode<Object> nodeG = new BniNetNode<Object>();
        nodeG.setName("G");
        NetNode<Object> nodeH = new BniNetNode<Object>();
        nodeH.setName("H");

        nodeQ.adoptChild(nodeE, blTotal - sumQR - blY);
        nodeJ1.adoptChild(nodeF, blTotal - sumQR - blY - blJ1);
        nodeJ1.adoptChild(nodeG, blTotal - sumQR - blY - blJ1);
        nodeR.adoptChild(nodeH, blTotal - sumQR - blY);



        return trueNetwork;
    }

    public Network<Object> getNewNetworkR2_2() {
        List<Double> blPool = new ArrayList<>();

        for (int i = 0; i < 12; i++)
            blPool.add(0.7 + _random.nextDouble() * (1.3 - 0.7));

        double blTotal = 8.0;
        double blV = blPool.get(0);
        double blX = blPool.get(1);
        double blW = blPool.get(2);
        double blY = blPool.get(3);

        double blL = blPool.get(4);
        double blO = blPool.get(5);
        double blP = blPool.get(6);
        double blM = blPool.get(7);

        double blJ = 0.2;
        double blK = blL + blO + blP - blJ;

        double blR = 0.2;

        double prob1 = 0.35;
        double prob2 = 1 - prob1;
        double prob3 = 0.2;
        double prob4 = 1 - prob3;

        double sumXY = (blV + blX + blW + blY) / 2;
        double sumX = blV + blX;
        double sumY = blW + blY;
        blV = blV / sumX * sumXY;
        blX = blX / sumX * sumXY;
        blW = blW / sumY * sumXY;
        blY = blY / sumY * sumXY;

        Network<Object> trueNetwork = Networks.readNetwork("Root;");

        NetNode<Object> nodeV = new BniNetNode<Object>();
        trueNetwork.getRoot().adoptChild(nodeV, blV);
        NetNode<Object> nodeT = new BniNetNode<Object>();
        nodeT.setName("T");
        nodeV.adoptChild(nodeT, blTotal - blV);
        NetNode<Object> nodeX = new BniNetNode<Object>();
        nodeV.adoptChild(nodeX, blX);

        NetNode<Object> nodeW = new BniNetNode<Object>();
        trueNetwork.getRoot().adoptChild(nodeW, blW);
        NetNode<Object> nodeU = new BniNetNode<Object>();
        nodeU.setName("U");
        nodeW.adoptChild(nodeU, blTotal - blW);
        NetNode<Object> nodeY = new BniNetNode<Object>();
        nodeW.adoptChild(nodeY, blY);

        NetNode<Object> nodeI = new BniNetNode<Object>();
        nodeX.adoptChild(nodeI, 0.0);
        nodeY.adoptChild(nodeI, 0.0);
        nodeI.setParentProbability(nodeX, prob1);
        nodeI.setParentProbability(nodeY, prob2);

        NetNode<Object> nodeJ = new BniNetNode<Object>();
        nodeI.adoptChild(nodeJ, blJ);
        NetNode<Object> nodeK = new BniNetNode<Object>();
        nodeJ.adoptChild(nodeK, blK);
        NetNode<Object> nodeE = new BniNetNode<Object>();
        nodeE.setName("E");
        nodeK.adoptChild(nodeE, blTotal - blK - blJ - sumXY);
        NetNode<Object> nodeF = new BniNetNode<Object>();
        nodeF.setName("F");
        nodeJ.adoptChild(nodeF, blTotal - blJ - sumXY);

        NetNode<Object> nodeM = new BniNetNode<Object>();
        nodeY.adoptChild(nodeM, blM);
        NetNode<Object> nodeG = new BniNetNode<Object>();
        nodeG.setName("G");
        nodeM.adoptChild(nodeG, blTotal - blM - sumXY);
        NetNode<Object> nodeH = new BniNetNode<Object>();
        nodeH.setName("H");
        nodeM.adoptChild(nodeH, blTotal - blM - sumXY);

        NetNode<Object> nodeL = new BniNetNode<Object>();
        nodeX.adoptChild(nodeL, blL);
        NetNode<Object> nodeN = new BniNetNode<Object>();
        nodeN.setName("N");
        nodeL.adoptChild(nodeN, blTotal - blL - sumXY);
        NetNode<Object> nodeO = new BniNetNode<Object>();
        nodeL.adoptChild(nodeO, blO);

        NetNode<Object> nodeA = new BniNetNode<Object>();
        nodeA.setName("A");
        nodeO.adoptChild(nodeA, blTotal - blO - blL - sumXY);
        NetNode<Object> nodeP = new BniNetNode<Object>();
        nodeO.adoptChild(nodeP, blP);
        NetNode<Object> nodeB = new BniNetNode<Object>();
        nodeB.setName("B");
        nodeP.adoptChild(nodeB, blTotal - blP - blO - blL - sumXY);


        NetNode<Object> nodeQ = new BniNetNode<Object>();
        nodeP.adoptChild(nodeQ, 0.0);
        nodeK.adoptChild(nodeQ, 0.0);
        nodeQ.setParentProbability(nodeP, prob3);
        nodeQ.setParentProbability(nodeK, prob4);

        NetNode<Object> nodeR = new BniNetNode<Object>();
        nodeQ.adoptChild(nodeR, blR);
        NetNode<Object> nodeC = new BniNetNode<Object>();
        nodeC.setName("C");
        nodeR.adoptChild(nodeC, blTotal - blR - blP - blO - blL - sumXY);
        NetNode<Object> nodeD = new BniNetNode<Object>();
        nodeD.setName("D");
        nodeR.adoptChild(nodeD, blTotal - blR - blP - blO - blL - sumXY);



        return trueNetwork;
    }

    public Network<Object> getNewNetwork2() { //8 parental trees
        List<Double> blPool = new ArrayList<>();
        for (int i = 0; i < 12; i++)
            blPool.add(0.7 + _random.nextDouble() * (1.3 - 0.7));

        blPool.set(6, blPool.get(0) + blPool.get(1) + blPool.get(2) - blPool.get(4) - blPool.get(5));
        double bl2 = blPool.get(0) + blPool.get(1) + blPool.get(2);
        double bl6 = blPool.get(4) + blPool.get(5) + blPool.get(6);
        blPool.set(4, blPool.get(4) / bl6 * bl2);
        blPool.set(5, blPool.get(5) / bl6 * bl2);
        blPool.set(6, blPool.get(6) / bl6 * bl2);

        double blTotal = 9.0;
        double blK = 0.2;
        double blAB = blTotal - blPool.get(0) - blPool.get(1) - blPool.get(2) - blPool.get(3);
        double blG = blTotal - blPool.get(0) - blPool.get(1);
        double blH = blTotal - blPool.get(0);
        double blI = blTotal - blPool.get(4);
        double blJ = blTotal - blPool.get(4) - blPool.get(5);
        double blEF = blTotal - blPool.get(4) - blPool.get(5) - blPool.get(6) - blPool.get(7);
        double blD = blTotal - blPool.get(0) - blPool.get(1) - blPool.get(2) - blK;
        double blY = 0.1;
        double blCZ = blD - blY;

        double prob1 = 0.35;
        double prob2 = 1 - prob1;
        String newick = "(((((A:" + blAB + ",B:" + blAB + "):" + blPool.get(3) + ",(((C:" + blCZ + ",Z:" + blCZ + ")Y:" + blY +",D:" + blD + ")K:" + blK + ")R#H1:0.0::" + prob1 + "):" + blPool.get(2) + ",G:" + blG + "):" + blPool.get(1) + ",H:" + blH + "):" + blPool.get(0) + ",((((E:" + blEF + ",F:" + blEF + "):" + blPool.get(7) + ",R#H1:0.0::" + prob2 + "):" + blPool.get(6) + ",J:" + blJ + "):" + blPool.get(5) + ",I:" + blI + "):" + blPool.get(4) + ");";

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork(newick);

        return trueNetwork;
    }

    public Network<Object> getNewNetwork() {

        Network<Object> network = getNewNetworkR2_1();

        for(NetNode<Object> node: network.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                if(node.getParentDistance(parent) != NetNode.NO_DISTANCE && node.getParentDistance(parent) < 0)
                    throw new RuntimeException("Negative branch length!");
            }
        }

        System.out.println(network.toString());

        return network;
    }

    public DataGenerator() {
        _numNetwork = 1;
        _setsPerNetwork = 1;
        _simulator = getSimulator(null);
        _random = new Random();
        _sizes = new int[]{1000};
    }

    public DataGenerator(int numNetwork, int [] sizes, int setsPerNetwork) {
        _numNetwork = numNetwork;
        _setsPerNetwork = setsPerNetwork;
        _simulator = getSimulator(null);
        _random = new Random();
        _sizes = sizes;
    }

    public void generateFiles(String name) {
        new File("./" + name).mkdir();
        try {
            PrintWriter out = new PrintWriter("./" + name + "/" + name + ".txt");
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        for(int i = 0 ; i < _numNetwork ; i++) {
            Network<Object> trueNetwork;
            trueNetwork = getNewNetwork();

            for(int j = 0 ; j < _setsPerNetwork ; j++) {

                String path = "./" + name + "/" + i + "-" + j;
                new File(path).mkdir();

                try {
                    PrintWriter out = new PrintWriter(path + "/network.txt");
                    out.println(trueNetwork.toString());
                    out.close();
                } catch (IOException ioe) {
                    ioe.printStackTrace();
                }

                List<List<MutableTuple<Tree, Double>>> simulatedGTs = null;
                for (int sizeIndex = 0; sizeIndex < _sizes.length; sizeIndex++) {
                    int size = _sizes[sizeIndex];

                    if (sizeIndex == 0) {
                        simulatedGTs = _simulator.execute(trueNetwork, size);
                    } else {
                        simulatedGTs.addAll(_simulator.execute(trueNetwork, size - _sizes[sizeIndex - 1]));
                    }

                    try {
                        PrintWriter out = new PrintWriter(path + "/" + size + ".txt");
                        for (int k = 0; k < size; k++)
                            out.println(simulatedGTs.get(k).get(0).Item1);
                        out.close();
                    }catch(IOException ioe) {
                        ioe.printStackTrace();
                    }

                }
            }

        }
    }

    void generateParentalTreesNexus(){
        Network network = getNewNetworkR2_2();

    }


    protected Func2<Network, Integer, List> getSimulator(final Map<String, List<String>> species2alleles) {
        return new Func2<Network, Integer, List>() {
            public List execute(Network network, Integer numGTs) {
                SimGTInNetworkByMS simulator = new SimGTInNetworkByMS();
                //SimGTInNetwork simulator = new SimGTInNetwork();
                List<List<MutableTuple>> gts = new ArrayList<List<MutableTuple>>();
                for (Tree tr : simulator.generateGTs(network, species2alleles, numGTs, _msdir)) {
                    //for (Tree tr : simulator.generateGTs(network, species2alleles, numGTs)) {
                    gts.add(Arrays.asList(new MutableTuple(tr, 1.0)));
                }

                /*SimGTInNetwork simulator = new SimGTInNetwork();
                List<List<MutableTuple>> gts = new ArrayList<List<MutableTuple>>();
                for (Tree tr : simulator.generateGTs(network, species2alleles, numGTs)) {
                    gts.add(Arrays.asList(new MutableTuple(tr, 1.0)));
                }*/
                return gts;
            }
        };
    }
}
