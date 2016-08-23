package edu.rice.cs.bioinfo.programs.phylonet.algos.clustering;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import javafx.util.Pair;

import java.util.*;
import java.util.concurrent.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 7/21/16
 * Time: 10:57 AM
 * To change this template use File | Settings | File Templates.
 */
public class ParentalTreeOperation {
    private STITree getSubtree(STITree tree, List<String> leaves) {
        Network<Object> network = Networks.readNetwork(tree.toNewick());

        boolean done = false;
        while(!done) {
            done = true;
            for (NetNode<Object> node : network.dfs()) {
                if (!leaves.contains(node.getName()) && node.getChildCount() == 0) {
                    for (NetNode<Object> parent : node.getParents()) {
                        parent.removeChild(node);
                        done = false;
                    }
                }
            }
        }

        STITree subtree = null;
        try {
            subtree = new STITree(network.toString());
        } catch(Exception e) {
            e.printStackTrace();
        }
        return subtree;
    }

    private void removeBinaryNodes(Network<Double> net)
    {
        // Find all binary nodes.
        List<NetNode<Double>> binaryNodes = new LinkedList<NetNode<Double>>();
        for (NetNode<Double> node : net.bfs()) {
            if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                binaryNodes.add(node);
            }
        }

        // Remove them.
        for (NetNode<Double> node : binaryNodes) {
            NetNode<Double> child = node.getChildren().iterator().next();	// Node's only child.
            if(child.getIndeg() != 1){
                continue;
            }
            NetNode<Double> parent = node.getParents().iterator().next();	// Node's only parent.
            double distance = node.getParentDistance(parent) + child.getParentDistance(node);
            double gamma1 = node.getParentProbability(parent)==NetNode.NO_PROBABILITY?1.0:node.getParentProbability(parent);
            double gamma2 = child.getParentProbability(node)==NetNode.NO_PROBABILITY?1.0:child.getParentProbability(node);
            double gamma =  gamma1 * gamma2;
            parent.removeChild(node);
            node.removeChild(child);
            parent.adoptChild(child, distance);
            child.setParentProbability(parent, gamma);
        }
    }

    public List<Tree> getParentalTrees(Network<Object> network){
        Network<Double> net = Networks.readNetwork(network.toString());
        STITree mulTree;
        List<String> stTaxa = new ArrayList<>();
        Map<String,Integer> nname2tamount = new TreeMap<>();
        removeBinaryNodes(net);
        mulTree = new STITree<Double>();
        ((STINode<Double>)(mulTree.getRoot())).setData(1.0);
        ((STINode<Double>)(mulTree.getRoot())).setName("root");
        Queue<NetNode<Double>> source = new LinkedList<NetNode<Double>>();
        Queue<TMutableNode> dest = new LinkedList<TMutableNode>();
        source.offer(net.getRoot());
        dest.offer((TMutableNode) mulTree.getRoot());
        int nameid = 1;
        //long nameid = 0;
        while(!source.isEmpty()){
            NetNode<Double> parent = source.poll();
            TMutableNode peer = dest.poll();
            for (NetNode<Double> child : parent.getChildren()) {
                TMutableNode copy;
                if (child.getName().equals(NetNode.NO_NAME)) {
                    child.setName("hnode" + (nameid++));
                }
                String name = child.getName();
                if(child.isNetworkNode()){
                    name = child.getName()+"TO"+parent.getName();
                }
                Integer amount = nname2tamount.get(name);
                if(amount==null){
                    amount = 0;
                }
                nname2tamount.put(name, ++amount);
                String newname = name + "_" + amount;
                copy = peer.createChild(newname);
                if(child.isLeaf()) {
                    stTaxa.add(newname);
                }

                // Update the distance and data for this child.
                double distance = child.getParentDistance(parent);
                if (distance == NetNode.NO_DISTANCE) {
                    //copy.setParentDistance(TNode.NO_DISTANCE);
                    copy.setParentDistance(0);
                }
                else {
                    copy.setParentDistance(distance);
                }

                double gamma = child.getParentProbability(parent);
                gamma = gamma==NetNode.NO_PROBABILITY?1.0:gamma;
                ((STINode<Double>)copy).setData(gamma);

                // Continue to iterate over the children of nn and tn.
                source.offer(child);
                dest.offer(copy);
                //index ++;
            }
        }

        List<String> leaves = new ArrayList<>();
        Stack<Pair<String, Integer>> stack = new Stack<>();
        for(Map.Entry<String, Integer> entry : nname2tamount.entrySet()) {
            if(stTaxa.contains(entry.getKey() + "_1")) {
                if (entry.getValue() > 1) {
                    stack.push(new Pair<String, Integer>(entry.getKey(), 1));
                }
                leaves.add(entry.getKey() + "_1");
            }
        }

        List<STITree<Double>> parentalTrees = new ArrayList<>();

        int k = stack.size() - 1;
        while(k >= 0) {
            parentalTrees.add(getSubtree(mulTree, leaves));
            leaves.remove(stack.get(k).getKey() + "_" + stack.get(k).getValue());
            stack.set(k, new Pair<String, Integer>(stack.get(k).getKey(), stack.get(k).getValue() + 1));
            leaves.add(stack.get(k).getKey() + "_" + stack.get(k).getValue());
            while(k > 0 && stack.get(k).getValue() > nname2tamount.get(stack.get(k).getKey())) {
                leaves.remove(stack.get(k).getKey() + "_" + stack.get(k).getValue());
                stack.set(k, new Pair<String, Integer>(stack.get(k).getKey(), 1));
                leaves.add(stack.get(k).getKey() + "_" + stack.get(k).getValue());
                k--;
                leaves.remove(stack.get(k).getKey() + "_" + stack.get(k).getValue());
                stack.set(k, new Pair<String, Integer>(stack.get(k).getKey(), stack.get(k).getValue() + 1));
                leaves.add(stack.get(k).getKey() + "_" + stack.get(k).getValue());
            }
            if(k == 0 && stack.get(k).getValue() > nname2tamount.get(stack.get(k).getKey()))
                break;
            k = stack.size() - 1;
        }

        List<Tree> newParentalTrees = new ArrayList<>();

        for(STITree<Double> tree : parentalTrees) {
            for(STINode<Double> node : tree.getNodes()) {
                if(node.isLeaf()) {
                    String name = node.getName();
                    node.setName(name.substring(0, name.lastIndexOf("_")));
                }
            }
            newParentalTrees.add(Trees.readTree(tree.toNewick()));
        }

        return newParentalTrees;
    }

    public void test() {
        int number = 100;
        int correctness = 0;
        double remainRatio = 0.5;

        DataGenerator dataGenerator = new DataGenerator();
        ParentalTreeOperation parentalTreeOperation = new ParentalTreeOperation();
        InferNetworkFromParentalTrees inferNetworkFromParentalTrees = new InferNetworkFromParentalTrees();
        StringBuilder result = new StringBuilder();

        for(int i = 0 ; i < number ; i++) {
            Network<Object> trueNetwork = dataGenerator.getNewNetworkR2_2();

            for(NetNode<Object> node: trueNetwork.dfs()){
                for(NetNode<Object> parent: node.getParents()){
                    node.setParentDistance(parent,NetNode.NO_DISTANCE);
                    if(node.isNetworkNode()){
                        if(node.getParentProbability(parent) == NetNode.NO_PROBABILITY)
                            node.setParentProbability(parent, 0.5);
                    }
                }
            }
            System.out.println("True network: " + trueNetwork.toString());

            List<Tree> parentalTrees = parentalTreeOperation.getParentalTrees(trueNetwork);
            Collections.shuffle(parentalTrees);
            parentalTrees = parentalTrees.subList(0, (int)(parentalTrees.size() * remainRatio));

            System.out.println("Number of input parental trees: " + parentalTrees.size());

            for(Tree tree : parentalTrees) {

                for(TNode node : tree.postTraverse()) {
                    STINode stiNode = (STINode) node;
                    stiNode.setParentDistance(TNode.NO_DISTANCE);
                    if(!stiNode.isLeaf()) {
                        stiNode.setName("");
                    }
                }

                System.out.println(tree);
            }

            List<List<MutableTuple<Tree, Double>>> parentalTrees0 = new ArrayList<>();
            for(Tree tree : parentalTrees) {
                parentalTrees0.add(Arrays.asList(new MutableTuple<Tree, Double>(Trees.readTree(tree.toNewick()), new Double(1.0))));

            }

            try {
                ExecutorService executor = Executors.newFixedThreadPool(1);

                Network<Object> inferredNetwork = null;


                Future<Network<Object>> future = executor.submit(new Callable<Network<Object>>() {
                    @Override
                    public Network<Object> call() {
                        Network<Object> inferredNetwork = inferNetworkFromParentalTrees.inferNetwork(parentalTrees0);
                        return inferredNetwork;
                    }
                });
                executor.shutdown();

                inferredNetwork = future.get(10, TimeUnit.SECONDS);

                System.out.println("Inferred Network: " + inferredNetwork.toString());

                LinkedList<Tree> trees1 = new LinkedList<Tree>();
                LinkedList<Tree> trees2 = new LinkedList<Tree>();


                for (NetworkTree<Object> nt : Networks.getTrees(trueNetwork)) {
                    trees1.add(nt.makeTree());
                }

                for (NetworkTree<Object> nt : Networks.getTrees(inferredNetwork)) {
                    trees2.add(nt.makeTree());
                }

                double dist[] = Networks.computeTreeDistance(trees1, trees2);
                System.out.println("\nThe tree-based distance between two networks: " + dist[0] + " " + dist[1] + " " + dist[2] + "\n");

                if (dist[2] < 1e-6)
                    correctness++;
                else
                    System.out.println("Incorrect!");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.out.println("Correctness " + correctness);

    }
}
