package edu.rice.cs.bioinfo.programs.phylonet.algos.clustering;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
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
import org.apache.commons.math3.util.Pair;

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
                    stack.push(new Pair<>(entry.getKey(), 1));
                }
                leaves.add(entry.getKey() + "_1");
            }
        }

        List<STITree<Double>> parentalTrees = new ArrayList<>();

        int k = stack.size() - 1;
        while(k >= 0) {
            parentalTrees.add(getSubtree(mulTree, leaves));
            //if(parentalTrees.size() > 1000) break;
            leaves.remove(stack.get(k).getKey() + "_" + stack.get(k).getValue());
            stack.set(k, new Pair<>(stack.get(k).getKey(), stack.get(k).getValue() + 1));
            leaves.add(stack.get(k).getKey() + "_" + stack.get(k).getValue());
            while(k > 0 && stack.get(k).getValue() > nname2tamount.get(stack.get(k).getKey())) {
                leaves.remove(stack.get(k).getKey() + "_" + stack.get(k).getValue());
                stack.set(k, new Pair<>(stack.get(k).getKey(), 1));
                leaves.add(stack.get(k).getKey() + "_" + stack.get(k).getValue());
                k--;
                leaves.remove(stack.get(k).getKey() + "_" + stack.get(k).getValue());
                stack.set(k, new Pair<>(stack.get(k).getKey(), stack.get(k).getValue() + 1));
                leaves.add(stack.get(k).getKey() + "_" + stack.get(k).getValue());
            }
            if(k == 0 && stack.get(k).getValue() > nname2tamount.get(stack.get(k).getKey()))
                break;
            k = stack.size() - 1;
        }

        List<Tree> newParentalTrees = new ArrayList<>();
        Set<String> newickStrings = new HashSet<>();

        for(STITree<Double> tree : parentalTrees) {
            for(STINode<Double> node : tree.getNodes()) {
                if(node.isLeaf()) {
                    String name = node.getName();
                    node.setName(name.substring(0, name.lastIndexOf("_")));
                } else {
                    node.setName("");
                }
                node.setParentDistance(TNode.NO_DISTANCE);
            }
            Trees.removeBinaryNodes(tree);
            Trees.convertToLexicographicTree(tree);
            if(newickStrings.contains(tree.toNewick()))
                continue;
            newickStrings.add(tree.toNewick());
            newParentalTrees.add(Trees.readTree(tree.toNewick()));
        }


        return newParentalTrees;
    }

    public void test() {
        int number = 100;
        int correctness = 0;
        double remainRatio = 1;


        List<Network<Object>> incorrectNetworks = new ArrayList<>();
        List<Network<Object>> exceptionNetworks = new ArrayList<>();
        StringBuilder result = new StringBuilder();

        for(int i = 0 ; i < number ; i++) {
            DataGenerator dataGenerator = new DataGenerator();
            ParentalTreeOperation parentalTreeOperation = new ParentalTreeOperation();
            Network<Object> trueNetwork = dataGenerator.getRandomNetwork(10);

            /*double gamma = 1 - 0.05;

            double y = 0.1;
            double x = 10000;

            trueNetwork = Networks.readNetwork("((((b:1000.0,c:1000.0)I4:" + y + ")I3#H1:0.0::" + (1 - gamma) + ",a:" + (1000 + y) + ")I1:" + x + ",(I3#H1:0.0::" + gamma + ",d:" + (1000 + y) + ")I2:"+ x +")I0;");
*/
            //trueNetwork = Networks.readNetwork("(((B)I1#H1,A),(((C)I2#H2,D),(I1#H1,I2#H2)));");
            //List<Tree> pt = getParentalTrees(trueNetwork);
            //trueNetwork = Networks.readNetwork("((6,(4,(10,7)I6)I4)I2,(2,(((8)I8#H1:::0.5,(((9)I12#H2:::0.5,1)I10,((I12#H2:::0.5,5)I11,I8#H1:::0.5)I9)I7)I5,3)I3)I1)I0;");
            //trueNetwork = Networks.readNetwork("(((5,6)I5,((((((((10,8)I10,4)I8,1)I6)I3#H2:::0.5,2)I9)I12#H1:::0.5,3)I11,(I12#H1:::0.5,7)I7)I4)I2,(I3#H2:::0.5,9)I1)I0;");
            //trueNetwork = Networks.readNetwork("(((3,1)I5,((((5,8)I10,9)I9)I8#H1,10)I4)I2,(2,(((6)I12#H2,I8#H1)I7,((I12#H2,4)I11,7)I6)I3)I1)I0;");
            //trueNetwork = Networks.readNetwork("((10,7)I1,(((((2)I11#H1:::0.5,(8)I6#H2:::0.5)I12,9)I4,(((1,6)I10,3)I9,(I6#H2:::0.5,5)I8)I7)I5,(I11#H1:::0.5,4)I3)I2)I0;");
            //trueNetwork = Networks.readNetwork("((10,4)I1,(7,(((9)I8#H1:::0.5,((8,(1)I11#H2:::0.5)I9,(I8#H1:::0.5,(2,I11#H2:::0.5)I12)I10)I7)I5,((3,6)I6,5)I4)I3)I2)I0;");
            //trueNetwork = Networks.readNetwork("((((1,(4,2)I5)I3,((6)I11#H1:::0.5,(((9)I9#H2:::0.5,((I11#H1:::0.5,(I9#H2:::0.5,3)I10)I12,7)I8)I7,10)I6)I4)I2,8)I1,5)I0;");
            //trueNetwork = Networks.readNetwork("((((7,6)I10)I6#H1:::0.5,((((5)I11#H2:::0.5,1)I12,I6#H1:::0.5)I9,(I11#H2:::0.5,10)I8)I5)I2,((8,2)I4,((9,3)I7,4)I3)I1)I0;");
            //trueNetwork = Networks.readNetwork("(((10,(((((((1)I5#H2,2)I10,7)I9,4)I8,6)I7,3)I6,9)I4)I2)I11#H1,((I11#H1,(I5#H2,5)I3)I12,8)I1)I0;");
            //trueNetwork = Networks.readNetwork("(((((((9)I8#H2:::0.5,3)I9,2)I5,(I8#H2:::0.5,(1,6)I7)I4)I2)I12#H1:::0.5,(I12#H1:::0.5,(4,((10,7)I10,8)I6)I3)I11)I1,5)I0;");
            Networks.removeAllParameters(trueNetwork);
            //System.out.println(trueNetwork.toString());
            //System.out.println("Level-" + Networks.computeLevel(trueNetwork));
            //System.out.println("IsGalledNetwork: " + Networks.isGalledNetwork(trueNetwork));
            //System.out.println("Num of PTs: " + pt.size());

            //trueNetwork = Networks.readNetwork("(((((3,((8,10)I10,7)I9)I7)I5#H1:::0.5,(((I5#H1:::0.5,1)I8,6)I6,4)I4)I3,(5,2)I2)I1,9)I0;");
            //trueNetwork = Networks.readNetwork("((((((6,3)I6)I4#H1:::0.5,((I4#H1:::0.5,(7,8)I9)I8,(5,1)I10)I7)I5,(4,2)I3)I2,9)I1,10)I0;");
            //trueNetwork = Networks.readNetwork("(((((2,1)I9,10)I7,3)I4)I2#H1:::0.5,(((I2#H1:::0.5,(8,(6,7)I10)I8)I6,(4,5)I5)I3,9)I1)I0;"); //TODO: From root edge
            //trueNetwork = Networks.readNetwork("(((((2,(3,5)I9)I8)I10#H1:::0.5,8)I7,(I10#H1:::0.5,7)I4)I2,(((10,4)I6,(9,1)I5)I3,6)I1)I0;");
            //trueNetwork = Networks.readNetwork("(((5,((2,10)I7,((3,(1,6)I10)I9,(7,8)I8)I6)I5)I4)I2#H1,((I2#H1,9)I3,4)I1)I0;");
            //trueNetwork = Networks.readNetwork("(((((8,3)I10,(5,10)I8)I7,(1,(4,9)I6)I4)I2)I5#H1:::0.5,((I5#H1:::0.5,6)I9,(7,2)I3)I1)I0;");
            //trueNetwork = Networks.readNetwork("(((6:1.0,7:1.0)I5:1.0,(8:1.0,(((4:1.0,2:1.0)I10:1.0)I8#H1:0.0,((I8#H1:0.0,9:1.0)I9:1.0,1:1.0)I7:1.0)I6:1.0)I4:1.0)I2:1.0,((5:1.0,3:1.0)I3:1.0,10:1.0)I1:1.0)I0;");
            //trueNetwork = Networks.readNetwork("((((((10,4)I8)I6#H1,7)I4,(I6#H1,((9,(5,(8,6)I10)I9)I7,2)I5)I3)I2,1)I1,3)I0;");
            //trueNetwork = Networks.readNetwork("((8,(((((4,10)I9)I10#H1:::0.5,9)I4,(I10#H1:::0.5,1)I8)I7,((2,7)I6,(6,5)I5)I3)I2)I1,3)I0;");
            //trueNetwork = Networks.readNetwork("(((((((4)I12#H3:::0.5)I17#H2:::0.5,(1,10)I9)I18)I16#H1:::0.5,((I16#H1:::0.5,I12#H3:::0.5)I11,2)I5)I15,(I17#H2:::0.5,((5)I4#H5:::0.5)I13#H4:::0.5)I6)I2,(I4#H5:::0.5,((9,7)I8,(((I13#H4:::0.5,8)I14,3)I10,6)I7)I3)I1)I0;");
            for(NetNode<Object> node: trueNetwork.dfs()){
                for(NetNode<Object> parent: node.getParents()){
                    node.setParentDistance(parent,NetNode.NO_DISTANCE);
                    if(node.isNetworkNode()){
                        if(node.getParentProbability(parent) == NetNode.NO_PROBABILITY)
                            node.setParentProbability(parent, 0.5);
                    }
                }
            }

            if(! (Networks.computeLevel(trueNetwork) == 1))
                continue;

            System.out.println("True network: " + trueNetwork.toString());

            List<Tree> parentalTrees = parentalTreeOperation.getParentalTrees(trueNetwork);

            if(true) {
                System.out.println("Number of parental trees: " + parentalTrees.size());
                continue;
            }
            //if(parentalTrees.size() > 128) {correctness++;continue;}

            Collections.shuffle(parentalTrees);
            parentalTrees = parentalTrees.subList(0, (int)(parentalTrees.size() * remainRatio));
            //parentalTrees.clear();
            //parentalTrees.add(Trees.readTree("((((((1,8),2),9),((3,6),5)),7),(10,4));"));
            //parentalTrees.add(Trees.readTree("(((((1,8),(2,9)),((3,6),5)),7),(10,4));"));
            //parentalTrees.add(Trees.readTree("((((((1,2),8),9),((3,6),5)),7),(10,4));"));
            //parentalTrees.add(Trees.readTree("((((((1,2),9),8),((3,6),5)),7),(10,4));"));

            System.out.println("Number of input parental trees: " + parentalTrees.size());
            //if(!Networks.isGalledNetwork(trueNetwork) /*|| !Networks.isNestedNetwork(trueNetwork)*/) {
            //    correctness++;
            //    continue;
            //}

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
                        InferNetworkFromParentalTrees inferNetworkFromParentalTrees = new InferNetworkFromParentalTrees();
                        Network<Object> inferredNetwork = inferNetworkFromParentalTrees.inferNetworkHeuristc(parentalTrees0);
                        return inferredNetwork;
                    }
                });
                executor.shutdown();

                inferredNetwork = future.get(1000, TimeUnit.SECONDS);

                executor.awaitTermination(1000, TimeUnit.SECONDS);

                System.out.println("Inferred Network: " + inferredNetwork.toString());

                List<Tree> newParentalTrees = parentalTreeOperation.getParentalTrees(inferredNetwork);
                Set<String> newParentalTreesNewick = new HashSet<>();
                Set<String> originalParentalTreesNewick = new HashSet<>();
                for(Tree t : newParentalTrees) {
                    newParentalTreesNewick.add(t.toNewick());
                }
                for(Tree t : parentalTrees) {
                    originalParentalTreesNewick.add(t.toNewick());
                }

                LinkedList<Tree> trees1 = new LinkedList<Tree>();
                LinkedList<Tree> trees2 = new LinkedList<Tree>();


                for (NetworkTree<Object> nt : Networks.getTrees(trueNetwork)) {
                    trees1.add(nt.makeTree());
                }

                for (NetworkTree<Object> nt : Networks.getTrees(inferredNetwork)) {
                    trees2.add(nt.makeTree());
                }

                double dist[] = Networks.computeTreeDistance(trees1, trees2);
                //System.out.println("\nThe tree-based distance between two networks: " + dist[0] + " " + dist[1] + " " + dist[2] + "\n");

                if(newParentalTreesNewick.containsAll(originalParentalTreesNewick))
                //if (dist[2] < 1e-6)
                    correctness++;
                else {
                    System.out.println("Incorrect!");
                    incorrectNetworks.add(trueNetwork);
                    for(String s : newParentalTreesNewick) {
                        System.out.println("New parental tree: " + s);
                    }
                    for(String s : originalParentalTreesNewick) {
                        if(!newParentalTreesNewick.contains(s))
                        System.out.println("Unshown parental tree: " + s);
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
                exceptionNetworks.add(trueNetwork);
            }
        }
        System.out.println("Correctness " + correctness);
        System.out.println("IncorrectNetworks:");
        for(Network net : incorrectNetworks)  {
            System.out.println(net.toString());
        }
        System.out.println("ExceptionNetworks:");
        for(Network net : exceptionNetworks)  {
            System.out.println(net.toString());
        }
    }

}
