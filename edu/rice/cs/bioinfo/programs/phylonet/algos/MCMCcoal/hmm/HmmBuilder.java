package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.RecombinationRate;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Default1dDict;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Default2dDict;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.util.Pair;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.util.*;
import java.util.stream.DoubleStream;


/**
 * Build coal-HMM from model species tree.
 * Created by Xinhao Liu on 12/2/19.
 */
public class HmmBuilder {

    private STITree<TreeNodeInfo> speciesTree;
    private RecombinationRate recombRate;

    private Map<String, Integer> _speciesName2MSName;
    private Map<Integer, String> _msName2SpeciesName;

    public HmmBuilder(STITree<TreeNodeInfo> tree, RecombinationRate recombRate) {
        this.speciesTree = tree;
        this.recombRate = recombRate;
    }

    /**
     * This function is to generate ms command given ModelTree._tree and ModelTree._recombRate
     */
    private String generateMSCommand() {
//        int N0 = 10000; // move to util? seuqencelength also move to util?
        // take ModelTree._tree and ModelTree._recombRate to generate MS command
        _speciesName2MSName = new HashMap<>();
        _msName2SpeciesName = new HashMap<>();
        List<Tuple3<STINode<TreeNodeInfo>, Double, Double>> nodeTimePopsizeList = new ArrayList<>();
        processTree(speciesTree, _speciesName2MSName, _msName2SpeciesName, nodeTimePopsizeList);

        // basic command
        int nsam = speciesTree.getLeafCount();
        int nreps = 1;
        // spatial structure
        int npop = speciesTree.getLeafCount();

        String outputCommand = "ms";
        outputCommand = outputCommand + " " + nsam + " " + nreps + " -T";

        //TODO: JUST TEST
        Utils.sequenceLength = (int) (300 / (4 * Utils.N0 * recombRate.getRecombRate()));
        System.out.println(Utils.sequenceLength);
        // add recombination
        double rho = 4 * Utils.N0 * recombRate.getRecombRate() * Utils.sequenceLength;
        outputCommand = outputCommand + " -r" + " " + rho + " " + Utils.sequenceLength;

        outputCommand = outputCommand + " -I " + npop;
        for (int i = 0; i < npop; i++) outputCommand += " 1";

        // Specify size of present day populations relative to N0
        for (String leafName:speciesTree.getLeaves()) {
            int popIndex = _speciesName2MSName.get(leafName);
            STINode<TreeNodeInfo> leafNode = speciesTree.getNode(leafName);
            double ratio = leafNode.getData().getPopSize() / (double) Utils.N0;
            outputCommand += " -n " + popIndex + " " + ratio;
        }

        Map<BitSet, Integer> edge2population = new HashMap<>();

        for (Tuple3<STINode<TreeNodeInfo>, Double, Double> bundle:nodeTimePopsizeList) {
            STINode<TreeNodeInfo> node = bundle.Item1;
            int nodeIndex = node.getData().getIndex();    // Index of this node
            double coalescentHeight = bundle.Item2;   // Height of this node in unit of 4N0 generations
            double relativePopsize = bundle.Item3;    // Population size of this node relative to N0 (e.g. 5000 => 0.5 if N0=10000)

            if (node.isLeaf()) {
                STINode<TreeNodeInfo> parent = node.getParent();
                BitSet edge = new BitSet();
                edge.set(parent.getData().getIndex());
                edge.set(nodeIndex);
                edge2population.put(edge, _speciesName2MSName.get(node.getName()));
            } else {
                Iterator<STINode<TreeNodeInfo>> children = node.getChildren().iterator();
                STINode<TreeNodeInfo> child1 = children.next();
                STINode<TreeNodeInfo> child2 = children.next();

                BitSet edge1 = new BitSet();
                edge1.set(child1.getData().getIndex());
                edge1.set(node.getData().getIndex());
                int population1 = edge2population.get(edge1);
                BitSet edge2 = new BitSet();
                edge2.set(node.getData().getIndex());
                edge2.set(child2.getData().getIndex());
                int population2 = edge2population.get(edge2);

                // move lineages and set population size
                outputCommand += " -ej " + coalescentHeight + " " + population1 + " " + population2;
                outputCommand += " -en " + coalescentHeight + " " + population2 + " " + relativePopsize;

                if (node.isRoot()) {
                    continue; // TODO: or break? do postorder traversal in processTree?
                }

                STINode<TreeNodeInfo> parent = node.getParent();
                BitSet edge = new BitSet();
                edge.set(parent.getData().getIndex());
                edge.set(nodeIndex);
                edge2population.put(edge, population2);
            }
        }

        //add seed
//        int seed1 = ThreadLocalRandom.current().nextInt(500, 50000);
//        int seed2 = ThreadLocalRandom.current().nextInt(500, 50000);
//        int seed3 = ThreadLocalRandom.current().nextInt(500, 50000);
        int seed1 = Randomizer.getIntRange(500, 50000);
        int seed2 = Randomizer.getIntRange(500, 50000);
        int seed3 = Randomizer.getIntRange(500, 50000);
        outputCommand += " -seeds " + seed1 + " " + seed2 + " " + seed3;
        return outputCommand;
    }

    /**
     * This function is to gather information of the species tree for generating ms command
     *
     * @param tree, the species tree to generate ms command
     * @param speciesName2MSName, a map from species name to ms population index; initially empty
     * @param msName2SpeciesName, a map from ms population index to species name; initially empty
     * @param nodeTimePopsizeList, should be a post-order traversal of nodes; initially empty
     */
    private void processTree(STITree<TreeNodeInfo> tree, Map<String, Integer> speciesName2MSName, Map<Integer,String> msName2SpeciesName, List<Tuple3<STINode<TreeNodeInfo>, Double, Double>> nodeTimePopsizeList) {
        int nodeIndex = 0;
        int leafIndex = 1;


//        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
//            node.getData().setIndex(nodeIndex++);
//            if (node.isLeaf()) {
//                int populationIndex = leafIndex;
//                leafIndex++;
//                speciesName2MSName.put(node.getName(), populationIndex);
//                msName2SpeciesName.put(populationIndex, node.getName());
//            }
//        }
        // Set index of each node and assign a population number to each leaf
        for (TNode node:tree.postTraverse()) {
            STINode<TreeNodeInfo> stiNode = (STINode<TreeNodeInfo>) node;
            stiNode.getData().setIndex(nodeIndex++);
            if (stiNode.isLeaf()) {
                int populationIndex = leafIndex;
                leafIndex++;
                speciesName2MSName.put(stiNode.getName(), populationIndex);
                msName2SpeciesName.put(populationIndex, stiNode.getName());
            }
        }

        for (TNode node:tree.postTraverse()) {
            STINode<TreeNodeInfo> stiNode = (STINode<TreeNodeInfo>) node;
            double generationHeight = stiNode.getNodeHeight();
            double coalescentHeight = generationHeight / (4 * Utils.N0);
            int popSize = stiNode.getData().getPopSize();
            double relativePopsize = (double) popSize / Utils.N0;

            Tuple3<STINode<TreeNodeInfo>, Double, Double> bundle = new Tuple3<>(stiNode, coalescentHeight, relativePopsize);
            nodeTimePopsizeList.add(bundle);
        }
    }

    private void buildGTNodeHeight(STITree<TreeNodeInfo> gt) {
        for (TNode node:gt.postTraverse()) {
            if (node.isLeaf()) {
                node.setNodeHeight(0);
            } else {
                double height = 0;
                for (TNode child:node.getChildren()) {
                    height = Math.max(height, child.getNodeHeight() + child.getParentDistance());
                }
                node.setNodeHeight(height);
            }
        }
    }

    /**
     * Returns the coalescent history h of a given gene tree as described in Degnan & Salter 2005
     * Returns a STITree along the way.
     * This version of the function works with ms coalescent trees.
     */
    public Tuple<STITree<TreeNodeInfo>, List<Integer>> getCoalescentHistory(String newickString) {
        List<Integer> h = new LinkedList<>();
        STITree<TreeNodeInfo> gt = null;
        try {
            gt = new STITree<>(newickString); //TODO: abused TreeNodeInfo here! should actually be Tree?
            Trees.convertToLexicographicTree(gt);
            buildGTNodeHeight(gt);
            for (TNode gtNode:gt.postTraverse()) {
                // for each internal node of gene tree in post order, excluding root
                if (!gtNode.isLeaf() && !gtNode.isRoot()) {
                    Set<String> clade = new HashSet<>();
                    for (TNode leaf:gtNode.getLeaves()) {
                        clade.add(_msName2SpeciesName.get(Integer.parseInt(leaf.getName())));
                    }
                    for (TNode stNode:speciesTree.postTraverse()) {
                        if (!stNode.isLeaf()) {
                            Set<String> leaves = new HashSet<>();
                            for (TNode leaf:stNode.getLeaves()) {
                                leaves.add(leaf.getName());
                            }
                            // coalescences between two lineages must occur at least as anciently on the species tree
                            // as the most recent common ancestor of the lineages coalescing
                            if (leaves.containsAll(clade)) {
                                STINode<TreeNodeInfo> stiStNode = (STINode<TreeNodeInfo>) stNode;
                                if (stNode.isRoot()) {
                                    h.add(stiStNode.getData().getLabel());
                                } else {
                                    STINode<TreeNodeInfo> p = stiStNode.getParent();
                                    double coalescentThisNodeHeight = stiStNode.getNodeHeight() / (double) (4 * Utils.N0);
                                    double coalescentParentNodeHeight = p.getNodeHeight() / (double) (4 * Utils.N0);
                                    if (coalescentThisNodeHeight <= gtNode.getNodeHeight() && gtNode.getNodeHeight() < coalescentParentNodeHeight) {
                                        h.add(stiStNode.getData().getLabel());
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } catch (IOException|ParseException e) {
            e.printStackTrace();
        }
        return new Tuple<>(gt, h);
    }

    public Tuple3<List<List<Tuple<STITree<TreeNodeInfo>, List<Integer>>>>, Default1dDict, Default2dDict> computeCounts() {
        Runtime rt = Runtime.getRuntime();
        String command = generateMSCommand();

//        System.out.println(command);

        List<List<Tuple<STITree<TreeNodeInfo>, List<Integer>>>> coalescentHistories = new LinkedList<>();
        Default1dDict stateCount = new Default1dDict(0);
        Default2dDict transitionCount = new Default2dDict();
        try {
            Process pr = rt.exec(command);
            //pr.waitFor();

            BufferedReader result = new BufferedReader(new InputStreamReader(pr.getInputStream()));
            String line;
            while (!(line = result.readLine()).equals("//")) {
//                System.out.println(line);
            }

            // begin main logic
            int prevState = -1;
            while ((line = result.readLine()) != null) {
                String[] parts = line.split("]");
                int length = Integer.parseInt(parts[0].substring(1));
                String treeString = parts[1];
//                System.out.println(length);
//                System.out.println(treeString);
                Tuple<STITree<TreeNodeInfo>, List<Integer>> info = getCoalescentHistory(treeString);
                boolean matched = false;
                for (int i = 0; i < coalescentHistories.size(); i++) {
                    List<Tuple<STITree<TreeNodeInfo>, List<Integer>>> pool = coalescentHistories.get(i);
                    STITree<TreeNodeInfo> representativeTree = pool.get(0).Item1;
                    List<Integer> representativeH = pool.get(0).Item2;
                    if (Trees.haveSameRootedTopology(representativeTree, info.Item1) && representativeH.equals(info.Item2)) {
                        Tuple<STITree<TreeNodeInfo>, List<Integer>> newMember = new Tuple<>(info.Item1, info.Item2);
                        pool.add(newMember);
                        stateCount.put(i, stateCount.get(i) + length);
                        if (prevState != -1) {
                            transitionCount.put(prevState, i, transitionCount.get(prevState, i) + 1);
                        }
                        transitionCount.put(i, i, transitionCount.get(i, i) + length - 1);
                        prevState = i;
                        matched = true;
                        break;
                    }
                }
                if (!matched) {
                    List<Tuple<STITree<TreeNodeInfo>, List<Integer>>> newPool = new LinkedList<>();
                    newPool.add(new Tuple<>(info.Item1, info.Item2));
                    coalescentHistories.add(newPool);
                    int i = coalescentHistories.size() - 1;
                    stateCount.put(i, stateCount.get(i) + length);
                    if (prevState != -1) {
                        transitionCount.put(prevState, i, transitionCount.get(prevState, i) + 1);
                    }
                    transitionCount.put(i, i, transitionCount.get(i, i) + length - 1);
                    prevState = i;
                }
            }
            pr.waitFor();
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
        return new Tuple3<>(coalescentHistories, stateCount, transitionCount);
    }


    /**
     * Summarize each coalescent history (average branch lengths).
     */
    public List<HiddenState> summarizeHiddenStates(List<List<Tuple<STITree<TreeNodeInfo>, List<Integer>>>> coalescentHistories) {
        List<HiddenState> hiddenStates = new ArrayList<>();
        for (int stateIndex = 0; stateIndex < coalescentHistories.size(); stateIndex++) {
            List<Tuple<STITree<TreeNodeInfo>, List<Integer>>> pool = coalescentHistories.get(stateIndex);
            STITree<TreeNodeInfo> representativeTree = new STITree<>(pool.get(0).Item1);
            buildGTNodeHeight(representativeTree);
            for (int i = 1; i < pool.size(); i++) {
                STITree<TreeNodeInfo> tree = pool.get(i).Item1;
                Map<TNode, TNode> sameTopoMap = Trees.mapTwoTopologies(representativeTree, tree);
                for (TNode node:representativeTree.postTraverse()) {
                    node.setNodeHeight(node.getNodeHeight() + sameTopoMap.get(node).getNodeHeight());
                }
            }
            for (TNode node:representativeTree.postTraverse()) {
                node.setNodeHeight(node.getNodeHeight() / pool.size());
            }
            for (TNode node:representativeTree.postTraverse()) {
                if (!node.isRoot()) {
                    node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
                }
            }
            // reset gene tree leaf names to species names
            for (String leafName:representativeTree.getLeaves()) {
                TMutableNode leafNode = representativeTree.getNode(leafName);
                leafNode.setName(_msName2SpeciesName.get(Integer.parseInt(leafName)));
            }
            HiddenState state = new HiddenState(representativeTree, stateIndex);
            hiddenStates.add(state);
        }
        return hiddenStates;
    }


    public HmmCore build() {
        Tuple3<List<List<Tuple<STITree<TreeNodeInfo>, List<Integer>>>>, Default1dDict, Default2dDict> metaInfo = computeCounts();
        List<List<Tuple<STITree<TreeNodeInfo>, List<Integer>>>> coalescentHistories = metaInfo.Item1;
        Default1dDict stateCount = metaInfo.Item2;
        Default2dDict transitionCount = metaInfo.Item3;
        int numStates = coalescentHistories.size();

        List<HiddenState> hiddenStates = summarizeHiddenStates(coalescentHistories);
        double[] pi = new double[numStates];
        double[][] a = new double[numStates][numStates];
        for (int i = 0; i < numStates; i++) {
            pi[i] = stateCount.get(i) / (double) Utils.sequenceLength;
        }
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                a[i][j] = (transitionCount.get(i, j) + Utils.hmmSmoothingParam) / (double) (stateCount.get(i) + Utils.hmmSmoothingParam * numStates);
            }
        }
        // normalize
        for (int i = 0; i < numStates; i++) {
            double rowSum = DoubleStream.of(a[i]).sum();
            for (int j = 0; j < numStates; j++) {
                a[i][j] = a[i][j] / rowSum;
            }
        }

        return new HmmCore(pi, a, hiddenStates);
    }

    public static ModelTree getHCGModel() {
        RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(30000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(40000);
                node.setNodeHeight(160000);
            }
        }
        tree.getRoot().getData().setPopSize(40000);
//        tree.getRoot().getData().setPopSize(400000);
        tree.getRoot().setNodeHeight(220000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGModelBad() {
        // Set HC and HCG population size to 60000, change HC divergence time to 100000 gens
        RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(30000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(60000);
                node.setNodeHeight(100000);
            }
        }
        tree.getRoot().getData().setPopSize(60000);
//        tree.getRoot().getData().setPopSize(400000);
        tree.getRoot().setNodeHeight(400000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGModelWrongRecomb(double recombRate) {
        RecombinationRate dummyRecombRate = new RecombinationRate(recombRate);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(30000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(40000);
                node.setNodeHeight(160000);
            }
        }
        tree.getRoot().getData().setPopSize(40000);
//        tree.getRoot().getData().setPopSize(400000);
        tree.getRoot().setNodeHeight(220000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static void main(String[] args) {
//        RecombinationRate dummyRecombRate = new RecombinationRate(1E-8);
        //RecombinationRate dummyRecombRate = new RecombinationRate(1E-6);
//        ModelTree model = new ModelTree("((A,B),(C,D));", dummyRecombRate);
        ModelTree model = getHCGModel();
//        ModelTree model = new ModelTree("(((D,(F,E)),(C,(B,A))),(J,(K,(I,G))));", dummyRecombRate);
        //System.out.println(generateMSCommand(model.getTree(), model.getRecombRate(), 10000));
        HmmBuilder builder = new HmmBuilder(model.getTree(), model.getRecombRate());
        //builder.build();
        //builder.generateMSCommand();
//        builder.getCoalescentHistory("((1:37.422,2:37.422):36.297,(3:37.819,4:37.819):35.899);");
        //List<Integer> h = builder.getCoalescentHistory("((3:9.622,4:9.622):9.764,(1:10.322,2:10.322):9.064);");
        //System.out.println(Arrays.toString(h.toArray()));
//        System.out.println(builder._speciesName2MSName);
//        System.out.println(builder._msName2SpeciesName);


//        Tuple3<List<List<Tuple<STITree<TreeNodeInfo>, List<Integer>>>>, Default1dDict, Default2dDict> result = builder.computeCounts();
//        List<List<Tuple<STITree<TreeNodeInfo>, List<Integer>>>> coalescentHistory = result.Item1;
//        System.out.println(coalescentHistory.size());
//        for (List<Tuple<STITree<TreeNodeInfo>, List<Integer>>> lst:coalescentHistory) {
//            System.out.println(lst.size());
//            for (Tuple<STITree<TreeNodeInfo>, List<Integer>> tup:lst) {
//                System.out.println(tup.Item1.toNewick());
//                System.out.println(Arrays.toString(tup.Item2.toArray()));
//            }
//        }

//        List<HiddenState> hiddenStates = builder.summarizeHiddenStates(coalescentHistory);
//        for (HiddenState state:hiddenStates) {
//            System.out.println(state.getTree().toNewick());
//        }

        HmmCore hmm = builder.build();
        System.out.println(Arrays.toString(hmm.getPi()));
        System.out.println(Arrays.deepToString(hmm.getA()));

        for (HiddenState state:hmm.getStates()) {
            System.out.println("============");
            System.out.println(state.getIndex());
            System.out.println(state.getTree().toNewick());
//            Tree tree = state.getTree();
//            TMutableNode node = (TMutableNode) tree.getNode("1");
//            node.setName("hhhhh");
//            System.out.println(node.getName());
//            System.out.println(state.getTree().toNewick());
//            for (String leafName:tree.getLeaves()) {
//                TMutableNode leafNode = (TMutableNode) tree.getNode(leafName);
//                leafNode.setName(builder._msName2SpeciesName.get(Integer.parseInt(leafName)));
//            }
//            System.out.println(state.getTree().toNewick());
        }
    }
}
