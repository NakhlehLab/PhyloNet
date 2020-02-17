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
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.DoubleStream;
import org.apache.commons.io.FileUtils;


/**
 * Build coal-HMM from model species tree using simulator fastsimcoal2.
 * Created by Xinhao Liu on 01/16/20.
 */
public class HmmBuilderFSC {

    private STITree<TreeNodeInfo> speciesTree;
    private RecombinationRate recombRate;

    private Map<String, Integer> _speciesName2DemeName;
    private Map<Integer, String> _demeName2SpeciesName;

    private String simulationName = "simulation";

    public HmmBuilderFSC(STITree<TreeNodeInfo> tree, RecombinationRate recombRate) {
        this.speciesTree = tree;
        this.recombRate = recombRate;
    }

    /**
     * This function is to generate fsc parameter file given ModelTree._tree and ModelTree._recombRate
     */
    private void generateFSCParFile(String fileName) {
        List<String> lines = new ArrayList<>();

        // TODO: take ModelTree._tree and ModelTree._recombRate to generate fsc command
        _speciesName2DemeName = new HashMap<>();
        _demeName2SpeciesName = new HashMap<>();
        List<Tuple3<STINode<TreeNodeInfo>, Double, Double>> nodeTimePopsizeList = new ArrayList<>();
        processTree(speciesTree, _speciesName2DemeName, _demeName2SpeciesName, nodeTimePopsizeList);

        // Number of population samples (demes)
        int npop = speciesTree.getLeafCount();
        lines.add("//Number of population samples (demes)");
        lines.add(String.valueOf(npop));

        // Population effective sizes (number of genes)
        lines.add("//Population effective sizes (number of genes)");
        for (int i = 0; i < npop; i++) {
            lines.add(String.valueOf(speciesTree.getNode(_demeName2SpeciesName.get(i)).getData().getPopSize()));
        }

        // Sample sizes
        lines.add("//Sample sizes");
        for (int i = 0; i < npop; i++) {
            lines.add("1");
        }

        // Growth rates: negative growth implies population expansion
        lines.add("//Growth rates: negative growth implies population expansion");
        for (int i = 0; i < npop; i++) {
            lines.add("0");
        }

        // Number of migration matrices : 0 implies no migration between demes
        lines.add("//Number of migration matrices: 0 implies no migration between demes");
        lines.add("0");

        // historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix
        lines.add("//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix ");
        lines.add((speciesTree.getNodeCount() - speciesTree.getLeafCount()) + " historical event");

        Map<BitSet, Tuple<Integer, Double>> edge2DemeAndSize = new HashMap<>();

        for (Tuple3<STINode<TreeNodeInfo>, Double, Double> bundle:nodeTimePopsizeList) {
            STINode<TreeNodeInfo> node = bundle.Item1;
            int nodeIndex = node.getData().getIndex();    // Index of this node
            double height = bundle.Item2;   // Height of this node in unit of generations
            double popSize = bundle.Item3;    // Population size of this node in unit of individuals
            //assert

            if (node.isLeaf()) {
                STINode<TreeNodeInfo> parent = node.getParent();
                BitSet edge = new BitSet();
                edge.set(parent.getData().getIndex());
                edge.set(nodeIndex);
                edge2DemeAndSize.put(edge, new Tuple<>(_speciesName2DemeName.get(node.getName()), popSize));
            } else {
                Iterator<STINode<TreeNodeInfo>> children = node.getChildren().iterator();
                STINode<TreeNodeInfo> child1 = children.next();
                STINode<TreeNodeInfo> child2 = children.next();

                BitSet edge1 = new BitSet();
                edge1.set(child1.getData().getIndex());
                edge1.set(node.getData().getIndex());
                int deme1 = edge2DemeAndSize.get(edge1).Item1;
                BitSet edge2 = new BitSet();
                edge2.set(node.getData().getIndex());
                edge2.set(child2.getData().getIndex());
                int deme2 = edge2DemeAndSize.get(edge2).Item1;

                // time, source, sink, migrants
                String event = (int) height + " " + deme1 + " " + deme2 + " " + "1" + " ";

                double relativeNewSize = popSize /  edge2DemeAndSize.get(edge2).Item2;
                // new size, new growth rate, migr. matrix
                event += relativeNewSize + " " + "0" + " " + "0";
                lines.add(event);

                if (node.isRoot()) {
                    continue;
                }

                STINode<TreeNodeInfo> parent = node.getParent();
                BitSet edge = new BitSet();
                edge.set(parent.getData().getIndex());
                edge.set(nodeIndex);
                edge2DemeAndSize.put(edge, new Tuple<>(deme2, popSize));
            }
        }

        //Number of independent loci [chromosome]
        lines.add("//Number of independent loci [chromosome]");
        lines.add("1 0");

        //Per chromosome: Number of linkage blocks
        lines.add("//Per chromosome: Number of linkage blocks");
        lines.add("1");

        //per Block: data type, num loci, rec. rate and mut rate + optional parameters
        lines.add("//per Block: data type, num loci, rec. rate and mut rate + optional parameters");
        String dataType = "DNA" + " " + Utils.sequenceLength + " " + recombRate.toString() + " 0 0.33";
        lines.add(dataType);
        lines.add("");

        Path file = Paths.get(fileName);
        try {
            Files.write(file, lines, StandardCharsets.UTF_8);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * This function is to gather information of the species tree for generating fsc parameter file
     *
     * @param tree, the species tree to generate fsc command
     * @param speciesName2DemeName, a map from species name to fsc deme index; initially empty
     * @param demeName2SpeciesName, a map from fsc deme index to species name; initially empty
     * @param nodeTimePopsizeList, should be a post-order traversal of nodes; initially empty
     */
    private void processTree(STITree<TreeNodeInfo> tree, Map<String, Integer> speciesName2DemeName, Map<Integer,String> demeName2SpeciesName, List<Tuple3<STINode<TreeNodeInfo>, Double, Double>> nodeTimePopsizeList) {
        int nodeIndex = 0;
        int leafIndex = 0;

        // Set index of each node and assign a deme index to each leaf
        for (TNode node:tree.postTraverse()) {
            STINode<TreeNodeInfo> stiNode = (STINode<TreeNodeInfo>) node;
            stiNode.getData().setIndex(nodeIndex++);
            if (stiNode.isLeaf()) {
                int demeIndex = leafIndex;
                leafIndex++;
                speciesName2DemeName.put(stiNode.getName(), demeIndex);
                demeName2SpeciesName.put(demeIndex, stiNode.getName());
            }
        }

        for (TNode node:tree.postTraverse()) {
            STINode<TreeNodeInfo> stiNode = (STINode<TreeNodeInfo>) node;
            double height = stiNode.getNodeHeight();
            double popSize = stiNode.getData().getPopSize();

            Tuple3<STINode<TreeNodeInfo>, Double, Double> bundle = new Tuple3<>(stiNode, height, popSize);
            nodeTimePopsizeList.add(bundle);
        }
    }

    private String generateFSCCommandLine() {
        String fileName = simulationName + ".par";
        generateFSCParFile(fileName);
//        int seed = Randomizer.getRandomInt(100000);
//        return "./fsc26 -i " + fileName + " -n 1 -T -r " + seed;
        return "./fsc26 -i " + fileName + " -n 1 -T";
    }

    private void buildGTNodeHeight(Tree gt) {
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
     * Returns a Tree representation of the gene tree along the way.
     * This version of the function works with fsc26 coalescent trees.
     */
    public Tuple<Tree, List<Integer>> getCoalescentHistory(String newickString) {
        List<Integer> h = new LinkedList<>();
        STITree<TreeNodeInfo> gt = null;
        try {
            gt = new STITree<>(newickString); //TODO: abused TreeNodeInfo here!
            Trees.convertToLexicographicTree(gt);
            buildGTNodeHeight(gt);
            for (TNode gtNode:gt.postTraverse()) {
                // for each internal node of gene tree in post order, excluding root
                if (!gtNode.isLeaf() && !gtNode.isRoot()) {
                    Set<String> clade = new HashSet<>();
                    for (TNode leaf:gtNode.getLeaves()) {
                        int demeName = Integer.parseInt(leaf.getName().split("\\.")[1]) - 1;
                        clade.add(_demeName2SpeciesName.get(demeName));
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
                                    double thisNodeHeight = stiStNode.getNodeHeight();
                                    double parentNodeHeight = p.getNodeHeight();
                                    if (thisNodeHeight <= gtNode.getNodeHeight() && gtNode.getNodeHeight() < parentNodeHeight) {
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

    public Tuple3<List<List<Tuple<Tree, List<Integer>>>>, Default1dDict, Default2dDict> computeCounts() {
        Runtime rt = Runtime.getRuntime();

        String command = generateFSCCommandLine();
        //System.out.println(command);
        List<List<Tuple<Tree, List<Integer>>>> coalescentHistories = new LinkedList<>();
        Default1dDict stateCount = new Default1dDict(0);
        Default2dDict transitionCount = new Default2dDict();
        try {
            Process pr = rt.exec(command);
            // deal with standard output here
            BufferedReader reader = new BufferedReader(new InputStreamReader(pr.getInputStream()));
            String stdout;
            while ((stdout = reader.readLine()) != null) {
            }
            pr.waitFor();

            File treeFile = new File(simulationName + "/" + simulationName + "_1_true_trees.trees");
            BufferedReader br = new BufferedReader(new FileReader(treeFile));
            String line;
            while (!(line = br.readLine()).equals("")) {
//                System.out.println(line);
            }

            List<Tuple<Integer, String>> pairIdxTreeLst = new ArrayList<>();
            while (!(line = br.readLine()).equals("end;")) {
                String[] parts = line.split(" = \\[&U\\] ");
//                System.out.println(Arrays.toString(parts));
                //String startIndex = parts[0].substring(parts[0].lastIndexOf("_") + 1);
                int startIndex = Integer.parseInt(parts[0].substring(parts[0].lastIndexOf("_") + 1));
//                System.out.println(startIndex);
//                System.out.println(parts[1]);
                pairIdxTreeLst.add(new Tuple<>(startIndex, parts[1]));
            }

            // begin main logic
            int prevState = -1;
            for (int sectionIdx = 0; sectionIdx < pairIdxTreeLst.size(); sectionIdx++) {
                int sectionStart = pairIdxTreeLst.get(sectionIdx).Item1;
                int sectionEnd = sectionIdx == pairIdxTreeLst.size() - 1 ? Utils.sequenceLength : pairIdxTreeLst.get(sectionIdx + 1).Item1;
                int length = sectionEnd - sectionStart;
                String treeString = pairIdxTreeLst.get(sectionIdx).Item2;

                Tuple<Tree, List<Integer>> info = getCoalescentHistory(treeString);
                boolean matched = false;
                for (int i = 0; i < coalescentHistories.size(); i++) {
                    List<Tuple<Tree, List<Integer>>> pool = coalescentHistories.get(i);
                    Tree representativeTree = pool.get(0).Item1;
                    List<Integer> representativeH = pool.get(0).Item2;
                    if (Trees.haveSameRootedTopology(representativeTree, info.Item1) && representativeH.equals(info.Item2)) {
                        Tuple<Tree, List<Integer>> newMember = new Tuple<>(info.Item1, info.Item2);
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
                    List<Tuple<Tree, List<Integer>>> newPool = new LinkedList<>();
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
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
        return new Tuple3<>(coalescentHistories, stateCount, transitionCount);
    }


    /**
     * Summarize each coalescent history (average branch lengths).
     */
    public List<HiddenState> summarizeHiddenStates(List<List<Tuple<Tree, List<Integer>>>> coalescentHistories) {
        List<HiddenState> hiddenStates = new ArrayList<>();
        for (int stateIndex = 0; stateIndex < coalescentHistories.size(); stateIndex++) {
            List<Tuple<Tree, List<Integer>>> pool = coalescentHistories.get(stateIndex);
            STITree<TreeNodeInfo> representativeTree = new STITree<>(pool.get(0).Item1);
            buildGTNodeHeight(representativeTree);
            for (int i = 1; i < pool.size(); i++) {
                Tree tree = pool.get(i).Item1;
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
            HiddenState state = new HiddenState(representativeTree, stateIndex);
            hiddenStates.add(state);
        }
        return hiddenStates;
    }


    public HmmCore build() {
        Tuple3<List<List<Tuple<Tree, List<Integer>>>>, Default1dDict, Default2dDict> metaInfo = computeCounts();
        List<List<Tuple<Tree, List<Integer>>>> coalescentHistories = metaInfo.Item1;
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

        // delete fsc generated files
        try {
            FileUtils.deleteDirectory(new File(simulationName));
            File par = new File(simulationName + ".par");
            par.delete();
            File seed = new File("seed.txt");
            seed.delete();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return new HmmCore(pi, a, hiddenStates);
    }


    private static ModelTree getHCGModel() {
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
        tree.getRoot().setNodeHeight(220000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static void main(String[] args) {
        RecombinationRate dummyRecombRate = new RecombinationRate(1E-8);
        //RecombinationRate dummyRecombRate = new RecombinationRate(1E-6);
//        ModelTree model = new ModelTree("((A,B),(C,D));", dummyRecombRate);
//        ModelTree model = new ModelTree("(((D,(F,E)),(C,(B,A))),(J,(K,(I,G))));", dummyRecombRate);
        ModelTree model = getHCGModel();
        //System.out.println(generateMSCommand(model.getTree(), model.getRecombRate(), 10000));
        HmmBuilderFSC builder = new HmmBuilderFSC(model.getTree(), model.getRecombRate());
        //builder.build();
        //builder.generateMSCommand();
//        System.out.println(builder._speciesName2MSName);
//        System.out.println(builder._msName2SpeciesName);


//        Tuple3<List<List<Tuple<Tree, List<Integer>>>>, Default1dDict, Default2dDict> result = builder.computeCounts();
//        List<List<Tuple<Tree, List<Integer>>>> coalescentHistory = result.Item1;
//        System.out.println(coalescentHistory.size());
//        for (List<Tuple<Tree, List<Integer>>> lst:coalescentHistory) {
//            System.out.println(lst.size());
//            for (Tuple<Tree, List<Integer>> tup:lst) {
//                System.out.println(tup.Item1.toNewick());
//                System.out.println(Arrays.toString(tup.Item2.toArray()));
//            }
//        }
//
//        List<HiddenState> hiddenStates = builder.summarizeHiddenStates(coalescentHistory);
//        for (HiddenState state:hiddenStates) {
//            System.out.println(state.getTree().toNewick());
//        }

//        HmmCore hmm = builder.build();
//        System.out.println(Arrays.toString(hmm.getPi()));
//        System.out.println(Arrays.deepToString(hmm.getA()));

        //TODO: new tests start here
//        ModelTree model = getHCGModel();
//        HmmBuilderFSC builder = new HmmBuilderFSC(model.getTree(), model.getRecombRate());
//        builder.generateFSCParFile("test.par");
//        System.out.println(builder._speciesName2DemeName);
//        System.out.println(builder._demeName2SpeciesName);
//        System.out.println(System.getProperty("user.dir"));

//        String command = builder.generateFSCCommandLine();
//        Runtime rt = Runtime.getRuntime();
//        try {
//            Process pr = rt.exec(command);
//            pr.waitFor();
//
//            File treeFile = new File(builder.simulationName + "/" + builder.simulationName + "_1_true_trees.trees");
//            BufferedReader br = new BufferedReader(new FileReader(treeFile));
//            String line;
//            while ((line = br.readLine()) != null) {
//                System.out.println(line);
//            }
//
//            FileUtils.deleteDirectory(new File(builder.simulationName));
//            File par = new File(builder.simulationName + ".par");
//            par.delete();
//            File seed = new File("seed.txt");
//            seed.delete();
//        } catch (IOException | InterruptedException e) {
//            e.printStackTrace();
//        }
        HmmCore hmm = builder.build();
        System.out.println(Arrays.toString(hmm.getPi()));
        System.out.println(Arrays.deepToString(hmm.getA()));

        for (HiddenState state:hmm.getStates()) {
            System.out.println("============");
            System.out.println(state.getIndex());
            System.out.println(state.getTree().toNewick());
        }


    }
}
