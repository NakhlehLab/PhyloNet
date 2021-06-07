package edu.rice.cs.bioinfo.programs.phylonet.algos.clustering;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.mast.AmirKeselmanMAST;
import edu.rice.cs.bioinfo.programs.phylonet.algos.mast.SteelWarnowMAST;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.HillClimberBase;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.SimpleHillClimbing;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing.SimulatedAnnealingSalterPearL;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomTopologyNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NonUltrametricNetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.File;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 7/4/16
 * Time: 4:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferNetworkFromParentalTrees {
    protected long _maxExaminations = 1000000;
    protected int _maxFailure = 100;
    protected int _numRuns = 10;
    protected int _numThreads = 1;
    protected int _moveDiameter = -1;
    protected int _reticulationDiameter = -1;
    protected double[] _topologyOperationWeights = {0.1,0.1,0.15,0.55,0.15,0.15};
    protected double[] _topologyVsParameterOperation = {1,0};
    protected Long _seed = new Long(12345);
    protected Set<String> _fixedHybrid = new HashSet<String>();
    protected File _logFile = null;
    protected File _intermediateResultFile = null;
    protected Set<String> _inputParentalTreesNewick = null;



    private boolean _verbose = true;

    public void removeBinaryNodes(Network<Object> net)
    {
        // Find all binary nodes.
        List<NetNode<Object>> binaryNodes = new LinkedList<NetNode<Object>>();
        for (NetNode<Object> node : net.bfs()) {
            if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                binaryNodes.add(node);
            }
        }

        // Remove them.
        for (NetNode<Object> node : binaryNodes) {
            NetNode<Object> child = node.getChildren().iterator().next();	// Node's only child.
            if(child.getIndeg() != 1){
                continue;
            }
            NetNode<Object> parent = node.getParents().iterator().next();	// Node's only parent.
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

    public static Map<String, int[]> countTripleFromSubtree(Tree gt){
        List<String> taxaList = new ArrayList<>();
        Map<String,Integer> pairwiseDepths = new HashMap<>();
        Map<TNode, Integer> node2depth = new HashMap<>();
        Map<TNode, List<String>> node2leaves = new HashMap<>();
        for(TNode node: gt.postTraverse()){
            int depth = 0;
            List<String> leaves = new ArrayList<>();
            if(node.isLeaf()){
                leaves.add(node.getName());
                taxaList.add(node.getName());
            }
            else{
                List<List<String>> childLeavesList = new ArrayList<>();
                for(TNode child: node.getChildren()){
                    depth = Math.max(depth, node2depth.get(child));
                    List<String> childLeaves = new ArrayList<>();
                    childLeaves.addAll(node2leaves.get(child));
                    childLeavesList.add(childLeaves);
                    leaves.addAll(childLeaves);
                }
                depth++;
                for(int i=0; i<childLeavesList.size(); i++){
                    List<String> childLeaves1 = childLeavesList.get(i);
                    for(int j=i+1; j<childLeavesList.size(); j++){
                        List<String> childLeaves2 = childLeavesList.get(j);
                        for(String leaf1: childLeaves1){
                            for(String leaf2: childLeaves2){
                                pairwiseDepths.put(leaf1+"&"+leaf2,depth);
                                pairwiseDepths.put(leaf2+"&"+leaf1,depth);
                            }
                        }
                    }

                }

            }
            node2depth.put(node, depth);
            node2leaves.put(node, leaves);
        }

        String[] taxa = taxaList.toArray(new String[0]);
        Arrays.sort(taxa);
        int numTaxa = taxa.length;
        Map<String, int[]> triple2counts = new HashMap<>();

        for(int i=0; i<numTaxa; i++){
            String taxonI = taxa[i];
            for(int j=i+1; j<numTaxa; j++){
                String taxonJ = taxa[j];
                String pair = taxonI+"&"+taxonJ;
                int ij = pairwiseDepths.get(pair);
                for(int k=j+1; k<numTaxa; k++){
                    String taxonK = taxa[k];
                    String triplet = pair+"&"+taxonK;
                    int[] freq = new int[3];
                    triple2counts.put(triplet, freq);
                    int minIndex = -1;
                    int ik = pairwiseDepths.get(taxonI+"&"+taxonK);
                    int jk = pairwiseDepths.get(taxonJ+"&"+taxonK);
                    // a&b&c
                    if(ij<ik && ij<jk){ //((a,b),c)
                        minIndex = 0;
                    }else if(ik<ij && ik<jk) //((a,c),b)
                    {
                        minIndex = 1;
                    }
                    else if(jk<ij && jk<ik){ //((b,c),a)
                        minIndex = 2;
                    }
                    if(minIndex!=-1){
                        freq[minIndex] = 1;
                    }
                }
            }
        }

        return triple2counts;
    }

    public int getAllTriplets(Tree tree, List<String> leaves, Map<String, int[]> triple2counts) {
        Map<TNode, Integer> marks = new HashMap<>();

        int OUT = 1;
        int IN = 2;

        for(TNode node : tree.postTraverse()) {
            marks.put(node, IN);
        }

        marks.put(tree.getRoot(), OUT);

        for(String leafName : tree.getLeaves()) {
            if(!leaves.contains(leafName)) {
                TNode node = tree.getNode(leafName);
                while(!node.isRoot()) {
                    marks.put(node, OUT);
                    node = node.getParent();
                }
            }
        }

        int count = 0;
        for(TNode node : tree.postTraverse()) {
            if(!node.isRoot()) {
                if(!marks.get(node).equals(marks.get(node.getParent()))) {
                    count++;
                    Map<String, int[]> curTripletCount = countTripleFromSubtree(Trees.readTree(node.toString()));
                    for(String triplet : curTripletCount.keySet()) {
                        if(!triple2counts.containsKey(triplet))
                            triple2counts.put(triplet, curTripletCount.get(triplet));
                        else {
                            int[] sum = triple2counts.get(triplet);
                            int[] cur = curTripletCount.get(triplet);
                            for(int m = 0 ; m < 3 ; m++)
                                sum[m] += cur[m];
                        }
                    }
                }
            }
        }

        return count;
    }

    public int getSubnetworkRetain(Network<Object> net, List<String> leaves) {

        Map<NetNode<Object>, Integer> marks = new HashMap<>();

        int OUT = 1;
        int IN = 2;

        for(NetNode<Object> node : net.bfs()) {
            marks.put(node, IN);
        }

        marks.put(net.getRoot(), OUT);

        for(NetNode<Object> leaf : net.getLeaves()) {
            if(!leaves.contains(leaf.getName())) {
                NetNode<Object> node = leaf;
                while(!node.isRoot()) {
                    marks.put(node, OUT);
                    for(NetNode<Object> parent : node.getParents()) {
                        node = parent;
                        break;
                    }
                }
            }
        }

        int count = 0;
        for(NetNode<Object> node : net.bfs()) {
            if(!node.isRoot()) {
                for (NetNode<Object> parent : node.getParents()) {
                    if(!marks.get(node).equals(marks.get(parent))) {
                        count++;
                    }
                }
            }
        }

        boolean done = false;
        while(!done) {
            done = true;
            for (NetNode<Object> node : net.dfs()) {
                if (!leaves.contains(node.getName()) && node.getChildCount() == 0) {
                    for (NetNode<Object> parent : node.getParents()) {
                        parent.removeChild(node);
                        done = false;
                    }
                }
            }
        }

        return count;
    }

    public NetNode<Object> getMRCA(Network<Object> net, List<String> leaves) {

        Map<NetNode<Object>, Integer> marks = new HashMap<>();

        int OUT = 1;
        int IN = 2;

        for(NetNode<Object> node : net.bfs()) {
            marks.put(node, OUT);
        }

        marks.put(net.getRoot(), IN);

        for(NetNode<Object> leaf : net.getLeaves()) {
            if(leaves.contains(leaf.getName())) {
                NetNode<Object> node = leaf;
                while(!node.isRoot()) {
                    marks.put(node, IN);
                    for(NetNode<Object> parent : node.getParents()) {
                        node = parent;
                        break;
                    }
                }
            }
        }

        NetNode<Object> mrcaRoot = net.getRoot();
        while(!mrcaRoot.isLeaf()) {
            int count = 0;
            NetNode<Object> inChild = null;
            for(NetNode<Object> child : mrcaRoot.getChildren()) {
                if(marks.get(child) == IN) {
                    count++;
                    inChild = child;
                }
            }
            if(count == 1)
                mrcaRoot = inChild;
            else
                break;
        }

        return mrcaRoot;
    }

    public int getDeep(Network<Object> net) {

        Map<NetNode<Object>, Integer> marks = new HashMap<>();
        for(NetNode<Object> node : net.bfs()) {
            marks.put(node, -1);
        }

        NetNode<Object> mrca = net.getRoot();
        while(mrca.getChildCount() == 1) {
            mrca = mrca.getChildren().iterator().next();
        }

        marks.put(mrca, 0);



        int count = Integer.MIN_VALUE;
        for(NetNode<Object> node : net.bfs()) {
            if(marks.get(node) >= 0) {
                for(NetNode<Object> child : node.getChildren()) {
                    marks.put(child, marks.get(node) + 1);
                }
            }
        }

        for(NetNode<Object> leaf : net.getLeaves()) {
            if(marks.get(leaf) > count)
                count = marks.get(leaf);
        }

        return count;
    }

    public List<NetNode> getConnectNodes(Network net, List<String> outLeaves) {
        List<NetNode> connectNodes = new ArrayList<>();
        Map<NetNode, Boolean> marks = new HashMap<>();

        return connectNodes;
    }

    class MountPoint {
        Tree _tree;
        MutableTuple<String, Integer> _position;

        MountPoint(Tree tree, String leaf, int up) {
            _tree = tree;
            _position = new MutableTuple<>(leaf, up);
        }

        @Override
        public boolean equals(Object candidate)
        {
            try
            {
                return equals((MountPoint) candidate);
            }
            catch (ClassCastException e)
            {
                return false;
            }
        }

        public boolean equals(MountPoint candidate)
        {
            if(candidate == null)
                return false;

            if(candidate._tree != _tree)
                return false;

            return isCongruence(candidate);
        }


        public boolean isCongruence(MountPoint candidate)
        {
            if(candidate == null)
                return false;

            if(candidate._tree != _tree)
                return false;

            TNode node1 = _tree.getNode(_position.Item1);
            for(int i = 0 ; i < _position.Item2 - 1 ; i++)
                node1 = node1.getParent();

            TNode node2 = _tree.getNode(candidate._position.Item1);
            for(int i = 0 ; i < candidate._position.Item2 - 1 ; i++)
                node2 = node2.getParent();

            return node1 == node2;
        }

        @Override
        public int hashCode()
        {
            TNode node = _tree.getNode(_position.Item1);
            for(int i = 0 ; i < _position.Item2 - 1 ; i++)
                node = node.getParent();
            return _tree.hashCode() ^ node.hashCode();
        }
    }

    private NetNode constructFromTripletsRecursively(Map<String, int[]> tripletCount, Set<String> names) {
        if(names.size() == 0) return null;
        if(names.size() == 1) {
            String name = names.iterator().next();
            NetNode node = new BniNetNode<>();
            node.setName(name);
            return node;
        }
        if(names.size() == 2) {
            Iterator<String> iter = names.iterator();
            NetNode root = new BniNetNode<>();
            NetNode child1 = new BniNetNode<>();
            child1.setName(iter.next());
            NetNode child2 = new BniNetNode<>();
            child2.setName(iter.next());
            root.adoptChild(child1, NetNode.NO_DISTANCE);
            root.adoptChild(child2, NetNode.NO_DISTANCE);
            return root;
        }

        DisjointSets<String> disjointSets = new DisjointSets<>(names);
        for(String triplet : tripletCount.keySet()) {
            String[] split = triplet.split("&");
            if(!names.containsAll(Arrays.asList(split)))
                continue;

            int[] nums = tripletCount.get(triplet);
            int nonzero = -1;
            for(int i = 0 ; i < 3 ; i++) {
                if(nums[i] > 0) {
                    if(nonzero == -1)
                        nonzero = i;
                    else
                        return null;
                }
            }
            if(nonzero == 0)
                disjointSets.union(split[0], split[1]);
        }

        Set<Set<String>> partition = disjointSets.getDisjointSets();
        if(partition.size() <= 1)
            return null;

        NetNode node = new BniNetNode<>();

        for(Set<String> group : partition) {
            NetNode child = constructFromTripletsRecursively(tripletCount, group);
            if(child == null)
                return null;
            node.adoptChild(child, NetNode.NO_DISTANCE);
        }

        return node;
    }

    public Network<Object> constructFromTriplets(Map<String, int[]> tripletCount) {
        Set<String> names = new HashSet<>();
        for(String triplet : tripletCount.keySet()) {
            String[] split = triplet.split("&");
            names.add(split[0]);
            names.add(split[1]);
            names.add(split[2]);
        }
        NetNode root = constructFromTripletsRecursively(tripletCount, names);
        if(root == null)
            return null;
        Network network = new BniNetwork((BniNetNode) root);

        return network;
    }

    protected Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }

    protected Func1<Network, Double> getScoreFunction(final Set<String> singleAlleleSpecies){
        return new Func1<Network, Double>() {
            public Double execute(Network speciesNetwork) {
                //System.out.println(speciesNetwork.toString());
                ParentalTreeOperation pto = new ParentalTreeOperation();
                List<Tree> parentalTrees = pto.getParentalTrees(speciesNetwork);
                Set<String> parentalTreesNewick = new HashSet<>();
                for(Tree t : parentalTrees) {
                    parentalTreesNewick.add(t.toNewick());
                }
                double score = 0;
                for(String s : _inputParentalTreesNewick) {
                    if(parentalTreesNewick.contains(s))
                        score += 1.0;
                }
                int count = 0;
                for(Object node : speciesNetwork.getNetworkNodes()) {
                    count++;
                }
                score += count > 0 ? 1.0 / count : 0;
                //System.out.println(score);
                return score;
            }
        };
    }

    public Network<Object> inferNetworkHeuristc(List<List<MutableTuple<Tree, Double>>> parentalTrees) {
        LinkedList<Tuple<Network,Double>> resultList = new LinkedList<>();
        _inputParentalTreesNewick = new HashSet<>();

        List<Tree> trees = new ArrayList<>();
        for(int i = 0 ; i < parentalTrees.size() ; i++) {
            STITree tree = null;
            try {
                tree = new STITree(parentalTrees.get(i).get(0).Item1.toNewick());
            } catch (Exception e){

            };
            Trees.removeBinaryNodes(tree);
            Trees.convertToLexicographicTree(tree);
            _inputParentalTreesNewick.add(tree.toNewick());
        }

        Set<String> singleAlleleSpecies = new HashSet<>();
        int maxReticulations = 3;
        int numSol = 1;

        String startingNetwork = inferNetwork(parentalTrees).toString();
        Network speciesNetwork = Networks.readNetwork(startingNetwork);
        for (Object nt : Networks.getTrees(speciesNetwork)) {
            speciesNetwork = Networks.readNetwork(((NetworkTree<Object>)nt).makeTree().toNewick());
            break;
        }

        NetworkRandomParameterNeighbourGenerator parameterGenerator = new NonUltrametricNetworkRandomParameterNeighbourGenerator(singleAlleleSpecies);
        NetworkRandomNeighbourGenerator allNeighboursStrategy = new NetworkRandomNeighbourGenerator(new NetworkRandomTopologyNeighbourGenerator(_topologyOperationWeights, maxReticulations, _moveDiameter, _reticulationDiameter, _fixedHybrid, _seed), _topologyVsParameterOperation[0], parameterGenerator, _topologyVsParameterOperation[1], _seed);
        Comparator<Double> comparator = getDoubleScoreComparator();

        HillClimberBase searcher;
        //searcher = new SimpleHillClimbing(comparator, allNeighboursStrategy);
        searcher = new SimulatedAnnealingSalterPearL(comparator, allNeighboursStrategy, _seed);

        searcher.setLogFile(_logFile);
        searcher.setIntermediateResultFile(_intermediateResultFile);
        Func1<Network, Double> scorer = getScoreFunction(singleAlleleSpecies);

        searcher.search(speciesNetwork, scorer, numSol, _numRuns, _maxExaminations, _maxFailure, true, resultList); // search starts here

        return resultList.get(0).Item1;
    }

    public Network<Object> inferNetwork(List<List<MutableTuple<Tree, Double>>> parentalTrees) {
        //Find backbone tree among parental trees

        Network<Object> network = null;
        List<Tree> trees = new ArrayList<>();
        for(int i = 0 ; i < parentalTrees.size() ; i++) {
            trees.add(parentalTrees.get(i).get(0).Item1);
        }

        //SteelWarnowMAST steelWarnowMAST = new SteelWarnowMAST();
        //Tree mast = steelWarnowMAST.computeRMAST(trees);

        AmirKeselmanMAST amirKeselmanMAST = new AmirKeselmanMAST();
        Tree mast = amirKeselmanMAST.computeRMAST(trees);
        mast = Trees.readTree(mast.toNewick());

        List<String> outLeaves = new ArrayList<>();
        for(String leaf : mast.getLeaves())
            outLeaves.add(leaf);

        System.out.println("MAST: " + mast.toNewick());

        network = Networks.readNetwork(mast.toNewick());

        Map<Set<String>, Integer> numberMountPointsPerGroup = new HashMap<>();
        Map<MountPoint, List<List<Set<String>>>> mountPointsWithLeaves = new HashMap<>();
        Map<MountPoint, List<Set<String>>> orderOnMountPoint = new HashMap<>();
        List<MutableTuple<Set<String>, Set<MountPoint>>> groupedInLeaves = new ArrayList<>();

        for(int i = 0 ; i < parentalTrees.size() ; i++) {
            Map<MountPoint, List<Set<String>>> currentMountPointsWithLeaves = new HashMap<>();
            Network<Object> net = Networks.readNetwork(parentalTrees.get(i).get(0).Item1.toNewick());

            for(NetNode<Object> node : net.bfs()) {
                if(!node.isLeaf()) {
                    node.setName("");
                }
                for(NetNode<Object> parent: node.getParents()){
                    node.setParentDistance(parent,NetNode.NO_DISTANCE);
                }
            }


            Map<NetNode<Object>, Integer> marks = new HashMap<>();

            int OUT = 1;
            int IN = 2;
            int DONE = 3;

            for(NetNode<Object> node : net.bfs()) {
                marks.put(node, IN);
            }

            marks.put(net.getRoot(), OUT);

            for(NetNode<Object> leaf : net.getLeaves()) {
                if(outLeaves.contains(leaf.getName())) {
                    NetNode<Object> node = leaf;
                    while(!node.isRoot()) {
                        marks.put(node, OUT);
                        for(NetNode<Object> parent : node.getParents()) {
                            node = parent;
                            break;
                        }
                    }
                }
            }

            //Find all group in current parental tree
            for(NetNode<Object> outNode : net.bfs()) {
                Set<String> currentGroup = new HashSet<>();
                Set<MountPoint> currentMountPoints = new HashSet<>();

                if(marks.get(outNode) == OUT) {
                    NetNode<Object> inNode = null;
                    for(NetNode<Object> child : outNode.getChildren()) {
                        if(marks.get(child) == IN) {
                            inNode = child;
                            break;
                        }
                    }
                    if(inNode == null)
                        continue;

                    //Find all leaves in current group in current parental tree
                    Queue<NetNode<Object>> queue = new LinkedList<>();
                    queue.offer(inNode);
                    marks.put(inNode, DONE);

                    while(queue.size() > 0) {
                        NetNode<Object> node = queue.poll();
                        if(node.isLeaf()) {
                            currentGroup.add(node.getName());
                        }
                        for(NetNode<Object> child : node.getChildren()) {
                            queue.offer(child);
                            marks.put(child, DONE);
                        }
                    }

                    //Find Mount Point of current group
                    int count = 1;
                    NetNode<Object> outLeaf = inNode.getParents().iterator().next();
                    while(!outLeaf.isLeaf()) {
                        int cc = 0;
                        for(NetNode<Object> son : outLeaf.getChildren()) {
                            if(marks.get(son) == OUT) {
                                outLeaf = son;
                                cc++;
                            }
                        }
                        //if out-degree > 1 in MAST
                        if(cc > 1)
                            count++;
                    }
                    MountPoint currentMountPoint = new MountPoint(mast, outLeaf.getName(), count);
                    currentMountPoints.add(currentMountPoint);

                    //Find corresponding Order of current Mount Point
                    List<Set<String>> order = null;
                    for(MountPoint mountPoint : currentMountPointsWithLeaves.keySet()) {
                        if(mountPoint.isCongruence(currentMountPoint)) {
                            order = currentMountPointsWithLeaves.get(mountPoint);
                            break;
                        }
                    }

                    if(order == null) {
                        currentMountPointsWithLeaves.put(currentMountPoint, new ArrayList<>());
                        order = currentMountPointsWithLeaves.get(currentMountPoint);
                    }

                    order.add(new HashSet<String>(currentGroup));

                    //register current Group
                    boolean inserted = false;
                    //combine sets of leaves with intersection
                    for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {
                        Set<String> group = tuple.Item1;
                        Set<String> intersection = new HashSet<String>(group);
                        intersection.retainAll(currentGroup);
                        if(intersection.size() > 0) {
                            inserted = true;
                            group.addAll(currentGroup);
                            tuple.Item2.addAll(currentMountPoints);
                        }
                    }

                    if(inserted) {
                        boolean done = false;
                        while(!done) {
                            done = true;
                            for(MutableTuple<Set<String>, Set<MountPoint>> tuple1 : groupedInLeaves) {
                                Set<String> group1 = tuple1.Item1;
                                for(MutableTuple<Set<String>, Set<MountPoint>> tuple2 : groupedInLeaves) {
                                    Set<String> group2 = tuple2.Item1;
                                    if(group1 != group2) {
                                        Set<String> intersection = new HashSet<String>(group1);
                                        intersection.retainAll(group2);
                                        if(intersection.size() > 0) {
                                            group1.addAll(group2);
                                            tuple1.Item2.addAll(tuple2.Item2);
                                            groupedInLeaves.remove(tuple2);
                                            done = false;
                                            break;
                                        }
                                    }
                                }
                                if(!done)
                                    break;
                            }
                        }
                    }

                    if(!inserted && currentGroup.size() > 0)
                        groupedInLeaves.add(new MutableTuple<Set<String>, Set<MountPoint>>(currentGroup, currentMountPoints));
                }
            }

            //combine current Order with Mount Point in current parental tree with global one
            for(MountPoint mountPoint1 : currentMountPointsWithLeaves.keySet()) {
                boolean found = false;
                for(MountPoint mountPoint2 : mountPointsWithLeaves.keySet()) {
                    if(mountPoint1.isCongruence(mountPoint2)) {
                        found = true;
                        mountPointsWithLeaves.get(mountPoint2).add(currentMountPointsWithLeaves.get(mountPoint1));
                        break;
                    }
                }

                if(!found) {
                    mountPointsWithLeaves.put(mountPoint1, new ArrayList<>());
                    mountPointsWithLeaves.get(mountPoint1).add(currentMountPointsWithLeaves.get(mountPoint1));
                }
            }


        }

        Map<String, Set<String>> groupFinder = new HashMap<>();
        boolean isLevel1 = true;

        //build Group finder
        for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {

            for(String leaf : tuple.Item1) {
                groupFinder.put(leaf, tuple.Item1);
            }

        }

        Map<Set<String>, Set<Set<String>>> subgroups = new IdentityHashMap<>();
        for(MountPoint mountPoint : mountPointsWithLeaves.keySet()) {
            //get number of SubGroup of one Group shown in current Mount Point
            //amount all Orders, keep the largest number
            Map<Set<String>, Integer> numberShown = new IdentityHashMap<>();
            for(List<Set<String>> order : mountPointsWithLeaves.get(mountPoint)) {
                Map<Set<String>, Integer> currentNumberShown = new IdentityHashMap<>();
                for(Set<String> set : order) {
                    Set<String> groupBelong = groupFinder.get(set.iterator().next());
                    if(currentNumberShown.get(groupBelong) == null)
                        currentNumberShown.put(groupBelong, 0);
                    currentNumberShown.put(groupBelong, currentNumberShown.get(groupBelong) + 1);
                }
                for(Set<String> group : currentNumberShown.keySet()) {
                    if(numberShown.get(group) == null)
                        numberShown.put(group, 0);
                    if(numberShown.get(group) < currentNumberShown.get(group))
                        numberShown.put(group, currentNumberShown.get(group));
                }

            }

            //initialize the SubGroups of one Group shown in current Mount Point
            //paritalOrder is the graph of Orders
            Map<Set<String>, Map<Set<String>, Integer>> partialOrder = new IdentityHashMap<>();
            Map<Set<String>, List<Set<String>>> subGroupOnMountPoint = new IdentityHashMap<>();
            for(Set<String> group : numberShown.keySet()) {
                subGroupOnMountPoint.put(group, new ArrayList<>());
                for(int i = 0 ; i < numberShown.get(group) ; i++) {
                    subGroupOnMountPoint.get(group).add(new HashSet<>());
                    partialOrder.put(subGroupOnMountPoint.get(group).get(i), new IdentityHashMap<>());
                }

            }

            //get all SubGroups of one Group shown in current Mount Point
            for(List<Set<String>> order : mountPointsWithLeaves.get(mountPoint)) {
                Map<Set<String>, Integer> currentNumberShown = new IdentityHashMap<>();
                Map<Set<String>, List<Set<String>>> curSubGroupOnMountPoint = new IdentityHashMap<>();

                //get all SubGroups of one Group shown in current Order in current Mount Point
                for(Set<String> set : order) {
                    Set<String> groupBelong = groupFinder.get(set.iterator().next());
                    if(currentNumberShown.get(groupBelong) == null)
                        currentNumberShown.put(groupBelong, 0);
                    currentNumberShown.put(groupBelong, currentNumberShown.get(groupBelong) + 1);
                    if(curSubGroupOnMountPoint.get(groupBelong) == null)
                        curSubGroupOnMountPoint.put(groupBelong, new ArrayList<>());
                    curSubGroupOnMountPoint.get(groupBelong).add(set);
                }

                //when number of SubGroup of one Group is correct in this Order -> update
                //else try to update
                for(Set<String> group : numberShown.keySet()) {
                    if(numberShown.get(group).equals(currentNumberShown.get(group))) {
                        for(int i = 0 ; i < numberShown.get(group) ; i++)
                            subGroupOnMountPoint.get(group).get(i).addAll(curSubGroupOnMountPoint.get(group).get(i));
                    } else if (curSubGroupOnMountPoint.containsKey(group)){
                        for(int i = 0 ; i < curSubGroupOnMountPoint.get(group).size() ; i++) {
                            int j = 0;
                            while(j < numberShown.get(group)) {
                                if(curSubGroupOnMountPoint.get(group).get(i).containsAll(subGroupOnMountPoint.get(group).get(j))) {
                                    subGroupOnMountPoint.get(group).get(j).addAll(curSubGroupOnMountPoint.get(group).get(i));
                                    break;
                                }
                                j++;
                            }
                        }
                    }

                }

                //when number of SubGroup of one Group is correct in this Order -> put into graph of Orders
                Map<Set<String>, Integer> currentScanNumberShown = new IdentityHashMap<>();
                Set<String> prevset = null;
                for(Set<String> set : order) {
                    Set<String> groupBelong = groupFinder.get(set.iterator().next());
                    if(currentScanNumberShown.get(groupBelong) == null)
                        currentScanNumberShown.put(groupBelong, 0);
                    Set<String> curset = subGroupOnMountPoint.get(groupBelong).get(currentScanNumberShown.get(groupBelong));
                    currentScanNumberShown.put(groupBelong, currentScanNumberShown.get(groupBelong) + 1);

                    if (numberShown.get(groupBelong).equals(currentNumberShown.get(groupBelong))) {

                        if(prevset != null) {
                            partialOrder.get(prevset).put(curset, 1);
                        }
                        prevset = curset;
                    }
                }




            }

            for(Set<String> group : subGroupOnMountPoint.keySet()) {
                if(subgroups.get(group) == null)
                    subgroups.put(group, new HashSet<>());
                subgroups.get(group).addAll(subGroupOnMountPoint.get(group));
            }

            //do topological sorting to get the correct Order
            List<Set<String>> currentOrder = new ArrayList<>();
            int totalSubMountPoints = partialOrder.size();
            Map<Set<String>, Integer> indegree = new IdentityHashMap<>();

            for(Set<String> set : partialOrder.keySet()) {
                indegree.put(set, 0);
            }

            for(Set<String> set : partialOrder.keySet()) {
                for(Set<String> nextset : partialOrder.get(set).keySet()) {
                    indegree.put(nextset, indegree.get(nextset) + 1);
                }
            }

            int subMountPointsCount = totalSubMountPoints;
            while(subMountPointsCount > 0) {
                boolean found = false;

                for(Set<String> set : partialOrder.keySet()) {
                    if(indegree.get(set) == 0) {
                        found = true;
                        currentOrder.add(set);
                        for(Set<String> nextset : partialOrder.get(set).keySet()) {
                            indegree.put(nextset, indegree.get(nextset) - 1);
                        }
                        indegree.put(set, -1);
                        subMountPointsCount--;
                        break;
                    }
                }

                if(subMountPointsCount == 0) break;
                //when there is a conflict, do breaking
                if(!found) {
                    for(Set<String> set : partialOrder.keySet()) {
                        if(indegree.get(set) > 0) {
                            found = true;
                            Set<String> newset = new HashSet<>(set);
                            indegree.put(newset, 0);
                            partialOrder.put(newset, new HashMap(partialOrder.get(set)));
                            partialOrder.get(set).clear();
                            Set<String> groupBelong = groupFinder.get(set.iterator().next());
                            numberShown.put(groupBelong, numberShown.get(groupBelong) + 1);
                            subMountPointsCount++;
                            totalSubMountPoints++;
                            break;
                        }
                    }
                }
            }

            for(Set<String> subgroup : currentOrder) {
                Set<String> groupBelong = groupFinder.get(subgroup.iterator().next());
                if(numberMountPointsPerGroup.get(groupBelong) == null)
                    numberMountPointsPerGroup.put(groupBelong, 0);
                numberMountPointsPerGroup.put(groupBelong, numberMountPointsPerGroup.get(groupBelong) + 1);
            }
            orderOnMountPoint.put(mountPoint, currentOrder);


        }

        System.out.println("Groups detected: " + groupedInLeaves.size());


        //check if every Group is associated with two Mount Point
        for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {

            if(numberMountPointsPerGroup.get(tuple.Item1) > 2)
                isLevel1 = false;
        }

        if(!isLevel1) {
            System.out.println("I don't think it is a level-1 network.");
        }

        if(false) {
            //remove duplicated mount points associated with one set of leaves
            for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {
                boolean done = false;
                while(!done) {
                    done = true;
                    for (MountPoint mountPoint1 : tuple.Item2) {
                        for (MountPoint mountPoint2 : tuple.Item2) {
                            if (mountPoint1 != mountPoint2 && mountPoint2.isCongruence(mountPoint1)) {
                                tuple.Item2.remove(mountPoint2);
                                done = false;
                                break;
                            }
                        }
                        if(!done)
                            break;
                    }
                }

            }

            //combine groups with identical set of mount points
            boolean combineDone = false;
            while(!combineDone) {
                combineDone = true;
                for (MutableTuple<Set<String>, Set<MountPoint>> tuple1 : groupedInLeaves) {
                    for (MutableTuple<Set<String>, Set<MountPoint>> tuple2 : groupedInLeaves) {
                        if(tuple1 != tuple2) {
                            boolean setCongruence = true;
                            for (MountPoint mountPoint1 : tuple1.Item2) {
                                boolean found = false;
                                for (MountPoint mountPoint2 : tuple2.Item2) {
                                    if (mountPoint1.isCongruence(mountPoint2)) {
                                        found = true;
                                        break;
                                    }
                                }
                                if (!found) {
                                    setCongruence = false;
                                    break;
                                }
                            }
                            if (setCongruence) {
                                combineDone = false;
                                tuple1.Item1.addAll(tuple2.Item1);
                                groupedInLeaves.remove(tuple2);
                                break;
                            }
                        }
                    }
                    if(!combineDone)
                        break;
                }
            }

            for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {
                List<String> currentInLeaves = new ArrayList<>(tuple.Item1);

                Network<Object> reticulatePart = null;
                Network<Object> mostFreqPart = null;
                double mostFreq = -1;

                for (int i = 0; i < parentalTrees.size(); i++) {
                    Network<Object> net = Networks.readNetwork(parentalTrees.get(i).get(0).Item1.toNewick());
                    int count = getSubnetworkRetain(net, currentInLeaves);
                    removeBinaryNodes(net);

                    if (count == 1) {
                        if (reticulatePart == null) {
                            reticulatePart = net;
                            removeBinaryNodes(reticulatePart);
                        }
                    } else {
                        if (mostFreq < parentalTrees.get(i).get(0).Item2) {
                            mostFreq = parentalTrees.get(i).get(0).Item2;
                            mostFreqPart = net;
                            removeBinaryNodes(mostFreqPart);
                        }
                    }
                }

                if (reticulatePart == null)
                    reticulatePart = mostFreqPart;

                List<NetNode<Object>> mountNodes = new ArrayList<>();

                for(MountPoint mountPoint : tuple.Item2) {
                    NetNode<Object> node = network.findNode(mountPoint._position.Item1);
                    for(int i = 0 ; i < mountPoint._position.Item2 - 1; i++) {
                        for(NetNode<Object> parent : node.getParents()) {
                            node = parent;
                            break;
                        }
                    }
                    mountNodes.add(node);
                }

                for(NetNode<Object> mountNode : mountNodes) {
                    NetNode<Object> newnode = (NetNode<Object>) new BniNetNode<Object>();
                    NetNode<Object> parent = mountNode.getParents().iterator().next();
                    parent.removeChild(mountNode);
                    parent.adoptChild(newnode, NetNode.NO_DISTANCE);
                    newnode.adoptChild(mountNode, NetNode.NO_DISTANCE);
                    newnode.adoptChild(reticulatePart.getRoot(), NetNode.NO_DISTANCE);
                }


            }
        } else {

            Map<NetNode<Object>, Integer> marks = new HashMap<>();
            Map<MountPoint, NetNode<Object>> mountNodes = new HashMap<>();
            int ORIGINAL = 0;
            int NEW = 1;
            for(NetNode<Object> node : network.bfs()) {
                marks.put(node, ORIGINAL);
            }

            for(MountPoint mountPoint : orderOnMountPoint.keySet()) {
                NetNode<Object> node = network.findNode(mountPoint._position.Item1);
                for(int i = 0 ; i < mountPoint._position.Item2 - 1; i++) {
                    do {
                        node = node.getParents().iterator().next();
                    } while(marks.get(node) == NEW);
                }

                for(int i = 0 ; i < orderOnMountPoint.get(mountPoint).size() ; i++) {
                    NetNode<Object> newnode = (NetNode<Object>) new BniNetNode<Object>();
                    if(node == network.getRoot()) {
                        newnode.adoptChild(node, NetNode.NO_DISTANCE);
                        network.resetRoot(newnode);
                    } else {
                        NetNode<Object> parent = node.getParents().iterator().next();
                        parent.removeChild(node);
                        parent.adoptChild(newnode, NetNode.NO_DISTANCE);
                        newnode.adoptChild(node, NetNode.NO_DISTANCE);
                        //mountNodes.put()
                    }
                    marks.put(newnode, NEW);
                }
            }


            for(MutableTuple<Set<String>, Set<MountPoint>> tuple : groupedInLeaves) {
                List<String> currentInLeaves = new ArrayList<>(tuple.Item1);

                Network<Object> reticulatePart = null;
                Network<Object> bestFairPart = null;
                double mostFreq = -1;
                int minDeep = Integer.MAX_VALUE;

                if(_verbose) {
                    for (String s : tuple.Item1)
                        System.out.print(s + ' ');
                    System.out.println();

                    System.out.println("Subgroup:");
                    for (Set<String> set : subgroups.get(tuple.Item1)) {
                        for (String s : set) {
                            System.out.print(s + ' ');
                        }
                        System.out.println();
                    }
                }

                List<Tree> subTreesForGroup = new ArrayList<>();
                //find the best reticulated part of one Group
                for (int i = 0; i < parentalTrees.size(); i++) {
                    if(_verbose)
                        System.out.println("Parental Tree #" + i + " " + parentalTrees.get(i).get(0).toString());

                    Network<Object> net = Networks.readNetwork(parentalTrees.get(i).get(0).Item1.toNewick());
                    for(NetNode<Object> node : net.bfs()) {
                        if(!node.isLeaf()) {
                            node.setName("");
                        }
                        for(NetNode<Object> parent: node.getParents()){
                            node.setParentDistance(parent,NetNode.NO_DISTANCE);
                        }
                    }

                    int count = getSubnetworkRetain(net, currentInLeaves);
                    removeBinaryNodes(net);


                    if(_verbose)
                        System.out.println("--> " + net.toString() + " " + count);


                    if (count == 1) {
                        if (reticulatePart == null) {
                            reticulatePart = net;
                        }
                        subTreesForGroup.add(Trees.readTree(net.toString()));

                    } else {
                        int deep = 0;
                        for(Set<String> subgroup : subgroups.get(tuple.Item1)) {
                            Network<Object> subnet = Networks.readNetwork(net.toString());
                            getSubnetworkRetain(subnet, new ArrayList<>(subgroup));
                            deep += getDeep(subnet);
                        }
                        if(_verbose)
                            System.out.println("Score: " + deep);
                        if (deep < minDeep) {
                            minDeep = deep;
                            bestFairPart = net;
                        }
                    }


                }
                Map<String, int[]> tripletCount = new HashMap<>();

                for (int i = 0; i < parentalTrees.size(); i++) {
                    if(_verbose)
                        System.out.println("Parental Tree #" + i + " " + parentalTrees.get(i).get(0).toString());

                    Tree curPTree = Trees.readTree(parentalTrees.get(i).get(0).Item1.toNewick());


                    int count = getAllTriplets(curPTree, currentInLeaves, tripletCount);

                }

                if(tuple.Item1.size() > 2)
                    reticulatePart = constructFromTriplets(tripletCount);
                if(_verbose)
                    System.out.println("Compute from triplets: " + (reticulatePart != null ? reticulatePart.toString() : "null"));

                if(tuple.Item1.size() > 2 && reticulatePart == null && subTreesForGroup.size() > 0) {
                    AmirKeselmanMAST submastObj = new AmirKeselmanMAST();
                    Tree submast = submastObj.computeRMAST(subTreesForGroup);
                    submast = Trees.readTree(submast.toNewick());
                    if (submast.getLeaves().length < tuple.Item1.size()) {
                        List<List<MutableTuple<Tree, Double>>> subpts = new ArrayList<>();
                        for (Tree tree : subTreesForGroup) {
                            subpts.add(Arrays.asList(new MutableTuple<Tree, Double>(Trees.readTree(tree.toNewick()), new Double(1.0))));

                        }
                        reticulatePart = inferNetwork(subpts);
                    }
                }

                if (reticulatePart == null)
                    reticulatePart = bestFairPart;

                if(_verbose)
                    System.out.println("Best reticulate part: " + reticulatePart.toString());

                /*NetNode<Object> tmpRoot = new BniNetNode<>();
                tmpRoot.adoptChild(network.getRoot(), NetNode.NO_DISTANCE);
                tmpRoot.adoptChild(reticulatePart.getRoot(), NetNode.NO_DISTANCE);
                network.resetRoot(tmpRoot);*/

                List<Tuple> mountInst = new ArrayList<>();
                //mount best reticulated part to MAST
                for(MountPoint mountPoint : tuple.Item2) {
                    NetNode<Object> node = network.findNode(mountPoint._position.Item1);
                    for(int i = 0 ; i < mountPoint._position.Item2 - 1; i++) {
                        do {
                            node = node.getParents().iterator().next();
                        } while(marks.get(node) == NEW);
                    }

                    int subMountPoint = orderOnMountPoint.get(mountPoint).size();
                    for(Set<String> subgroup : orderOnMountPoint.get(mountPoint)) {
                        if(tuple.Item1.containsAll(subgroup)) {
                            NetNode<Object> connectNode = node;
                            for(int i = 0 ; i < subMountPoint ; i++)
                                connectNode = connectNode.getParents().iterator().next();



                            /*if(reticulatePart.getRoot().getParentCount() == 2) {
                                NetNode<Object> newroot = new BniNetNode<>();
                                newroot.adoptChild(reticulatePart.getRoot(), NetNode.NO_DISTANCE);
                                reticulatePart.resetRoot(newroot);
                            }*/


                            if(tuple.Item1.equals(subgroup) && reticulatePart.getRoot() != network.getRoot()) {
                                NetNode<Object> newRetiRoot = new BniNetNode<>();
                                NetNode<Object> curRetiRoot = reticulatePart.getRoot();
                                //NetNode<Object> retiParent = curRetiRoot.getParents().iterator().next();
                                //retiParent.removeChild(curRetiRoot);
                                //retiParent.adoptChild(newRetiRoot, NetNode.NO_DISTANCE);
                                newRetiRoot.adoptChild(curRetiRoot, NetNode.NO_DISTANCE);
                                reticulatePart.resetRoot(newRetiRoot);
                                mountInst.add(new Tuple(connectNode, reticulatePart.getRoot()));
                                //connectNode.adoptChild(reticulatePart.getRoot(), NetNode.NO_DISTANCE);
                            } else {

                                NetNode<Object> newConnectNode = new BniNetNode<>();
                                NetNode<Object> reticulateNode = getMRCA(reticulatePart, new ArrayList<>(subgroup));
                                if (reticulateNode.isRoot()) {
                                    newConnectNode.adoptChild(reticulateNode, NetNode.NO_DISTANCE);
                                    reticulatePart.resetRoot(newConnectNode);
                                } else {
                                    NetNode<Object> parent = reticulateNode.getParents().iterator().next();
                                    parent.removeChild(reticulateNode);
                                    parent.adoptChild(newConnectNode, NetNode.NO_DISTANCE);
                                    newConnectNode.adoptChild(reticulateNode, NetNode.NO_DISTANCE);
                                }
                                mountInst.add(new Tuple(connectNode, newConnectNode));
                                //connectNode.adoptChild(newConnectNode, NetNode.NO_DISTANCE);
                            }
                        }
                        subMountPoint--;
                    }
                }

                //Deal with root of reticulate part first : otherwise Maybe bug when encounter partial SubGroup first
                for(Tuple<NetNode, NetNode> curMountInst : mountInst) {
                    if(curMountInst.Item2 == reticulatePart.getRoot()) {
                        curMountInst.Item1.adoptChild(curMountInst.Item2, NetNode.NO_DISTANCE);
                        mountInst.remove(curMountInst);
                        break;
                    }
                }
                for(Tuple<NetNode, NetNode> curMountInst : mountInst) {
                    curMountInst.Item1.adoptChild(curMountInst.Item2, NetNode.NO_DISTANCE);
                }
                //tmpRoot.removeChild(reticulatePart.getRoot());
                //network.resetRoot(tmpRoot.getChildren().iterator().next());
                //tmpRoot.removeChild(network.getRoot());
            }
        }

        removeBinaryNodes(network);
        return network;
    }

    public static void main(String[] args) {

        Network trueNetwork = Networks.readNetwork("((((((C)#H1,D)#H2,E),F)#H3,B), (((A,#H2),#H1),#H3));");
        ParentalTreeOperation parentalTreeOperation = new ParentalTreeOperation();
        InferNetworkFromParentalTrees inferNetworkFromParentalTrees = new InferNetworkFromParentalTrees();
        List<Tree> parentalTrees = parentalTreeOperation.getParentalTrees(trueNetwork);
        List<List<MutableTuple<Tree, Double>>> parentalTrees0 = new ArrayList<>();
        for(Tree tree : parentalTrees) {
            parentalTrees0.add(Arrays.asList(new MutableTuple<Tree, Double>(Trees.readTree(tree.toNewick()), new Double(1.0))));

        }
        Network inferredNetwork = inferNetworkFromParentalTrees.inferNetwork(parentalTrees0);
        System.out.println(inferredNetwork);

        List<Tree> parentalTreesI = parentalTreeOperation.getParentalTrees(inferredNetwork);
        List<List<MutableTuple<Tree, Double>>> parentalTreesI0 = new ArrayList<>();
        for(Tree tree : parentalTreesI) {
            parentalTreesI0.add(Arrays.asList(new MutableTuple<Tree, Double>(Trees.readTree(tree.toNewick()), new Double(1.0))));

        }

        RECOMB_CG_16_JZ fu = new RECOMB_CG_16_JZ();
        System.out.println(fu.getBestMatchDistanceBetweenTrees(parentalTrees0, parentalTreesI0));
    }
}
