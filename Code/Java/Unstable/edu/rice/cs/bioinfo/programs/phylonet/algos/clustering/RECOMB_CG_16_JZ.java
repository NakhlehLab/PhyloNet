package edu.rice.cs.bioinfo.programs.phylonet.algos.clustering;

import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MajorityConsensusInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.mast.SteelWarnowMAST;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.BipartiteGraph;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.TimeUnit;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 7/21/16
 * Time: 11:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class RECOMB_CG_16_JZ {

    private String _msdir = "/scratch/jz55/Luay/msdir/ms";
    private String _treedistPath = "/scratch/jz55/Luay/phylip-3.695/src";


    public List<List<MutableTuple<Tree, Double>>> loadGTs(String filename) {
        List<List<MutableTuple<Tree, Double>>> result = new ArrayList<>();
        try {
            File file = new File(filename);
            Scanner scanner = new Scanner(file);

            while(scanner.hasNextLine()) {
                String s = scanner.nextLine();
                if(s.equals(""))
                    break;
                result.add(Arrays.asList(new MutableTuple(Trees.readTree(s), 1.0)));
            }

            scanner.close();

        }catch (Exception e) {
            e.printStackTrace();
        }
        return result;
    }

    public Network<Object> loadNetwork(String filename) {
        Network<Object> result = null;
        try {
            File file = new File(filename);
            Scanner scanner = new Scanner(file);

            String s = scanner.nextLine();
            result = Networks.readNetwork(s);

            scanner.close();

        }catch (Exception e) {
            e.printStackTrace();
        }
        return result;
    }

    public void loadMatrix(String filename, int n, double [][] matrix) {
        try {
            File file = new File(filename);
            Scanner scanner = new Scanner(file);

            int total = scanner.nextInt();
            if(total < n)
                throw new RuntimeException("Matrix is not big enough!");

            for(int i = 0 ; i < n ; i++) {
                for (int j = 0; j < n; j++)
                    matrix[i][j] = scanner.nextDouble();
                scanner.nextLine();
            }

            scanner.close();
        }catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getDistancePHYLIP(List<List<MutableTuple<Tree, Double>>> gts, double[][] distMatrix) {
        String suffix = Long.toString(System.nanoTime()) + "_" + Thread.currentThread().getId() + "_" + (int)(Math.random() * 1000000);
        String path = _treedistPath;
        int n = gts.size();
        try {
            PrintWriter out = new PrintWriter("intree" + suffix);
            for(int i = 0 ; i < n ; i++)
                out.println(gts.get(i).get(0).Item1.toString());
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        try {
            File f = new File("outfile" + suffix);
            f.delete();
        } catch (Exception e) {
        }

        try {
            Process proc = Runtime.getRuntime().exec(path + "/treedist " + "intree" + suffix + " outfile" + suffix, null, null);
            proc.waitFor();
            File file = new File("outfile" + suffix);
            Scanner scanner = new Scanner(file);

            for(int k = 0 ; k < n * n ; k++) {
                int i = scanner.nextInt();
                int j = scanner.nextInt();
                double dist = scanner.nextDouble();
                distMatrix[i-1][j-1] = dist;
            }

            scanner.close();

            new File("outfile" + suffix).delete();
            new File("intree" + suffix).delete();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public double getDistance(Tree tree1, Tree tree2) {
        SymmetricDifference symmetricDifference = new SymmetricDifference();
        symmetricDifference.computeDifference(tree1, tree2, true);
        double diff = symmetricDifference.getWeightedAverage();
        return diff;
    }

    public double getDistance(Tree tree1, Tree tree2, boolean rooted) {
        SymmetricDifference symmetricDifference = new SymmetricDifference();
        symmetricDifference.computeDifference(tree1, tree2, rooted);
        double diff = symmetricDifference.getWeightedAverage();
        return diff;
    }

    public double getBestMatchDistanceBetweenTrees(List<List<MutableTuple<Tree, Double>>> trees1, List<List<MutableTuple<Tree, Double>>> trees2) {
        int size1 = trees1.size();
        int size2 = trees2.size();

        double distMatrix[][] = new double[size1][size2];

        for(int i = 0 ; i < size1 ; i++)
            for(int j = 0 ; j < size2 ; j++)
                distMatrix[i][j] = getDistance(trees1.get(i).get(0).Item1, trees2.get(j).get(0).Item1, true);

        HungarianMatching HungarianMatching = new HungarianMatching(distMatrix);
        int[] assignment = HungarianMatching.compute();
        double score = 0;
        for(int i = 0 ; i < assignment.length ; i++)
            if(assignment[i] != -1)
                score += distMatrix[i][assignment[i]];

        //double score2 = getBestMatchDistanceBetweenTreesRecursively(trees1, trees2, distMatrix, Math.min(size1, size2));
        return score;

    }

    public double getBestMatchDistanceBetweenTrees(List<Tree> trees1, List<Tree> trees2, boolean rooted) {
        int size1 = trees1.size();
        int size2 = trees2.size();

        double distMatrix[][] = new double[size1][size2];

        for(int i = 0 ; i < size1 ; i++)
            for(int j = 0 ; j < size2 ; j++)
                distMatrix[i][j] = getDistance(trees1.get(i), trees2.get(j), rooted);

        HungarianMatching HungarianMatching = new HungarianMatching(distMatrix);
        int[] assignment = HungarianMatching.compute();
        double score = 0;
        for(int i = 0 ; i < assignment.length ; i++)
            if(assignment[i] != -1)
                score += distMatrix[i][assignment[i]];

        //double score2 = getBestMatchDistanceBetweenTreesRecursively(trees1, trees2, distMatrix, Math.min(size1, size2));
        return score;

    }

    public void getDistanceMatrixBestTree(List<List<MutableTuple<Tree, Double>>> gts, String path, double[][] distMatrix) {
        int size = gts.size();
        List<Tree> bestTrees = new ArrayList<>();
        try {
            for(int i = 0 ; i < size ; i++) {

                File file = new File(path + "/seq/" + i + ".bestTree");
                Scanner scanner = new Scanner(file);
                String s = scanner.nextLine();
                Tree newtree = Trees.readTree(s);
                bestTrees.add(newtree);

                scanner.close();
            }


        } catch (IOException e) {
            e.printStackTrace();
        }


        for(int i = 0 ; i < size ; i++) {
            for(int j = i + 1 ; j < size ; j++) {
                double distance = getDistance(bestTrees.get(i), bestTrees.get(j), false);
                distMatrix[i][j] = distMatrix[j][i] = distance;
            }
        }
    }

    public void getDistanceMatrixBootstrap(List<List<MutableTuple<Tree, Double>>> gts, String path, double[][] distMatrix) {
        int size = gts.size();
        List<List<Tree>> bootstrapTrees = new ArrayList<>();
        try {
            for(int i = 0 ; i < size ; i++) {
                bootstrapTrees.add(new ArrayList<>());

                File file = new File(path + "/seq/" + i + ".bootstrap");
                Scanner scanner = new Scanner(file);
                while(scanner.hasNextLine()) {
                    String s = scanner.nextLine();
                    if(s.length() > 1) {
                        Tree newtree = Trees.readTree(s);
                        bootstrapTrees.get(i).add(newtree);
                    }
                }

                scanner.close();
            }


        } catch (IOException e) {
            e.printStackTrace();
        }


        for(int i = 0 ; i < size ; i++) {
            for(int j = i + 1 ; j < size ; j++) {
                //System.out.println("Computing gene " + i + "(" + bootstrapTrees.get(i).size() + ")" + " vs gene " + j + "(" + bootstrapTrees.get(j).size() + ")");
                //double distance = Networks.computeTreeDistance(bootstrapTrees.get(i), bootstrapTrees.get(j))[2];
                double distance = getBestMatchDistanceBetweenTrees(bootstrapTrees.get(i), bootstrapTrees.get(j), false);
                distMatrix[i][j] = distMatrix[j][i] = distance;
            }
        }
    }

    public void getDistanceMatrixBootstrapConsensus(List<List<MutableTuple<Tree, Double>>> gts, String path, double[][] distMatrix) {
        int size = gts.size();
        List<List<Tree>> bootstrapTrees = new ArrayList<>();
        List<Tree> consensusTrees = new ArrayList<>();
        try {
            for(int i = 0 ; i < size ; i++) {
                List<List<MutableTuple<Tree, Double>>> treeList = new ArrayList<>();
                bootstrapTrees.add(new ArrayList<>());

                File file = new File(path + "/seq/" + i + ".bootstrap");
                Scanner scanner = new Scanner(file);
                while(scanner.hasNextLine()) {
                    String s = scanner.nextLine();
                    if(s.length() > 1) {
                        Tree newtree = Trees.readTree(s);
                        bootstrapTrees.get(i).add(newtree);
                        treeList.add(Arrays.asList(new MutableTuple<Tree, Double>(newtree, 1.0)));
                    }
                }

                scanner.close();

                InferTreeWrapper inferTreeWrapper = new InferTreeWrapper();
                List<MutableTuple<Tree, Double>> distinctTrees = new ArrayList<>();
                inferTreeWrapper.summarizeData(treeList, null, distinctTrees);

                MajorityConsensusInference majorityConsensusInference = new MajorityConsensusInference();
                consensusTrees.add(majorityConsensusInference.inferSpeciesTree(distinctTrees, false, 50));
            }


        } catch (IOException e) {
            e.printStackTrace();
        }


        for(int i = 0 ; i < size ; i++) {
            for(int j = i + 1 ; j < size ; j++) {

                double distance = getDistance(consensusTrees.get(i), consensusTrees.get(j), false);
                distMatrix[i][j] = distMatrix[j][i] = distance;
            }
        }
    }

    String getSebastienTreeWithFourTaxa(List<List<MutableTuple<Tree, Double>>> originalGTs, Map<String,List<String>> species2alleles) {
        //List<List<MutableTuple>> gts = new ArrayList<List<MutableTuple>>();
        //for (List<MutableTuple> gtList : originalGTs) {
        //    gts.add(Arrays.asList(new MutableTuple(Trees.readTree(gtList.get(0).Item1.toString()), 1.0)));
        //}

        List<Tree> gts = new ArrayList<>();
        for (List<MutableTuple<Tree, Double>> gtList : originalGTs) {
            gts.add(Trees.readTree(gtList.get(0).Item1.toString()));
        }

        //Tree tmp = getGLASSTree(gts, null, false);


        LinkedList<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
        ArrayList<STITreeCluster> allClusters = new ArrayList<STITreeCluster>();

        List<String> temp = new LinkedList<String>();
        Tree t1 = null;
        for(Tree tr:gts){
            if(t1 == null){
                t1 = tr;
                for(TNode node:tr.getRoot().getLeaves()){
                    temp.add(node.getName());
                }
            }
            else{
                if(!Trees.leafSetsAgree(t1, tr)){
                    for(TNode node:tr.getRoot().getLeaves()){
                        if(!temp.contains(node.getName())){
                            temp.add(node.getName());
                        }
                    }
                }
            }
        }

        String[] taxa = new String[temp.size()];
        int index = 0;
        for(String taxon:temp){
            taxa[index++] = taxon;
        }

        Arrays.sort(taxa);

        for (String taxon: taxa){
            STITreeCluster cluster = new STITreeCluster(taxa);
            cluster.addLeaf(taxon);
            allClusters.add(cluster);
        }

        Map<STITreeCluster,ArrayList<Double>> taxonPairDM = new HashMap<>();
        Map<STITreeCluster,Double> taxonPairD = new HashMap<>();
        ArrayList<STITreeCluster<Double>> taxonPairs = new ArrayList<STITreeCluster<Double>>();

        for(int i=0;i<taxa.length;i++){
            for(int j=i+1;j<taxa.length;j++){
                STITreeCluster<Double> cl = new STITreeCluster<Double>(taxa);
                BitSet bs = new BitSet(taxa.length);
                bs.set(i);
                bs.set(j);
                cl.setCluster(bs);
                cl.setData(-1.0);
                taxonPairs.add(cl);
            }
        }

        for(STITreeCluster<Double> clwt : taxonPairs){
            STITreeCluster cl = new STITreeCluster(taxa);
            cl.setCluster(clwt.getCluster());
            taxonPairDM.put(cl, new ArrayList<>());
            taxonPairD.put(cl, clwt.getData());
        }

        for (Tree gt : gts) {
            Hashtable<TNode,Integer> id_lookup = new Hashtable<>();
            TNode[] node_lookup = new TNode[taxa.length];
            double [][] dist = Trees.getLeafDistanceMatrix(gt, id_lookup, node_lookup);

            for(STITreeCluster<Double> pair : taxonPairs) {
                STITreeCluster cl = new STITreeCluster(taxa);
                cl.setCluster(pair.getCluster());
                String [] p = pair.getClusterLeaves();
                taxonPairDM.get(cl).add(dist[id_lookup.get(gt.getNode(p[0]))][id_lookup.get(gt.getNode(p[1]))]);
            }

        }

        for(STITreeCluster<Double> pair : taxonPairs) {
            STITreeCluster cl = new STITreeCluster(taxa);
            cl.setCluster(pair.getCluster());
            Collections.sort(taxonPairDM.get(cl));
            ArrayList a = taxonPairDM.get(cl);
            double d = 0;
            if(a.size() % 2 == 0) {
                d = ((double)a.get(a.size() / 2 - 1) + (double)a.get(a.size() / 2)) / 2.0;
            }
            else {
                d = (double)a.get(a.size() / 2);
            }
            taxonPairD.put(cl, d);

        }

        MutableTree firstTree = new STITree<Object>();

        double minD = Double.MAX_VALUE;
        STITreeCluster minCluster = null;
        for(STITreeCluster<Double> pair : taxonPairs) {
            STITreeCluster cl = new STITreeCluster(taxa);
            cl.setCluster(pair.getCluster());
            if(minD > taxonPairD.get(cl)) {
                minD = taxonPairD.get(cl);
                minCluster = cl;
            }
        }

        firstTree.getRoot().createChild(minCluster.getClusterLeaves()[0]);
        firstTree.getRoot().createChild(minCluster.getClusterLeaves()[1]);

        for(String taxon : minCluster.getClusterLeaves()) {
            STITreeCluster cluster = new STITreeCluster(taxa);
            cluster.addLeaf(taxon);
            allClusters.remove(cluster);
        }

        String l_F1 = allClusters.get(0).getClusterLeaves()[0];
        String l_F2 = allClusters.get(1).getClusterLeaves()[0];
        String l_AUB = minCluster.getClusterLeaves()[0];
        allClusters.add(minCluster);

        STITreeCluster pair1 = new STITreeCluster(taxa);
        pair1.addLeaf(l_F1);
        pair1.addLeaf(l_AUB);

        STITreeCluster pair2 = new STITreeCluster(taxa);
        pair2.addLeaf(l_F2);
        pair2.addLeaf(l_AUB);

        STITreeCluster pair3 = new STITreeCluster(taxa);
        pair3.addLeaf(l_F2);
        pair3.addLeaf(l_F1);

        minD = taxonPairD.get(pair1);
        String thirdLeaf = l_F1;
        String fourthLeaf = l_F2;

        if(taxonPairD.get(pair1) > taxonPairD.get(pair2)) {
            minD = taxonPairD.get(pair2);
            thirdLeaf = l_F2;
            fourthLeaf = l_F1;
        }

        MutableTree secondTree = new STITree<Object>();

        if(minD > taxonPairD.get(pair3)) {
            MutableTree tree = new STITree<Object>();
            tree.getRoot().createChild(thirdLeaf);
            tree.getRoot().createChild(fourthLeaf);
            secondTree.getRoot().createChild(tree.getRoot());
            secondTree.getRoot().createChild(firstTree.getRoot());
        }
        else {
            MutableTree tree = new STITree<Object>();
            tree.getRoot().createChild(firstTree.getRoot());
            tree.getRoot().createChild(thirdLeaf);
            secondTree.getRoot().createChild(tree.getRoot());
            secondTree.getRoot().createChild(fourthLeaf);
        }

        return secondTree.toNewick();

    }

    public void testSebastienTree(String [] args) {
        System.out.println("Input File: " + args[0]);
        System.out.println("Output File: " + args[1]);

        String filename = args[0];
        List<List<MutableTuple<Tree, Double>>> gts = new ArrayList<List<MutableTuple<Tree, Double>>>();
        try {
            File file = new File(filename);
            Scanner scanner = new Scanner(file);

            while(scanner.hasNextLine()) {
                String line = scanner.nextLine();
                gts.add(Arrays.asList(new MutableTuple(Trees.readTree(line), 1.0)));
            }

            scanner.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }


        String result = getSebastienTreeWithFourTaxa(gts, null);

        try {
            PrintWriter out = new PrintWriter(args[1]);
            out.println(result);
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

    }

    public void checkResultsFor3Methods(String [] args) { //The result of 3 methods tested: ASTRAL, Dendroscope, Sebastien
        List<Double> results = new ArrayList<>();

        try {

            File file = new File("/Users/zhujiafan/Documents/Luay/Dendroscope/speciestree_gamma005.txt");
            Scanner scanner = new Scanner(file);

            int s = 0;
            for (int i = 0; i < 16 * 102; i++) {
                String str = scanner.nextLine();
                if(str.charAt(0) == '(') {
                    Tree tree = Trees.readTree(str);
                    Tree tree1 = Trees.readTree("(((b,c),a),d);");
                    Tree tree2 = Trees.readTree("(((b,c),d),a);");
                    SymmetricDifference symmetricDifference1 = new SymmetricDifference();
                    symmetricDifference1.computeDifference(tree1, tree, true);
                    double diff1 = symmetricDifference1.getWeightedAverage();

                    SymmetricDifference symmetricDifference2 = new SymmetricDifference();
                    symmetricDifference2.computeDifference(tree2, tree, true);
                    double diff2 = symmetricDifference2.getWeightedAverage();

                    if(Math.min(diff1, diff2) == 0) {
                        s++;
                    }
                }
                else if(str.charAt(0) == '0' || str.charAt(0) == '1') {
                    if(str.charAt(0) == '1')
                        s++;
                }
                else {
                    if(i % 102 == 0 && i > 2)
                        results.add(s / 100.0);
                    System.out.println(str);
                    s = 0;
                }
            }
            results.add(s / 100.0);

            scanner.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        for(int i = 0 ; i < 4 ; i++) {
            for(int j = 0 ; j < 4 ; j++) {
                System.out.print("\t");
                System.out.print(results.get(i * 4 + j));
            }
            System.out.println("");
        }


    }

    public void generateDataFor3Methods(String [] args){
        int numTestSet = 0;
        int simulateSize = 2000;
        double targetEdgeInheriProb = 1 - 0.05;
        Func2<Network, Integer, List> simulator = getSimulator(null);

        for(int i = 0 ; i < numTestSet ; i++) {
            new File("./data/" + i).mkdir();
        }

        double branchLengths[] = new double[]{0.1, 0.2, 0.5, 1.0};
        //List<Double> branchLengths = new ArrayList<>();
        //for(int k = 1 ; k <= 30 ; k++)
        //branchLengths.add(k / 1000.0);

        for(double _targetEdgeBrlen : branchLengths) {
            Network<Object> trueNetwork;
            trueNetwork = Networks.readNetwork("((((b:1000.0,c:1000.0)I4:" + _targetEdgeBrlen + ")I3#H1:0.0::" + (1 - targetEdgeInheriProb) + ",a:" + (1000 + _targetEdgeBrlen) + ")I1:0.1,(I3#H1:0.0::" + targetEdgeInheriProb + ",d:" + (1000 + _targetEdgeBrlen) + ")I2:0.1)I0;");

            for (int num = 0; num < numTestSet; num++) {

                List<List<MutableTuple<Tree, Double>>> simulatedGTs = simulator.execute(trueNetwork, 25);
                Collections.shuffle(simulatedGTs);
                int sizes[] = new int[]{25, 50, 100, 200};

                for (int size : sizes)
                    try {
                        if(size != 25)
                            simulatedGTs.addAll(simulator.execute(trueNetwork, size / 2));

                        String GTsWithBLFilename = "data/" + num + "/Special4_" + ((int) Math.round(_targetEdgeBrlen * 10)) / 10.0 + "_" + size + "GTsWithBL";
                        PrintWriter GTsWithBLOut = new PrintWriter(GTsWithBLFilename);

                        String GTsWithoutBLFilename = "data/" + num + "/Special4_" + ((int) Math.round(_targetEdgeBrlen * 10)) / 10.0 + "_" + size + "GTsWithoutBL";
                        PrintWriter GTsWithoutBLOut = new PrintWriter(GTsWithoutBLFilename);

                        for (int j = 0; j < size; j++) {
                            Tree tree = (Tree) simulatedGTs.get(j).get(0).Item1;
                            GTsWithBLOut.println(tree);

                            Tree treeWithoutBL =  Trees.readTree(tree.toNewick());
                            for (TNode node : treeWithoutBL.postTraverse()) {
                                node.setParentDistance(TNode.NO_DISTANCE);
                            }

                            GTsWithoutBLOut.println(treeWithoutBL);
                        }

                        GTsWithBLOut.close();
                        GTsWithoutBLOut.close();
                    } catch (IOException ioe) {
                        ioe.printStackTrace();
                    }
            }
        }
    }

    public class SimulationConfiguration {
        public String _clusteringMethod;
        public String _inferringMethod;
        public String _distanceMethod;
        public int[] _sizes;
        public int _fixedClusterNumberTest;
        public int _lockedClusterNumber;
        public boolean _saveDistanceMatrixThenExit;
        public boolean _updateDistMatrixAfterMDS;
        public int _MDSdim;

        public String intArrayToString(int [] a) {
            StringBuilder sb = new StringBuilder();
            if(a != null) {
                for (int i = 0; i < a.length; i++) {
                    sb.append(a[i]);
                    sb.append(' ');
                }
            }
            return sb.toString();
        }

        @Override
        public String toString() {
            StringBuffer str = new StringBuffer("");
            str.append("_clusteringMethod " + _clusteringMethod + "\n");
            str.append("_inferringMethod " + _inferringMethod + "\n");
            str.append("_distanceMethod " + _distanceMethod + "\n");
            str.append("_sizes " + intArrayToString(_sizes) + "\n");
            str.append("_fixedClusterNumberTest " + _fixedClusterNumberTest + "\n");
            str.append("_lockedClusterNumber " + _lockedClusterNumber + "\n");
            str.append("_saveDistanceMatrixThenExit " + _saveDistanceMatrixThenExit + "\n");
            str.append("_updateDistMatrixAfterMDS " + _updateDistMatrixAfterMDS + "\n");
            str.append("_MDSdim " + _MDSdim + "\n");
            return str.toString();
        }
    }

    public class NewNetworkProcessor extends RecursiveTask<String> {
        private int _number = 0;
        private String _path;
        private boolean _useFile;
        private SimulationConfiguration _config;
        private InferNetworkClustering _inferNetworkClustering;
        private ParentalTreeOperation _parentalTreeOperation;
        private Random _random;
        private Func2<Network, Integer, List> _simulator;
        private InferTreeWrapper _inferTreeWrapper;

        public NewNetworkProcessor(int number, String path, SimulationConfiguration config) {
            _number = number;
            _path = path;
            if(path == null) {
                _useFile = false;
            }
            else {
                _useFile = true;
            }
            _config = config;
            _inferNetworkClustering = new InferNetworkClustering();
            _parentalTreeOperation = new ParentalTreeOperation();
            _random = new Random();
            _simulator = getSimulator(null);
            _inferTreeWrapper = new InferTreeWrapper();
        }

        public Network<Object> getNewNetwork() {
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
            double blK = 0.2;
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

        protected String compute() {
            StringBuilder output = new StringBuilder();
            List<Map<String, String>> result = new ArrayList<>();
            Network<Object> trueNetwork = null;


            if(_useFile) {
                try {
                    File file = new File(_path + "/network.txt");
                    Scanner scanner = new Scanner(file);

                    trueNetwork = Networks.readNetwork(scanner.nextLine());

                    scanner.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            } else {
                trueNetwork = getNewNetwork();
            }


            List<Tree> parentalTrees = _parentalTreeOperation.getParentalTrees(trueNetwork);

            List<MutableTuple<Tree, Integer>> expectations = new ArrayList<>();
            for(Tree parentalTree : parentalTrees) {
                expectations.add(new MutableTuple(Trees.readTree(parentalTree.toNewick()), 0));
            }
            //expectations.add(new MutableTuple(Trees.readTree("(((((B:1.0,A:1.0):1.0,C:1.0):1.0,G:1.0):1.0,H:1.0):1.0,((((F:1.0,E:1.0):1.0,D:1.0):1.0,J:1.0):1.0,I:1.0):1.0);"), 0));
            //expectations.add(new MutableTuple(Trees.readTree("(((((B:1.0,A:1.0):1.0,D:1.0):1.0,G:1.0):1.0,H:1.0):1.0,((((F:1.0,E:1.0):1.0,C:1.0):1.0,J:1.0):1.0,I:1.0):1.0);"), 0));
            //expectations.add(new MutableTuple(Trees.readTree("(((((D:1.0,C:1.0):1.0,(B:1.0,A:1.0):1.0):1.0,G:1.0):1.0,H:1.0):1.0,(((F:1.0,E:1.0):1.0,J:1.0):1.0,I:1.0):1.0);"), 0));
            //expectations.add(new MutableTuple(Trees.readTree("(((((F:1.0,E:1.0):1.0,(D:1.0,C:1.0):1.0):1.0,J:1.0):1.0,I:1.0):1.0,(((B:1.0,A:1.0):1.0,G:1.0):1.0,H:1.0):1.0);"), 0));

            List<List<MutableTuple<Tree, Double>>> expectedGTs = new ArrayList<>();
            for(int i = 0 ; i < expectations.size() ; i++){
                expectedGTs.add(Arrays.asList(new MutableTuple(expectations.get(i).Item1, 1.0)));
            }


            int sizes[] = {50, 250, 500, 1000};
            if(_config._sizes != null)
                sizes = _config._sizes;

            List<List<MutableTuple<Tree, Double>>> simulatedGTs = null;

            for(int sizeIndex = 0 ; sizeIndex < sizes.length ; sizeIndex++) {
                int size = sizes[sizeIndex];

                if(_useFile) {
                    simulatedGTs = loadGTs(_path + "/" + size + ".txt");
                } else {
                    if(sizeIndex == 0) {
                        simulatedGTs = _simulator.execute(trueNetwork, size);
                    } else {
                        simulatedGTs.addAll(_simulator.execute(trueNetwork, size - sizes[sizeIndex - 1]));
                    }
                }

                for(int i = 0 ; i < expectations.size() ; i++){
                    expectations.get(i).Item2 = 0;
                }

                double distMatrix[][] = new double[size][size];
                if (_config._distanceMethod.equals("File"))
                    loadMatrix(_path + "/500.mat", size, distMatrix);
                else if(_config._distanceMethod.equals("BSD")) {
                    getDistancePHYLIP(simulatedGTs, distMatrix);
                }
                else if(_config._distanceMethod.equals("RF")){
                    for (int i = 0; i < size; i++)
                        for (int j = i + 1; j < size; j++)
                            distMatrix[i][j] = distMatrix[j][i] = getDistance(simulatedGTs.get(i).get(0).Item1, simulatedGTs.get(j).get(0).Item1);

                } else if(_config._distanceMethod.equals("BestTree")) {
                    getDistanceMatrixBestTree(simulatedGTs, _path, distMatrix);
                } else if(_config._distanceMethod.equals("Bootstrap")) {
                    getDistanceMatrixBootstrap(simulatedGTs, _path, distMatrix);
                } else if(_config._distanceMethod.equals("BootstrapConsensus")) {
                    getDistanceMatrixBootstrapConsensus(simulatedGTs, _path, distMatrix);
                } else {
                    throw new RuntimeException("Unknown distance method: " + _config._distanceMethod);
                }

                if(_config._saveDistanceMatrixThenExit) {
                    try {
                        PrintWriter out = new PrintWriter(new FileOutputStream(new File(_path + "/" + size + ".mat")));
                        out.println(size);
                        for (int i = 0; i < size; i++) {
                            for (int j = 0; j < size; j++)
                                out.print("\t" + distMatrix[i][j]);
                            out.println();
                        }
                        out.close();
                    } catch (IOException e){
                        e.printStackTrace();
                    }
                    return null;
                }

                int dim = _config._MDSdim;


                double coordinates[][] = new double[size][dim];
                _inferNetworkClustering.getMDScoordinates(size, dim, distMatrix, coordinates);
                if(_config._updateDistMatrixAfterMDS) {
                    for (int i = 0; i < size; i++)
                        for (int j = i + 1; j < size; j++)
                            distMatrix[i][j] = distMatrix[j][i] = _inferNetworkClustering.getEuclideanDistance(dim, coordinates[i], coordinates[j]);
                }


                List<Double> silhouetteIndices = new ArrayList<>();
                double maxSilhouetteIndex = Double.MIN_VALUE;
                int bestNumOfClusters = 2;
                int numOfClusters = 2;
                int numOfClustersLimit = 100;
                if(_config._fixedClusterNumberTest > 0)
                    numOfClustersLimit = _config._fixedClusterNumberTest;
                if(_config._lockedClusterNumber > 0) {
                    numOfClusters = _config._lockedClusterNumber;
                    bestNumOfClusters = _config._lockedClusterNumber;
                    numOfClustersLimit = 0;
                }


                while (numOfClusters <= numOfClustersLimit) {
                    List<List<Integer>> cluster = new ArrayList<>();


                    if(_config._clusteringMethod.equals("MDS+KMeans")) {
                        _inferNetworkClustering.KMeansClustering(cluster, coordinates, size, numOfClusters, dim);
                    }
                    else if(_config._clusteringMethod.equals("MDC+KMeans")) {
                        _inferNetworkClustering.KMeansClustering_MDC(cluster, simulatedGTs, size, numOfClusters);
                    }
                    else if(_config._clusteringMethod.equals("MDS+KMedoids")) {
                        double newdistMatrix[][] = new double[size][size];
                        for (int i = 0; i < size; i++)
                            for (int j = i + 1; j < size; j++)
                                newdistMatrix[i][j] = newdistMatrix[j][i] = _inferNetworkClustering.getEuclideanDistance(dim, coordinates[i], coordinates[j]);

                        _inferNetworkClustering.KMedoidsClustering(cluster, newdistMatrix, size, numOfClusters);
                    }
                    else
                        throw new RuntimeException("Unknown clustering method: " + _config._clusteringMethod);
                    /**/
                    double score = _inferNetworkClustering.getSilhouetteIndex(cluster, distMatrix, size, numOfClusters);
                    //System.out.println("numOfCluster = " + numOfClusters);
                    //System.out.println("SilhouetteIndex = " + score);
                    silhouetteIndices.add(score);

                    if (maxSilhouetteIndex < score) {
                        maxSilhouetteIndex = score;
                        bestNumOfClusters = numOfClusters;
                    } else if(_config._fixedClusterNumberTest <= 0){
                        bestNumOfClusters = numOfClusters - 1;
                        break;
                    }

                    numOfClusters++;


                }


                numOfClusters = bestNumOfClusters;
                List<List<Integer>> cluster = new ArrayList<>();

                List<String> results = new ArrayList<>();
                List<List<MutableTuple<Tree, Double>>> resultGTs = new ArrayList<>();
                List<Tree> inferredTrees = new ArrayList<>();


                if(_config._clusteringMethod.equals("MDS+KMeans")) {
                    _inferNetworkClustering.KMeansClustering(cluster, coordinates, size, numOfClusters, dim);
                }
                else if(_config._clusteringMethod.equals("MDC+KMeans")) {
                    _inferNetworkClustering.KMeansClustering_MDC(cluster, simulatedGTs, size, numOfClusters);
                }
                else if(_config._clusteringMethod.equals("MDS+KMedoids")) {
                    double newdistMatrix[][] = new double[size][size];
                    for (int i = 0; i < size; i++)
                        for (int j = i + 1; j < size; j++)
                            newdistMatrix[i][j] = newdistMatrix[j][i] = _inferNetworkClustering.getEuclideanDistance(dim, coordinates[i], coordinates[j]);

                    _inferNetworkClustering.KMedoidsClustering(cluster, newdistMatrix, size, numOfClusters);
                }
                else
                    throw new RuntimeException("Unknown clustering method: " + _config._clusteringMethod);

                for (int k = 0; k < numOfClusters; k++) {
                    List<Tree> gts = new ArrayList<>();
                    for (Integer i : cluster.get(k))
                        gts.add(simulatedGTs.get(i).get(0).Item1);
                    System.out.println("Cluster size = " + gts.size());
                    if (gts.size() == 0)
                        continue;
                    String speciesTree = _inferTreeWrapper.inferTreeByMethod(gts, null, new String(_config._inferringMethod));
                    System.out.println(speciesTree);
                    results.add(speciesTree);
                    inferredTrees.add(Trees.readTree(speciesTree));
                    resultGTs.add(Arrays.asList(new MutableTuple<Tree, Double>(Trees.readTree(speciesTree), 1.0)));
                    for (MutableTuple<Tree, Integer> t : expectations) {
                        if (getDistance(t.Item1, Trees.readTree(speciesTree)) < 1e-6) {
                            ++t.Item2;
                        }
                    }
                }

                double bestMatch2 = getBestMatchDistanceBetweenTrees(expectedGTs, resultGTs);
                double bestMatch = 0.0;//Networks.computeTreeDistance(inferredTrees, parentalTrees)[2];

                int count = 0;
                for (MutableTuple<Tree, Integer> t : expectations) {
                    if (t.Item2 > 0) {
                        ++count;
                    }
                }

                System.out.println("Test #" + _number);
                System.out.println("Size of GTs = " + size);
                System.out.println("Silhouette Index = " + maxSilhouetteIndex);
                System.out.println("Number of clusters = " + bestNumOfClusters);
                System.out.println("Number of correct parental trees = " + count);
                System.out.println("Best match error = " + bestMatch);
                System.out.println("Best match error2 = " + bestMatch2);



                output.append("Size of gene trees:\n" + size + "\n");
                output.append("Silhouette Index:\n" + maxSilhouetteIndex + "\n");
                output.append("Number of clusters:\n" + bestNumOfClusters + "\n");
                output.append("Number of correct parental trees:\n" + count + "\n");
                output.append("Best match error:\n" + bestMatch + "\n");
                output.append("Best match error2:\n" + bestMatch2 + "\n");
                output.append("Test #\n" + _number + "\n");
                output.append("Path:\n" + _path + "\n");
                output.append("Silhouette Indices: ");
                for(Double x : silhouetteIndices)
                    output.append(" " + x);
                output.append("\n");
                for (int k = 0; k < numOfClusters; k++) {
                    output.append("Cluster size = " + cluster.get(k).size() + "\n");
                    output.append(results.get(k) + "\n");
                }
                output.append("\n");

            }

            return output.toString();
        }

    }

    public void testNewNetworkMultithread(String[] args) {
        boolean useFile = true;
        boolean singleThread = false;
        ForkJoinPool pool = new ForkJoinPool();
        SimulationConfiguration config = new SimulationConfiguration();
        config._clusteringMethod = new String("MDS+KMeans");
        config._inferringMethod = new String("MDC");
        config._distanceMethod = "BestTree";
        config._sizes = new int[]{50, 250, 500, 1000};
        config._fixedClusterNumberTest = -1;
        config._lockedClusterNumber = -1;
        config._saveDistanceMatrixThenExit = false;
        config._updateDistMatrixAfterMDS = false;
        config._MDSdim = 3;

        System.out.println(config.toString());

        System.out.println("Clustering method:" + config._clusteringMethod);
        System.out.println("Inferring method: " + config._inferringMethod);

        if(singleThread)
            pool = new ForkJoinPool(1);

        String path = "newdata/";
        List<NewNetworkProcessor> processors = new ArrayList<>();
        int number = -1;

        int curPart = 0;
        int totalParts = 0;

        if(args.length == 2) {
            curPart = Integer.parseInt(args[0]);
            totalParts = Integer.parseInt(args[1]);
        }

        if(useFile) {
            List<String> paths = new ArrayList<>();
            for (File file : new File(path).listFiles()) {
                if (file.isDirectory()) {
                    paths.add(file.getPath());
                    number++;
                }
            }
            Collections.sort(paths);
            int totalSets = paths.size();

            if(totalParts > 0) {
                int begin = totalSets / totalParts * curPart;
                int end = totalSets / totalParts * (curPart + 1);
                if(curPart == totalParts - 1)
                    end = totalSets;
                paths = paths.subList(begin, end);
            }

            for(number = 0 ; number < paths.size() ; number++) {
                //if(number != 6) continue;
                System.out.println(paths.get(number));
                NewNetworkProcessor p = new NewNetworkProcessor(number, paths.get(number), config);
                pool.execute(p);
                processors.add(p);
            }

        } else {
            for(number = 0 ; number < 300 ; number++) {
                NewNetworkProcessor p = new NewNetworkProcessor(number, null, config);
                pool.execute(p);
                processors.add(p);
            }
        }

        long seconds = 0;
        do
        {
            System.out.println("Time Elapsed: " + seconds + " seconds.");
            System.out.printf("******************************************\n");
            System.out.printf("Main: Parallelism: %d\n", pool.getParallelism());
            System.out.printf("Main: Active Threads: %d\n", pool.getActiveThreadCount());
            System.out.printf("Main: Task Count: %d\n", pool.getQueuedTaskCount());
            System.out.printf("Main: Steal Count: %d\n", pool.getStealCount());
            System.out.printf("******************************************\n");
            seconds += 10;
            try
            {
                TimeUnit.SECONDS.sleep(10);
            } catch (InterruptedException e)
            {
                e.printStackTrace();
            }
        } while (!pool.isQuiescent());

        List<String> results = new ArrayList<>();

        for (NewNetworkProcessor p : processors){
            results.add(p.join());
        }

        try {
            String logFilename = "log.txt";
            if(totalParts > 0) {
                logFilename = "log_" + curPart + ".txt";
            }
            PrintWriter out = new PrintWriter(new FileOutputStream(new File(logFilename)));
            for(String result : results) {
                if(result == null) continue;
                out.println(result);
                out.println();
            }
            out.close();
        }catch (IOException ioe) {
            ioe.printStackTrace();
        }


    }


    public void testHungarianMatching(){
        double [][] costMatrix = {{82,83,69,92}, {77,37,49,92}, {11,69,5,86}, {8,9,98,23}};
        int workers = costMatrix.length;
        int jobs = costMatrix[0].length;
        HungarianMatching HungarianMatching = new HungarianMatching(costMatrix);
        int [] assignment = HungarianMatching.compute();
        System.out.println("Assignments:");
        double cost = 0;

        //Should be 2 1 0 3
        for(int i = 0 ; i < workers ; i++) {
            System.out.println("\t" + assignment[i]);
            if(assignment[i] != -1) {
                cost += costMatrix[i][assignment[i]];
            }
        }
        //Should be 140
        System.out.println("Cost: " + cost);

        BipartiteGraph bg = new BipartiteGraph(costMatrix.length, costMatrix[0].length);
        for(int i = 0 ; i < workers ; i++) {
            for(int j = 0 ; j < jobs ; j++) {
                bg.addEdge(i, j, costMatrix[i][j]);
            }
        }
        System.out.println("Cost: " + bg.getMinEdgeCoverWeight());

    }

    public void generateData() {
        //int sizes[] = {50, 250, 500, 1000};
        int sizes[] = {50, 250, 500, 1000, 20000};
        DataGenerator dataGenerator = new DataGenerator(10, sizes, 30);
        dataGenerator.generateFiles("R2_1_2long_0.5_0.5");
        //Network<Object> network = dataGenerator.getNewNetworkR2_2();
        //List<Tree> trees = getParentalTrees(network);
        //System.out.println(trees);
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

    public void printGeneTreeProbWithFourTaxa(Network<Object> trueNetwork) {
        List<Double> branchLengths = new ArrayList<>();
        GeneTreeProbability geneTreeProbability = new GeneTreeProbability();
        List<Double> gtprob = new ArrayList<>();
        List<Tree> gts = new ArrayList<>();
        gts.add(Trees.readTree("(((b, c), a), d);"));
        gts.add(Trees.readTree("(((b, c), d), a);"));
        gts.add(Trees.readTree("((a, b), (c, d));"));
        gts.add(Trees.readTree("((a, c), (b, d));"));
        gts.add(Trees.readTree("(((a, b), c), d);"));
        gts.add(Trees.readTree("(((a, c), b), d);"));
        gts.add(Trees.readTree("(a, (b, (c, d)));"));
        gts.add(Trees.readTree("(((b, d), c), a);"));
        gts.add(Trees.readTree("((a, d), (b, c));"));
        gts.add(Trees.readTree("(((a, b), d), c);"));
        gts.add(Trees.readTree("(b, (a, (c, d)));"));
        gts.add(Trees.readTree("(((a, d), b), c);"));
        gts.add(Trees.readTree("(((b, d), a), c);"));
        gts.add(Trees.readTree("(((a, c), d), b);"));
        gts.add(Trees.readTree("(((a, d), c), b);"));


        gtprob = geneTreeProbability.calculateGTDistribution((Network)trueNetwork, gts, null, false);
        //System.out.println(gts);
        for(double x : gtprob){
            System.out.print('\t');
            System.out.print(x);
        }
        System.out.println("");
        System.out.println("");
    }

    public void testTable1(){
        double gamma = 1 - 0.05;

        double y = 0.1;
        double x = 10000;

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork("((((b:1000.0,c:1000.0)I4:" + y + ")I3#H1:0.0::" + (1 - gamma) + ",a:" + (1000 + y) + ")I1:" + x + ",(I3#H1:0.0::" + gamma + ",d:" + (1000 + y) + ")I2:"+ x +")I0;");

        printGeneTreeProbWithFourTaxa(trueNetwork);
    }

    public static void main(String[] args) {

        /*SymmetricDifference sd = new SymmetricDifference();
        Tree t1 = Trees.readTree("(((b, c), a), d);");
        Tree t2 = Trees.readTree("(((b, c), d), a);");
        sd.computeDifference(t1, t2, true);
        System.out.println(sd.getUnweightedAverage());
        System.out.println(sd.getWeightedAverage());
        System.out.println(sd.getNumInternalEdges1());
        System.out.println(sd.getNumInternalEdges2());
*/




        ParentalTreeOperation parentalTreeOperation = new ParentalTreeOperation();
        parentalTreeOperation.test();

        if(true)
            return;

        RECOMB_CG_16_JZ recomb_cg_16_jz = new RECOMB_CG_16_JZ();
        recomb_cg_16_jz.testTable1();


        Network<Object> net = Networks.readNetwork("((A,B),C);");
        System.out.println(net.toString());

        List<Tree> trees = new ArrayList<>();
        trees.add(Trees.readTree("(((((((8,10),7))),(((((3)),1),6),4)),(5,2)),9);"));
        trees.add(Trees.readTree("(((((3)),(((((((8,10),7))),1),6),4)),(5,2)),9);"));
        trees.add(Trees.readTree("(((((((8)))),(((((3,((10),7))),1),6),4)),(5,2)),9);"));

        //trees.add(Trees.readTree("((((a,b),e),c),d);"));

        for(int i = 0 ; i < trees.size() ; i++) {
            InferNetworkFromParentalTrees inferNetworkFromParentalTrees = new InferNetworkFromParentalTrees();
            Network<Object> tmp = Networks.readNetwork(trees.get(i).toNewick());
            inferNetworkFromParentalTrees.removeBinaryNodes(tmp);
            trees.set(i, Trees.readTree(tmp.toString()));

        }
        SteelWarnowMAST steelWarnowMAST = new SteelWarnowMAST();
        //Tree mast = steelWarnowMAST.computeRMAST(trees);
        Tree mast = steelWarnowMAST.computeRMAST(trees.get(0), trees.get(1));
        System.out.println(mast);


    }
}
