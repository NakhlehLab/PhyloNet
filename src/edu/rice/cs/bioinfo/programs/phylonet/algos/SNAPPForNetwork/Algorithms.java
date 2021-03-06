package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;


import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.Splitting;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.FelsensteinAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import jeigen.DenseMatrix;
import org.apache.commons.math3.util.ArithmeticUtils;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * This class holds all the various algorithms from http://mbe.oxfordjournals.org/content/29/8/1917.
 */
public class Algorithms
{
    private static final boolean PRINT_DETAILS = false;
    public static boolean CORRECTION_AT_LEAVES = false;
    public static boolean HAS_DOMINANT_MARKER = false;
    public static final boolean SWITCH_FASTER_BIALLILE = true;
    public static final boolean SWITCH_EXP_APPROX = true;
    public static boolean SWITCH_APPROX_SPLIT = false;
    public static boolean SWITCH_FMATRIX_CACHE = false;
    public static List<List<Integer>> targetSplittingIndices = null;

    private static int[] mergeTwoSplittingIndices(int[] index1, int[] index2){
        for(int i=0; i<index1.length; i++) {
            if(index1[i] != index2[i] && index1[i] != 0 && index2[i] != 0)
                return null;
        }
        int[] newIndex = new int[index1.length];
        for(int i=0; i<index1.length; i++){
            if(index1[i] == index2[i]){
                newIndex[i] = index1[i];
            }
            else if(index1[i] == 0){
                newIndex[i] = index2[i];
            }
            else if(index2[i] == 0){
                newIndex[i] = index1[i];
            }
            else{
                return null;
            }
        }

        return newIndex;
    }




    /**
     * This transforms a nucleotide observation map of strings to characters to a "nucleotideIndex" map of strings to ints.
     *
     * @param obs The nucleotide observations.
     * @return The nucleotideIndex observations.
     */
    private static Map<String, R> getNucleotideIndexMap(NucleotideObservation obs, Map<String, String> allele2species)
    {
        Map<String, int[]> colorMap = new HashMap<String, int[]>();
        for (String allele : obs.getAlleles())
        {
            String species = allele;
            if(allele2species != null){
                species = allele2species.get(allele);
            }
            int[] counter = colorMap.get(species);
            if(counter==null){
                counter = new int[4];
                colorMap.put(species, counter);
            }
            counter[alleleToColor(obs.getObservationForAllele(allele))]++;

        }
        Map<String, R> leafNode2R = new HashMap<String, R>();
        for(Map.Entry<String, int[]> entry: colorMap.entrySet()){
            int[] Rvalue = new int[R.dims];
            int n = 0;
            for(int i = 0; i< R.dims; i++){
                Rvalue[i] = entry.getValue()[i];
                n += entry.getValue()[i];
            }
            n += entry.getValue()[R.dims];
            leafNode2R.put(entry.getKey(), new R(n, Rvalue));
        }
        return leafNode2R;
    }

    /**
     * Converts a nucleotide to a "nucleotideIndex".
     *
     * @param character The nucleotide to convert.
     * @return The "nucleotideIndex".
     */
    private static int alleleToColor(char character)
    {
        if (R.getNumberOfTypes() == 2)
        {
            switch (character)
            {
                case '0':
                    return 0;
                case '1':
                    return 1;
                default:
                    throw new RuntimeException("Invalid nucleotide only 0 or 1 allowed. Color was: " + character);
            }
        } else
        {
            int number = FelsensteinAlgorithm.nucleotides.indexOf(character);

            if (number == -1)
                throw new RuntimeException("Bad nucleotide was " + character);
            return number;
        }
    }

    /**
     * This is the only public method of this file.
     * It calculates the probabilty of observing a species network given some observations and a MatrixQ.
     * In other words, it calculates P(obs | speciesTree, Q)
     *
     * @param speciesNetwork The species network.
     * @param pattern         The nucleotide observations.
     * @param Q           The transition matrix for SNAPP.
     * @return The probability.
     */
    @SuppressWarnings("unchecked")
    public static double getProbabilityObservationGivenNetwork(Network<SNAPPData[]> speciesNetwork, RPattern pattern, QParameters Q, int siteID)
    {
        //Map<String, R> nucleotideIndexMap = getNucleotideIndexMap(obs, allele2species);
        Map<String, R> nucleotideIndexMap = pattern.getPattern();
        Set<String> articulationNodes = new HashSet<String>();
        for(Object node: Networks.getLowestArticulationNodes(cloneNetwork(speciesNetwork))){
            articulationNodes.add(((NetNode)node).getName());
        }
        int numReticulations = speciesNetwork.getReticulationCount();
        int reticulationID = 0;
        for (NetNode<SNAPPData[]> node : Networks.postTraversal(speciesNetwork))
        {
            boolean isArticulation = articulationNodes.contains(node.getName());
            processNode(Q, nucleotideIndexMap, node, isArticulation, numReticulations, reticulationID, siteID);
            if(node.isNetworkNode()){
                reticulationID++;
            }
        }

        return getProbabilityOfNetwork(speciesNetwork, Q, siteID);
    }

    public static double getProbabilityObservationGivenNetwork(Network<SNAPPData[]> speciesNetwork, Map<String,String> allele2species, NucleotideObservation obs, QParameters Q, int siteID)
    {
        Map<String, R> nucleotideIndexMap = getNucleotideIndexMap(obs, allele2species);
        Set<String> articulationNodes = new HashSet<String>();
        for(Object node: Networks.getLowestArticulationNodes(cloneNetwork(speciesNetwork))){
            articulationNodes.add(((NetNode)node).getName());
        }
        int numReticulations = speciesNetwork.getReticulationCount();
        int reticulationID = 0;
        for (NetNode<SNAPPData[]> node : Networks.postTraversal(speciesNetwork))
        {
            boolean isArticulation = articulationNodes.contains(node.getName());
            processNode(Q, nucleotideIndexMap, node, isArticulation, numReticulations, reticulationID, siteID);
            if(node.isNetworkNode()){
                reticulationID++;
            }
        }

        return getProbabilityOfNetwork(speciesNetwork, Q, siteID);
    }


    public static double updateProbabilityObservationGivenNetwork(Network<SNAPPData[]> speciesNetwork, QParameters Q, int siteID, RPattern pattern, Map<String, String> allele2species, Map<NetNode, Boolean> node2update)
    {
        //Map<String, R> nucleotideIndexMap = (obs != null) ? getNucleotideIndexMap(obs, allele2species) : null;
        Map<String, R> nucleotideIndexMap = pattern.getPattern();
        Set<String> articulationNodes = new HashSet<String>();
        for(Object node: Networks.getLowestArticulationNodes(cloneNetwork(speciesNetwork))){
            articulationNodes.add(((NetNode)node).getName());
        }
        int numReticulations = speciesNetwork.getReticulationCount();
        int reticulationID = 0;
        for (NetNode<SNAPPData[]> node : Networks.postTraversal(speciesNetwork))
        {
            /*
            if(node2update.containsKey(node)) {
                boolean isArticulation = articulationNodes.contains(node.getName());
                updateNode(Q, node, isArticulation, numReticulations, reticulationID, siteID, node2update.get(node));
            }
            */

            boolean isArticulation = articulationNodes.contains(node.getName());
            updateNode(speciesNetwork, Q, node, isArticulation, numReticulations, reticulationID, siteID, nucleotideIndexMap, node2update.get(node));


            if(node.isNetworkNode()){
                reticulationID++;
            }
        }

        return getProbabilityOfNetwork(speciesNetwork, Q, siteID);
    }

    public static double updateProbabilityObservationGivenNetwork(Network<SNAPPData[]> speciesNetwork, QParameters Q, int siteID, NucleotideObservation obs, Map<String, String> allele2species, Map<NetNode, Boolean> node2update)
    {
        Map<String, R> nucleotideIndexMap = (obs != null) ? getNucleotideIndexMap(obs, allele2species) : null;
        Set<String> articulationNodes = new HashSet<String>();
        for(Object node: Networks.getLowestArticulationNodes(cloneNetwork(speciesNetwork))){
            articulationNodes.add(((NetNode)node).getName());
        }
        int numReticulations = speciesNetwork.getReticulationCount();
        int reticulationID = 0;
        for (NetNode<SNAPPData[]> node : Networks.postTraversal(speciesNetwork))
        {
            /*
            if(node2update.containsKey(node)) {
                boolean isArticulation = articulationNodes.contains(node.getName());
                updateNode(Q, node, isArticulation, numReticulations, reticulationID, siteID, node2update.get(node));
            }
            */

            boolean isArticulation = articulationNodes.contains(node.getName());
            updateNode(speciesNetwork, Q, node, isArticulation, numReticulations, reticulationID, siteID, nucleotideIndexMap, node2update.get(node));


            if(node.isNetworkNode()){
                reticulationID++;
            }
        }

        return getProbabilityOfNetwork(speciesNetwork, Q, siteID);
    }



    /**
     * This method calculates the probabilty of data given a processed species network.
     *
     * @param speciesNetwork An already processed species network.
     * @param Q           A SNAPP transition matrix.
     * @return The probability
     */
    protected static double getProbabilityOfNetwork(Network<SNAPPData[]> speciesNetwork, QParameters Q, int siteID)
    {
        FMatrix rootFBot = speciesNetwork.getRoot().getData()[siteID].getFBottoms(null).iterator().next().Item1;

        double theta = speciesNetwork.getRoot().getRootPopSize();
        MatrixQ matQ;
        if(Q._gTheta != null) {
            matQ = Q._gMatrix;
        } else {
            matQ = new MatrixQ(Q._rModel, Q._M, theta);
        }
        DenseMatrix eq = matQ.getEquilibrium();
        double sum = 0;
        for (int n = 1; n <= rootFBot.mx; n++)
        {
            for (R r : R.loopOver(n))
            {
                //if(r.getIndex() < eq.rows)
                sum += rootFBot.get(r) * eq.get(r.getIndex(), 0);
            }
        }
        if(sum == 0.0)
            sum = 1e-6;
        return sum;
    }

    /**
     * Process a node for the snapp algorithm.
     *
     * @param Q                  The transition matrix for SNAPP.
     * @param nucleotideIndexMap The observations encoded as indices.
     * @param node               The node to process.
     */
    private static void processNode(QParameters Q, Map<String, R> nucleotideIndexMap, NetNode<SNAPPData[]> node, boolean isArticulation, int numReticulations, int reticulationID, int siteID)
    {
        if(PRINT_DETAILS){
            System.out.println("\nNode " + node.getName());
        }
        SNAPPData data = new SNAPPData();
        node.getData()[siteID] = data;

        if (node.isLeaf())
        {
            processLeaf(nucleotideIndexMap, node, data, numReticulations);
        } else if(node.isTreeNode())
        {
            processInternalTreeNode(node, data, isArticulation, siteID);
        } else if(node.isNetworkNode()){
            processNetworkNode(node, data, numReticulations, reticulationID, siteID);
            //examineBranchOfNode(Q, node, data);
        }
        if(PRINT_DETAILS){
            System.out.println("FBottoms:");
            data.printFMatrixs("b");
        }

        processTopBranchOfNode(Q, node, data);
        if(PRINT_DETAILS){
            System.out.println("FTop:");
            data.printFMatrixs("t");
        }
    }


    /**
     * Process a node for the snapp algorithm.
     *
     * @param Q                  The transition matrix for SNAPP.
     * @param node               The node to process.
     */
    private static void updateNode(Network<SNAPPData[]> speciesNetwork, QParameters Q, NetNode<SNAPPData[]> node, boolean isArticulation, int numReticulation, int reticulationID, int siteID, Map<String, R> nucleotideIndexMap, Boolean constructFBottom)
    {
        if(PRINT_DETAILS){
            System.out.println("\nNode " + node.getName());
        }
        SNAPPData data = node.getData()[siteID];
        if(constructFBottom==null){
            if (PRINT_DETAILS) {

                System.out.println("FBottoms:");
                data.printFMatrixs("b");
                System.out.println("FTop:");
                data.printFMatrixs("t");
            }
            return;
        }



        if(constructFBottom) {
            data.cleanFBottom();
            if (node.isLeaf()) {
                processLeaf(nucleotideIndexMap, node, data, numReticulation);
            } else if (node.isTreeNode()) {
                processInternalTreeNode(node, data, isArticulation, siteID);
            } else if (node.isNetworkNode()) {
                processNetworkNode(node, data, numReticulation, reticulationID, siteID);
            }
            if (PRINT_DETAILS) {
                System.out.println("FBottoms:");
                data.printFMatrixs("b");
            }

            data.cleanFTop();
            processTopBranchOfNode(Q, node, data);
        }
        if(PRINT_DETAILS){
            System.out.println("FTop:");
            data.printFMatrixs("t");
        }
    }


    private static void examineBranchOfNode(QParameters Q, NetNode<SNAPPData[]> node, SNAPPData data)
    {
        Iterator<NetNode<SNAPPData[]>> parentIt = node.getParents().iterator();
        NetNode<SNAPPData[]> parent1 = parentIt.next();
        NetNode<SNAPPData[]> parent2 = parentIt.next();

        int mx = data.getFBottoms(parent1).get(data.getFBottoms(parent1).size() - 1).Item1.mx;
        int index = 0;
        double sum = 0;

        Map<Integer, Double> sums = new TreeMap<>();

        for (index = 0 ; index < data.getFBottoms(parent1).size() ; index++) {
            Tuple<FMatrix, int[]> fBot1 = data.getFBottoms(parent1).get(index);
            Tuple<FMatrix, int[]> fBot2 = data.getFBottoms(parent2).get(index);
            int n = fBot1.Item1.mx + fBot2.Item1.mx;
            if(!sums.containsKey(n))
                sums.put(n, 0.0);
            sums.put(n, sums.get(n) + fBot1.Item1.getSum() * fBot2.Item1.getSum());
        }

        for(Integer n : sums.keySet()) {
            System.out.println(n + "\t" + sums.get(n));
        }
/*
        for (int n = 1; n <= mx; n++) {
            //System.out.println("n=" + n);
            sum = 0;
            for (R r : R.loopOver(n)) {

                //System.out.println(r);
                for (R[] splitRPair : splittingR(r)) {
                    //System.out.println(r + " -> " + splitRPair[0] + " + " + splitRPair[1]);

                    Tuple<FMatrix, int[]> fBot1 = data.getFBottoms(parent1).get(index);
                    Tuple<FMatrix, int[]> fBot2 = data.getFBottoms(parent2).get(index);
                    //System.out.println(fBot2.Item2[0] + "\t" + fBot1.Item1.getSum() * fBot2.Item1.getSum() + "\t" + r + " -> " + splitRPair[0] + " + " + splitRPair[1]);
                    sum += fBot1.Item1.getSum() * fBot2.Item1.getSum();
                    index++;
                }
            }
            System.out.println(n + "\t" + sum);
        }
        System.out.println(sum);*/
    }


    /**
     * Process the top branch of a node.
     *
     * @param Q    The transition matrix for SNAPP.
     * @param node The node to process.
     */
    private static void processTopBranchOfNode(QParameters Q, NetNode<SNAPPData[]> node, SNAPPData data)
    {

        for(NetNode<SNAPPData[]> parent: node.getParents()) {
            if (Double.isNaN(node.getParentDistance(parent)) || Double.isInfinite(node.getParentDistance(parent)))
                throw new RuntimeException("Snapp only works with finite branch distances: " + node.getParentDistance(parent));

            double theta = node.getParentSupport(parent);
            MatrixQ matQ;
            if(Q._gTheta != null) {
                matQ = Q._gMatrix;
            }
            else {
                matQ = new MatrixQ(Q._rModel, Q._M, theta, !SWITCH_EXP_APPROX);
            }
            double t = node.getParentDistance(parent);

            if(!SWITCH_EXP_APPROX) {
                int maxColumnSize = 0;
                for(Tuple<FMatrix,int[]> fBot: data.getFBottoms(parent)) {
                    maxColumnSize = Math.max(maxColumnSize, fBot.Item1.getArr().length);
                }
                matQ.setTime(t, maxColumnSize);
            }
            if(data.getFBottoms(parent).isEmpty()) {

            } else {
                int index = 1;
                for (Tuple<FMatrix, int[]> fBot : data.getFBottoms(parent)) {
                    FMatrix fTop = new FMatrix(fBot.Item1.mx, fBot.Item1.hasEmptyR);
                    if (fTop.mx != 0 && !fBot.Item1.isArrAllZero()) {
                        //System.out.println(Arrays.toString(fBot.Item1.getArr()) + ": " +  fBot.Item1.isArrAllZero());
                        boolean foundInCache = false;
                        Long cacheAccessor = null;
                        if(SWITCH_FMATRIX_CACHE) {

                            cacheAccessor = Double.doubleToLongBits(theta);
                            cacheAccessor = (cacheAccessor << 3) ^ Double.doubleToLongBits(t);
                            cacheAccessor = (cacheAccessor << 27) ^ fBot.Item1.hash();

                            /*StringBuilder sb = new StringBuilder();
                            sb.append(String.format("%.6f ", theta));
                            sb.append(String.format("%.6f ", t));
                            for(int i = 0 ; i < fBot.Item1.getArr().length ; i++) {
                                sb.append(String.format("%.6f ", fBot.Item1.getArr()[i]));
                            }
                            cacheAccessor = sb.toString();*/
                            double[] arr = FMatrix.cache.get(cacheAccessor);
                            FMatrix.cacheAccess.add(1.0);
                            if(arr != null) {
                                fTop.setMatrix(arr);
                                foundInCache = true;
                                FMatrix.cacheHit.add(1.0);
                            }

                        }
                        if(!foundInCache){

                            if (SWITCH_EXP_APPROX)
                                fTop.setMatrix(matQ.expQTtx(t, fBot.Item1.getArr(), fBot.Item1.mx));
                            else
                                fTop.setMatrix(matQ.getProbabilityForColumn(fBot.Item1.getArr()));

                            if(SWITCH_FMATRIX_CACHE) {
                                FMatrix.cache.put(cacheAccessor, fTop.getArr());
                            }
                        }
                    }
                    data.addFTop(parent, fTop, fBot.Item2);
                    //if(node.getName().equals("I1"))
                    //    System.out.println(index + "\t" + fTop.getSum());
                    index++;
                }
            }
        }
    }



    /**
     * Process an internal non-leaf node.
     *
     * @param node The node to process.
     * @param data The data for that node.
     */
    private static void processInternalTreeNode(NetNode<SNAPPData[]> node, SNAPPData data, boolean isArticulation, int siteID) {

        if (node.getChildCount() != 2)
            throw new RuntimeException("SNAPP does not work on networks with internal tree node with 3 or more children per node");

        Iterator<NetNode<SNAPPData[]>> children = node.getChildren().iterator();
        SNAPPData childData1 = children.next().getData()[siteID];
        SNAPPData childData2 = children.next().getData()[siteID];

        NetNode parent = null;
        if(!node.isRoot()) {
            parent = node.getParents().iterator().next();
        }

        //System.out.println(childData1.getFTops(node).size() + " " + childData2.getFTops(node).size());
        boolean added = false;
        if(childData1.getFTops(node).size() < 100 || childData2.getFTops(node).size() < 100) {
        //if(true){

            for (Tuple<FMatrix, int[]> fTop1 : childData1.getFTops(node)) {
                for (Tuple<FMatrix, int[]> fTop2 : childData2.getFTops(node)) {
                    int[] newIndex = mergeTwoSplittingIndices(fTop1.Item2, fTop2.Item2);
                    if (newIndex != null) {
                        data.addFBottom(parent, getFBottomNormal(fTop1.Item1, fTop2.Item1), newIndex);
                        added = true;
                    }
                }
            }
        }
//        else {
//            Map<SplittingIndexVector, Set<Integer>> fTop2Cache = new HashMap<>();
//            int fTop2Index = 0;
//            for(Tuple<FMatrix,int[]> fTop2: childData2.getFTops(node)) {
//                SplittingIndexVector siv = new SplittingIndexVector(fTop2.Item2);
//                if(!fTop2Cache.containsKey(siv))
//                    fTop2Cache.put(siv, new HashSet<>());
//                fTop2Cache.get(siv).add(fTop2Index);
//                fTop2Index++;
//            }
//
//            for(Tuple<FMatrix,int[]> fTop1: childData1.getFTops(node)) {
//                SplittingIndexVector siv = new SplittingIndexVector(fTop1.Item2);
//                for(SplittingIndexVector e : siv.getAllCompatible()) {
//                    if(fTop2Cache.containsKey(e)) {
//                        for(int index : fTop2Cache.get(e)) {
//                            Tuple<FMatrix, int[]> fTop2 = childData2.getFTops(node).get(index);
//                            int[] newIndex = mergeTwoSplittingIndices(fTop1.Item2, fTop2.Item2);
//                            if (newIndex != null) {
//                                data.addFBottom(parent, getFBottomNormal(fTop1.Item1, fTop2.Item1), newIndex);
//                            }
//                        }
//                    }
//                }
//            }
//        }
        else {
            List<Map<Integer, Set<Integer>>> fTop2Cache = new ArrayList<>();
            int reticulationNumber = -1;
            int fTop2Index = 0;
            for(Tuple<FMatrix,int[]> fTop2: childData2.getFTops(node)) {
                if(reticulationNumber == -1) {
                    reticulationNumber = fTop2.Item2.length;
                    for(int i = 0 ; i < reticulationNumber ; i++) {
                        fTop2Cache.add(new HashMap<>());
                    }
                }
                for(int i = 0 ; i < reticulationNumber ; i++) {
                    if(fTop2Cache.get(i).get(fTop2.Item2[i]) == null)
                        fTop2Cache.get(i).put(fTop2.Item2[i], new HashSet<>());
                    fTop2Cache.get(i).get(fTop2.Item2[i]).add(fTop2Index);
                }
                fTop2Index++;
            }

            int maxsize = 0;

            if(reticulationNumber == 0) {
                for (Tuple<FMatrix, int[]> fTop1 : childData1.getFTops(node)) {
                    for (Tuple<FMatrix, int[]> fTop2 : childData2.getFTops(node)) {
                        int[] newIndex = mergeTwoSplittingIndices(fTop1.Item2, fTop2.Item2);
                        if (newIndex != null) {
                            data.addFBottom(parent, getFBottomNormal(fTop1.Item1, fTop2.Item1), newIndex);
                        }
                    }
                }
            } else {
                for(Tuple<FMatrix,int[]> fTop1: childData1.getFTops(node)) {
                    Set<Integer> stfTop2 = new HashSet<>();
                    if(fTop1.Item2[0] == 0) {
                        for(int k : fTop2Cache.get(0).keySet())
                            stfTop2.addAll(fTop2Cache.get(0).get(k));
                    } else {
                        if (fTop2Cache.get(0).containsKey(0))
                            stfTop2.addAll(fTop2Cache.get(0).get(0));
                        if (fTop2Cache.get(0).containsKey(fTop1.Item2[0]))
                            stfTop2.addAll(fTop2Cache.get(0).get(fTop1.Item2[0]));
                    }

                    for(int i = 1 ; i < reticulationNumber ; i++) {
                        Set<Integer> currentSet = new HashSet<>();
                        if(fTop1.Item2[i] == 0) {
                            for(int k : fTop2Cache.get(i).keySet())
                                currentSet.addAll(fTop2Cache.get(i).get(k));
                        } else {
                            if (fTop2Cache.get(i).containsKey(0))
                                currentSet.addAll(fTop2Cache.get(i).get(0));
                            if (fTop2Cache.get(i).containsKey(fTop1.Item2[i]))
                                currentSet.addAll(fTop2Cache.get(i).get(fTop1.Item2[i]));
                        }
                        stfTop2.retainAll(currentSet);
                    }

                    if(maxsize < stfTop2.size()) {
                        maxsize = stfTop2.size();
                        //System.out.println("stfTop2 size " + maxsize);
                    }
                    for(int index : stfTop2) {
                        Tuple<FMatrix, int[]> fTop2 = childData2.getFTops(node).get(index);
                        int[] newIndex = mergeTwoSplittingIndices(fTop1.Item2, fTop2.Item2);
                        if (newIndex != null) {
                            data.addFBottom(parent, getFBottomNormal(fTop1.Item1, fTop2.Item1), newIndex);
                            added = true;
                        }
                    }
                }
            }
        }
        if(!added) {
            int splittingIndexDimension = -1;
            for (Tuple<FMatrix, int[]> fTop1 : childData1.getFTops(node)) {
                splittingIndexDimension = fTop1.Item2.length;
                break;
            }
            for (Tuple<FMatrix, int[]> fTop2 : childData2.getFTops(node)) {
                splittingIndexDimension = fTop2.Item2.length;
                break;
            }
            int newIndex[] = new int[splittingIndexDimension];
            data.addFBottom(parent, new FMatrix(0, true), newIndex);
        }
        if(isArticulation){
            if(PRINT_DETAILS){
                System.out.println("FBottoms:");
                data.printFMatrixs("b");
            }
            data.cleanFBottomSplittingIndices();
        }
    }


    /**
     * The slow normal calcluation for F bottom.
     *
     * @param fTop1  FTop from one child.
     * @param fTop2  FTop from another child.
     * @return The F Bottom matrix.
     */
    private static FMatrix getFBottomNormal(FMatrix fTop1, FMatrix fTop2)
    {
        if(fTop1.ifHasEmptyR() && fTop2.ifHasEmptyR()) {
            return new FMatrix(0, true);
        }

        if(fTop1.hasEmptyR) {
            double u2[] = fTop2.getArr().clone();
            return new FMatrix(fTop2.mx, u2, false);
        }

        if(fTop2.hasEmptyR) {
            double u1[] = fTop1.getArr().clone();
            return new FMatrix(fTop1.mx, u1, false);
        }



        FMatrix FBottom = new FMatrix(fTop1.mx + fTop2.mx, fTop1.ifHasEmptyR() && fTop2.ifHasEmptyR());

        Long cacheAccessor = null;
        double[] cacheResult = null;
        if(SWITCH_FMATRIX_CACHE) {
            FMatrix.cacheAccess.add(1.0);
            cacheAccessor = (long) FBottom.mx;
            cacheAccessor = cacheAccessor * 11 + fTop1.hash() + fTop2.hash();
            cacheResult = FMatrix.cache.get(cacheAccessor);
            if(cacheResult != null) {
                FMatrix.cacheHit.add(1.0);
            }
        }

        if(cacheResult != null) {
        //if(false){
            FBottom.setMatrix(cacheResult);
        } else {

            if (R.dims == 1 && SWITCH_FASTER_BIALLILE) {
                //faster implementation
                double u1[] = fTop1.getArr().clone();
                double u2[] = fTop2.getArr().clone();
                for (int n = 1; n <= fTop1.mx; n++) {
                    double b = 1.0;
                    for (int i = 0; i <= n; i++) {
                        u1[n * (n + 1) / 2 - 1 + i] *= b;
                        b *= 1.0 * (n - i) / (i + 1);
                    }
                }
                for (int n = 1; n <= fTop2.mx; n++) {
                    double b = 1.0;
                    for (int i = 0; i <= n; i++) {
                        u2[n * (n + 1) / 2 - 1 + i] *= b;
                        b *= 1.0 * (n - i) / (i + 1);
                    }
                }

                double fb[] = FBottom.getArr();
                for (int n1 = 1; n1 <= fTop1.mx; n1++) {
                    for (int i = 0; i <= n1; i++) {
                        double f11 = u1[n1 * (n1 + 1) / 2 - 1 + i];
                        for (int n2 = 1; n2 <= fTop2.mx; n2++) {
                            for (int j = 0; j <= n2; j++) {
                                fb[(n1 + n2) * (n1 + n2 + 1) / 2 - 1 + (i + j)] += f11 * u2[n2 * (n2 + 1) / 2 - 1 + j];
                            }
                        }
                    }
                }

                for (int n = 1; n <= FBottom.mx; n++) {
                    double b = 1.0;
                    for (int i = 0; i <= n; i++) {
                        fb[n * (n + 1) / 2 - 1 + i] = Math.max(0.0, fb[n * (n + 1) / 2 - 1 + i] / b);
                        b *= 1.0 * (n - i) / (i + 1);
                    }
                }


            } else {

                FBottom = new FMatrix(fTop1.mx + fTop2.mx, fTop1.ifHasEmptyR() && fTop2.ifHasEmptyR());
                for (int n = 1; n <= FBottom.mx; n++) {
                    for (R r : R.loopOver(n))
                        FBottom.set(r, getFBottom(n, r, fTop1, fTop2));
                }
            }

            /*if(cacheResult != null) {
                for(int i = 0 ; i < FBottom.getArr().length ; i++) {
                    if(Math.abs(FBottom.getArr()[i] - cacheResult[i]) > 1e-6)
                        throw new RuntimeException("");
                }
            }*/

            if(SWITCH_FMATRIX_CACHE && cacheAccessor != null) {
                FMatrix.cache.put(cacheAccessor, FBottom.getArr().clone());
            }
        }

        return FBottom;
    }

    /**
     * Calculates an element of the F bottom matrix for a non-leaf node
     *
     * @param nBot The number of lineages used to index the F matrix.
     * @param rBot The R vector used to index the F matrix.
     * @param fTop1  FTop from one child.
     * @param fTop2  FTop from another child.
     * @return The F matrix value for the given n and r for that node.
     */
    private static double getFBottom(int nBot, R rBot, FMatrix fTop1, FMatrix fTop2)
    {

        double sum = 0;
        for (int nTop = 0; nTop <= nBot; nTop++)
        {

            if (nTop > fTop1.mx || (nBot - nTop) > fTop2.mx)
                continue;

            for (R rTop : R.loopOver(nTop))
            {

                R otherRTop = rBot.subtract(rTop);

                if (otherRTop == null)
                    continue;

                sum += fTop1.get(rTop) * fTop2.get(otherRTop) *
                        rBot.getProbabilityOfSelecting(rTop);
                //System.out.println(leftChild.getData().FTop.get(rTop) + "*" + rightChild.getData().FTop.get(otherRTop) + "*" +rBot.getProbabilityOfSelecting(rTop));
            }
        }
        return sum;
    }

    private static void processNetworkNodeSampling(NetNode<SNAPPData[]> node, SNAPPData data, int numReticulations, int reticulationID, Splitting splitting, int siteID) {
        if (node.getChildCount() != 1)
            throw new RuntimeException("SNAPP does not work on networks with reticulation node with 2 or more children per node");

        int splittingIndexDimension = -1;
        double[] inheritanceProbs = new double[2];
        NetNode[] parents = new NetNode[2];
        int index = 0;
        for(NetNode parent: node.getParents()){
            inheritanceProbs[index] = node.getParentProbability(parent);
            parents[index++] = parent;
        }

        if(node.getChildren().iterator().next().getData()[siteID].getFTops(node).size() > 1) {
            throw new RuntimeException("Not in sampling mode!");
        }

        double weight_sum = 0.0;
        R sampled_r = null;

        Tuple<FMatrix,int[]> tuple = node.getChildren().iterator().next().getData()[siteID].getFTops(node).get(0);
        FMatrix fTop = tuple.Item1;
        if(fTop.mx == 0) {
            int[] splittingIndex = tuple.Item2.clone();
            //splittingIndex[reticulationID] = index++;
            FMatrix fm = new FMatrix(0, true);
            data.addFBottom(parents[0], fm, splittingIndex);
            data.addFBottom(parents[1], fm, splittingIndex);
        } else {
            int mx = fTop.mx;
            int k = 1;
            for (int n = 1; n <= mx; n++) {
                for (R r : R.loopOver(n)) {
                    double prob = fTop.get(r);

                    if(sampled_r == null) {
                        sampled_r = new R(r);;
                        weight_sum += prob;
                    } else {
                        weight_sum += prob;
                        double p = prob / weight_sum;
                        double j = Randomizer.getRandomDouble();
                        if (j <= p) {
                            sampled_r = new R(r);
                        }
                    }
                }
            }

            int n0 = 0, n1 = 0; // to parents[0]
            for(int i0 = 1 ; i0 <= sampled_r.getNum(0) ; i0++) {
                double p = Randomizer.getRandomDouble();
                if(p < inheritanceProbs[0]) {
                    n0++;
                    splitting.logWeight += Math.log(inheritanceProbs[0]);
                } else {
                    splitting.logWeight += Math.log(inheritanceProbs[1]);
                }
            }

            for(int i1 = 1 ; i1 <= sampled_r.getNum(1) ; i1++) {
                double p = Randomizer.getRandomDouble();
                if(p < inheritanceProbs[0]) {
                    n1++;
                    splitting.logWeight += Math.log(inheritanceProbs[0]);
                } else {
                    splitting.logWeight += Math.log(inheritanceProbs[1]);
                }
            }

            R[] splitRPair = new R[]{new R(n0 + n1, new int[]{n0}), new R(sampled_r.n - n0 - n1, new int[]{sampled_r.getNum(0) - n0})};
            double prob = fTop.get(sampled_r);

            splitting.logWeight += Math.log(prob / weight_sum);

            boolean probSet = false;
            for(int i=0; i<2; i++){
                FMatrix fm = new FMatrix(splitRPair[i].n, splitRPair[i].n==0);
                if(fm.mx != 0){
                    if(probSet){
                        fm.set(splitRPair[i],1);
                    }
                    else{
                        double weight = sampled_r.getProbabilityOfSelecting(splitRPair[0]) / calculateRWeight(sampled_r, splitRPair[0]);
                        //System.out.println(splitRPair[0].n + " vs. " + splitRPair[1].n + ": " + Math.pow(inheritanceProbs[0],splitRPair[0].n) + " * " + Math.pow(inheritanceProbs[1],splitRPair[1].n) + " * " + weight);
                        double newProb = prob / weight * Math.pow(inheritanceProbs[0],splitRPair[0].n) * Math.pow(inheritanceProbs[1],splitRPair[1].n);
                        //System.out.println(splitRPair[0].n + " vs. " + splitRPair[1].n + ": " + newProb);
                        //System.out.println(index);
                        fm.set(splitRPair[i],newProb);
                        probSet = true;
                    }
                }
                int[] splittingIndex = tuple.Item2.clone();
                splittingIndex[reticulationID] = 0;
                data.addFBottom(parents[i], fm, splittingIndex);
            }
        }

    }

    /**
     * Process a network node.
     *
     * @param node The node to process.
     * @param data The data for that node.
     */
    private static void processNetworkNode(NetNode<SNAPPData[]> node, SNAPPData data, int numReticulations, int reticulationID, int siteID)
    {

        if (node.getChildCount() != 1)
            throw new RuntimeException("SNAPP does not work on networks with reticulation node with 2 or more children per node");

        int added[] = new int[]{0, 0};
        int splittingIndexDimension = -1;
        double[] inheritanceProbs = new double[2];
        NetNode[] parents = new NetNode[2];
        int index = 0;
        for(NetNode parent: node.getParents()){
            inheritanceProbs[index] = node.getParentProbability(parent);
            parents[index++] = parent;
        }

        index = 1;
        for(Tuple<FMatrix,int[]> tuple: node.getChildren().iterator().next().getData()[siteID].getFTops(node)){
            index = 1;
            FMatrix fTop = tuple.Item1;
            if(splittingIndexDimension == -1)
                splittingIndexDimension = tuple.Item2.length;
            if(fTop.mx == 0) {
                int[] splittingIndex = tuple.Item2.clone();
                //splittingIndex[reticulationID] = index++;
                FMatrix fm = new FMatrix(0, true);
                data.addFBottom(parents[0], fm, splittingIndex);
                data.addFBottom(parents[1], fm, splittingIndex);
            }
            int mx = fTop.mx;
            for (int n = 1; n <= mx; n++)
            {
                //System.out.println("n=" + n);
                for (R r : R.loopOver(n)) {
                    double prob = fTop.get(r);
                    /*if(SWITCH_APPROX_SPLIT) {
                        if (prob < 1e-3)
                            continue;
                    }*/

                    //System.out.println(r);
                    for(R[] splitRPair: splittingR(r)){
                        //System.out.println(r.getNum(0) + "\t" + r.getNum(1) + "\t" + splitRPair[0].getNum(0) + "\t" + splitRPair[0].getNum(1) + "\t" + splitRPair[1].getNum(0) + "\t" + splitRPair[1].getNum(1) + "\t" + r + " -> " + splitRPair[0] + " + " + splitRPair[1]);

                        boolean probSet = false;
                        for(int i=0; i<2; i++){
                            FMatrix fm = new FMatrix(splitRPair[i].n, splitRPair[i].n==0);
                            if(fm.mx != 0){
                                if(probSet){
                                    fm.set(splitRPair[i],1);
                                }
                                else{
                                    //int weight = calculateRWeight(r, splitRPair[0]);
                                    double weight = r.getProbabilityOfSelecting(splitRPair[0]) / calculateRWeight(r, splitRPair[0]);
                                    //System.out.println(splitRPair[0].n + " vs. " + splitRPair[1].n + ": " + Math.pow(inheritanceProbs[0],splitRPair[0].n) + " * " + Math.pow(inheritanceProbs[1],splitRPair[1].n) + " * " + weight);
                                    double newProb = prob / weight * Math.pow(inheritanceProbs[0],splitRPair[0].n) * Math.pow(inheritanceProbs[1],splitRPair[1].n);
                                    //System.out.println(splitRPair[0].n + " vs. " + splitRPair[1].n + ": " + newProb);
                                    //System.out.println(index);
                                    fm.set(splitRPair[i],newProb);
                                    probSet = true;
                                }
                            }
                            int[] splittingIndex = tuple.Item2.clone();
                            splittingIndex[reticulationID] = index;
                            /*if(SWITCH_APPROX_SPLIT) {
                                double threshold = 1e-10;
                                if(numReticulations == 1)
                                    threshold = 1e-5;
                                else if(numReticulations == 2)
                                    threshold = 1e-4;
                                else if(numReticulations > 2)
                                    threshold = 1e-3;
                                threshold = 1e-5;
                                if (!fm.ifHasEmptyR() && fm.getSum() < threshold && added[i] > 0) {
                                    i = 2;
                                    continue;
                                }
                            }*/
                            added[i]++;
                            data.addFBottom(parents[i], fm, splittingIndex);
                        }
                        index++;
                    }
                    //System.out.println(r.n + ": " + totalweight);
                }
            }
        }

        data.pruneFBottom(targetSplittingIndices);
    }

    private static void processNetworkNode2(NetNode<SNAPPData[]> node, SNAPPData data, int numReticulations, int reticulationID, int siteID)
    {
        if (node.getChildCount() != 1)
            throw new RuntimeException("SNAPP does not work on networks with reticulation node with 2 or more children per node");

        int added[] = new int[]{0, 0};
        int splittingIndexDimension = -1;
        double[] inheritanceProbs = new double[2];
        NetNode[] parents = new NetNode[2];
        int index = 0;
        for(NetNode parent: node.getParents()){
            inheritanceProbs[index] = node.getParentProbability(parent);
            parents[index++] = parent;
        }

        index = 1;
        for(Tuple<FMatrix,int[]> tuple: node.getChildren().iterator().next().getData()[siteID].getFTops(node)){
            index = 1;
            FMatrix fTop = tuple.Item1;
            if(splittingIndexDimension == -1)
                splittingIndexDimension = tuple.Item2.length;
            if(fTop.mx == 0) {
                int[] splittingIndex = tuple.Item2.clone();
                //splittingIndex[reticulationID] = index++;
                FMatrix fm = new FMatrix(0, true);
                data.addFBottom(parents[0], fm, splittingIndex);
                data.addFBottom(parents[1], fm, splittingIndex);
            }
            int mx = fTop.mx;
            if(SWITCH_APPROX_SPLIT) {
                mx = Math.min(mx, 5);
            }
            for (int n = 1; n <= mx; n++)
            {
                //System.out.println("n=" + n);
                for (R r : R.loopOver(n)) {
                    double prob = fTop.get(r);
                    /*if(SWITCH_APPROX_SPLIT) {
                        if (prob < 1e-3)
                            continue;
                    }*/

                    //System.out.println(r);
                    for(R[] splitRPair: splittingR(r)){

                        boolean probSet = false;
                        for(int i=0; i<2; i++){
                            FMatrix fm = new FMatrix(splitRPair[i].n, splitRPair[i].n==0);
                            if(fm.mx != 0){
                                if(probSet){
                                    fm.set(splitRPair[i],1);
                                }
                                else{
                                    //int weight = calculateRWeight(r, splitRPair[0]);
                                    double weight = r.getProbabilityOfSelecting(splitRPair[0]) / calculateRWeight(r, splitRPair[0]);
                                    //System.out.println(splitRPair[0].n + " vs. " + splitRPair[1].n + ": " + Math.pow(inheritanceProbs[0],splitRPair[0].n) + " * " + Math.pow(inheritanceProbs[1],splitRPair[1].n) + " * " + weight);
                                    double newProb = prob / weight * Math.pow(inheritanceProbs[0],splitRPair[0].n) * Math.pow(inheritanceProbs[1],splitRPair[1].n);
                                    //System.out.println(splitRPair[0].n + " vs. " + splitRPair[1].n + ": " + newProb);
                                    fm.set(splitRPair[i],newProb);
                                    probSet = true;
                                }
                            }
                            int[] splittingIndex = tuple.Item2.clone();
                            splittingIndex[reticulationID] = index;
                            /*if(SWITCH_APPROX_SPLIT) {
                                double threshold = 1e-10;
                                if(numReticulations == 1)
                                    threshold = 1e-5;
                                else if(numReticulations == 2)
                                    threshold = 1e-4;
                                else if(numReticulations > 2)
                                    threshold = 1e-3;
                                threshold = 1e-5;
                                if (!fm.ifHasEmptyR() && fm.getSum() < threshold && added[i] > 0) {
                                    i = 2;
                                    continue;
                                }
                            }*/
                            added[i]++;
                            data.addFBottom(parents[i], fm, splittingIndex);
                        }
                        index++;
                    }
                    //System.out.println(r.n + ": " + totalweight);
                }
            }
        }
        data.cleanFBottomSplittingIndices2();

    }


    protected static int calculateRWeight(R fromR, R toR){
        int weight = 1;
        for(int type = 0; type<(R.dims + 1); type++){
            weight *= ArithmeticUtils.binomialCoefficient(fromR.getNum(type), toR.getNum(type));
        }
        return weight;
    }






    public static List<R[]> splittingR(R r){
        List<R[]> splittingRs = new ArrayList<R[]>();

        if(R.dims == 3) {
            for (int A = 0; A <= r.getNum(0); A++)
                for (int C = 0; C <= r.getNum(1); C++)
                    for (int T = 0; T <= r.getNum(2); T++)
                        for (int G = 0; G <= r.getNum(3); G++) {
                            int totalN = A + C + T + G;
                            int[] rValues = new int[3];
                            rValues[0] = A;
                            rValues[1] = C;
                            rValues[2] = T;
                            R[] splitR = new R[2];
                            splitR[0] = new R(totalN, rValues);
                            splitR[1] = r.subtract(splitR[0]);
                            splittingRs.add(splitR);
                        }
        } else if(R.dims == 1) {
            for(int zero = 0 ; zero <= r.getNum(0) ; zero ++) {
                for(int one = 0 ; one <= r.getNum(1) ; one++) {
                    int totalN = zero + one;
                    int[] rValues = new int[1];
                    rValues[0] = zero;
                    R[] splitR = new R[2];
                    splitR[0] = new R(totalN, rValues);
                    splitR[1] = r.subtract(splitR[0]);
                    splittingRs.add(splitR);
                }
            }
        }
        return splittingRs;
    }



    /**
     * Perform the initial processing of a leaf.
     *
     * @param nucleotideIndexMap The observations encoded as indices.
     * @param node               The node to process.
     * @param data               The data for that node.
     */
    protected static void processLeaf(Map<String, R> nucleotideIndexMap, NetNode<SNAPPData[]> node, SNAPPData data, int numReticulations)
    {
        NetNode parent = node.getParents().iterator().next();
        int[] splittingIndex = new int[numReticulations];
        R r = nucleotideIndexMap.get(node.getName());
        FMatrix fBot;

        if(r == null || r.n == 0) {
            data.addFBottom(parent, 0, true, splittingIndex);
            return;
        }

        if(!HAS_DOMINANT_MARKER)
            fBot = data.addFBottom(parent, r.n, false, splittingIndex);
        else
            fBot = data.addFBottom(parent, r.n * 2, false, splittingIndex);
        //fBot = data.addFBottom(parent, r.n, false, splittingIndex);

        int weight = 1;

        if(CORRECTION_AT_LEAVES) {
            int total = r.n;
            for(int i = 0; i<(R.dims + 1); i++){
                weight *= ArithmeticUtils.binomialCoefficient(total, r.getNum(i));
                total = total - r.getNum(i);
            }
        }

        if(HAS_DOMINANT_MARKER) {
            int N = r.n * 2;
            double p_r_k_n = 1.0; //p(0,0,n)=1
            int nReds = r.values[0];
            if(nReds > 0) {
                int n = r.n;
                for (int i = 1; i <= nReds; i++)
                    p_r_k_n = (p_r_k_n * 2.0 * (n - i + 1.0)) / (2.0 * n - i + 1.0);
                for (int k = nReds; k <= 2 * nReds; k++) {
                    if (k > nReds)
                        p_r_k_n = (p_r_k_n * (2.0 * nReds - k + 1) * k) / (2.0 * (k - nReds) * (2.0 * n - k + 1.0));
                    fBot.getArr()[N * (N + 1) / 2 - 1 + k] = p_r_k_n / weight;
                }
            } else {
                fBot.getArr()[N * (N + 1) / 2 - 1 + r.values[0]] = 1.0 / weight;
            }
        } else {

            if (R.dims == 1 && SWITCH_FASTER_BIALLILE)
                fBot.getArr()[r.n * (r.n + 1) / 2 - 1 + r.values[0]] = 1.0 / weight;
            else
                fBot.set(r, 1.0 / weight);
        }
    }


    private static Network cloneNetwork(Network net) {
        return Networks.readNetwork(net.toString());
    }




}




