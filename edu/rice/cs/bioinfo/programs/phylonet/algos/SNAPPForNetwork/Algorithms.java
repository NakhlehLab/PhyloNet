package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;


import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.FelsensteinAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import jeigen.DenseMatrix;
import org.apache.commons.math3.util.ArithmeticUtils;

import java.util.*;

/**
 * This class holds all the various algorithms from http://mbe.oxfordjournals.org/content/29/8/1917.
 */
public class Algorithms
{
    private static final boolean PRINT_DETAILS = false;


    private static int[] mergeTwoSplittingIndices(int[] index1, int[] index2){
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
            int[] Rvalue = new int[3];
            int n = 0;
            for(int i=0; i<3; i++){
                Rvalue[i] = entry.getValue()[i];
                n += entry.getValue()[i];
            }
            n += entry.getValue()[3];
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
     * @param obs         The nucleotide observations.
     * @param Q           The transition matrix for SNAPP.
     * @return The probability.
     */
    @SuppressWarnings("unchecked")
   public static double getProbabilityObservationGivenNetwork(Network<SNAPPData[]> speciesNetwork, Map<String,String> allele2species, NucleotideObservation obs, MatrixQ Q, int siteID)
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


    public static double updateProbabilityObservationGivenNetwork(Network<SNAPPData[]> speciesNetwork, MatrixQ Q, int siteID, Map<NetNode, Boolean> node2update)
    {

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
            updateNode(Q, node, isArticulation, reticulationID, siteID, node2update.get(node));


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
    private static double getProbabilityOfNetwork(Network<SNAPPData[]> speciesNetwork, MatrixQ Q, int siteID)
    {
        FMatrix rootFBot = speciesNetwork.getRoot().getData()[siteID].getFBottoms(null).iterator().next().Item1;

        DenseMatrix eq = Q.getEquilibrium();
        double sum = 0;
        for (int n = 1; n <= rootFBot.mx; n++)
        {
            for (R r : R.loopOver(n))
            {
                sum += rootFBot.get(r) * eq.get(r.getIndex(), 0);
            }
        }
        return sum;
    }

    /**
     * Process a node for the snapp algorithm.
     *
     * @param Q                  The transition matrix for SNAPP.
     * @param nucleotideIndexMap The observations encoded as indices.
     * @param node               The node to process.
     */
    private static void processNode(MatrixQ Q, Map<String, R> nucleotideIndexMap, NetNode<SNAPPData[]> node, boolean isArticulation, int numReticulations, int reticulationID, int siteID)
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
            processNetworkNode(node, data, reticulationID, siteID);
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
    private static void updateNode(MatrixQ Q, NetNode<SNAPPData[]> node, boolean isArticulation, int reticulationID, int siteID, Boolean constructFBottom)
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
            if (node.isTreeNode()) {
                processInternalTreeNode(node, data, isArticulation, siteID);
            } else if (node.isNetworkNode()) {
                processNetworkNode(node, data, reticulationID, siteID);
            }
            if (PRINT_DETAILS) {
                System.out.println("FBottoms:");
                data.printFMatrixs("b");
            }
        }
        data.cleanFTop();
        processTopBranchOfNode(Q, node, data);
        if(PRINT_DETAILS){
            System.out.println("FTop:");
            data.printFMatrixs("t");
        }
    }





    /**
     * Process the top branch of a node.
     *
     * @param Q    The transition matrix for SNAPP.
     * @param node The node to process.
     */
    private static void processTopBranchOfNode(MatrixQ Q, NetNode<SNAPPData[]> node, SNAPPData data)
    {

        for(NetNode<SNAPPData[]> parent: node.getParents()) {
            if (Double.isNaN(node.getParentDistance(parent)) || Double.isInfinite(node.getParentDistance(parent)))
                throw new RuntimeException("Snapp only works with finite branch distances: " + node.getParentDistance(parent));

            for(Tuple<FMatrix,int[]> fBot: data.getFBottoms(parent)) {
                FMatrix fTop = new FMatrix(fBot.Item1.mx, fBot.Item1.hasEmptyR);
                if(fTop.mx!=0 && !fBot.Item1.isArrAllZero()) {
                    //System.out.println(Arrays.toString(fBot.Item1.getArr()) + ": " +  fBot.Item1.isArrAllZero());
                    fTop.setMatrix(Q.getProbabilityForColumn(node.getParentDistance(parent), fBot.Item1.getArr()));
                }
                data.addFTop(parent, fTop, fBot.Item2);
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


        for(Tuple<FMatrix,int[]> fTop1: childData1.getFTops(node)){
            for(Tuple<FMatrix,int[]> fTop2: childData2.getFTops(node)){
                int[] newIndex = mergeTwoSplittingIndices(fTop1.Item2, fTop2.Item2);
                if(newIndex != null){
                    data.addFBottom(parent, getFBottomNormal(fTop1.Item1, fTop2.Item1), newIndex);
                }
            }
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
        FMatrix FBottom = new FMatrix(fTop1.mx + fTop2.mx, fTop1.ifHasEmptyR() && fTop2.ifHasEmptyR());
        for (int n = 1; n <= FBottom.mx; n++)
        {
            for (R r : R.loopOver(n))
                FBottom.set(r, getFBottom(n, r, fTop1, fTop2));
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



    /**
     * Process a network node.
     *
     * @param node The node to process.
     * @param data The data for that node.
     */
    private static void processNetworkNode(NetNode<SNAPPData[]> node, SNAPPData data, int reticulationID, int siteID)
    {
        if (node.getChildCount() != 1)
            throw new RuntimeException("SNAPP does not work on networks with reticulation node with 2 or more children per node");


        double[] inheritanceProbs = new double[2];
        NetNode[] parents = new NetNode[2];
        int index = 0;
        for(NetNode parent: node.getParents()){
            inheritanceProbs[index] = node.getParentProbability(parent);
            parents[index++] = parent;
        }

        index = 1;
        for(Tuple<FMatrix,int[]> tuple: node.getChildren().iterator().next().getData()[siteID].getFTops(node)){
            FMatrix fTop = tuple.Item1;
            for (int n = 1; n <= fTop.mx; n++)
            {
                for (R r : R.loopOver(n)) {
                    double prob = fTop.get(r);
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
                                    fm.set(splitRPair[i],newProb);
                                    probSet = true;

                                }
                            }
                            int[] splittingIndex = tuple.Item2.clone();
                            splittingIndex[reticulationID] = index;
                            data.addFBottom(parents[i], fm, splittingIndex);
                        }
                        index++;
                    }
                    //System.out.println(r.n + ": " + totalweight);
                }
            }
        }
    }


    private static int calculateRWeight(R fromR, R toR){
        int weight = 1;
        for(int type=0; type<4; type++){
            weight *= ArithmeticUtils.binomialCoefficient(fromR.getNum(type), toR.getNum(type));
        }
        return weight;
    }






    private static List<R[]> splittingR(R r){
        List<R[]> splittingRs = new ArrayList<R[]>();
        for(int A=0; A<=r.getNum(0); A++)
            for(int C=0; C<=r.getNum(1); C++)
                for(int T=0; T<=r.getNum(2); T++)
                    for(int G=0; G<=r.getNum(3); G++){
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
        return splittingRs;
    }



    /**
     * Perform the initial processing of a leaf.
     *
     * @param nucleotideIndexMap The observations encoded as indices.
     * @param node               The node to process.
     * @param data               The data for that node.
     */
    private static void processLeaf(Map<String, R> nucleotideIndexMap, NetNode<SNAPPData[]> node, SNAPPData data, int numReticulations)
    {
        NetNode parent = node.getParents().iterator().next();
        int[] splittingIndex = new int[numReticulations];
        R r = nucleotideIndexMap.get(node.getName());
        FMatrix fBot = data.addFBottom(parent, r.n, false, splittingIndex);
        int weight = 1;
        int total = r.n;
        for(int i=0; i<4; i++){
            weight *= ArithmeticUtils.binomialCoefficient(total, r.getNum(i));
            total = total - r.getNum(i);
        }
        fBot.set(r , 1.0/weight);
    }


    private static Network cloneNetwork(Network net) {
        return Networks.readNetwork(net.toString());
    }




}




