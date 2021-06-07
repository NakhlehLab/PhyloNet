package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.Splitting;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.Algorithms.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 3/2/19
 * Time: 3:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class AlgorithmsWithSampling {
    private static boolean PRINT_DETAILS = false;

    public static double getProbabilityObservationGivenNetwork(Network<SNAPPData[]> speciesNetwork, RPattern pattern, QParameters Q, Splitting splitting, int siteID)
    {
        Map<String, R> nucleotideIndexMap = pattern.getPattern();
        Set<String> articulationNodes = new HashSet<String>();
        for (Object node : Networks.getLowestArticulationNodes(speciesNetwork.clone())) {
            articulationNodes.add(((NetNode) node).getName());
        }
        int numReticulations = speciesNetwork.getReticulationCount();
        int reticulationID = 0;

        if(splitting.CANCELED_LIKELIHOOD) {
            double prob = 1.0;
            for (NetNode<SNAPPData[]> node : Networks.postTraversal(speciesNetwork)) {
                boolean isArticulation = articulationNodes.contains(node.getName());
                if(node.isLeaf())
                    prob *= processNode(Q, nucleotideIndexMap, node, isArticulation, numReticulations, reticulationID, splitting, siteID);
                if (node.isNetworkNode()) {
                    reticulationID++;
                }
            }

            return prob;
        } else {
            double prob = 1.0;

            for (NetNode<SNAPPData[]> node : Networks.postTraversal(speciesNetwork)) {
                boolean isArticulation = articulationNodes.contains(node.getName());
                prob *= processNode(Q, nucleotideIndexMap, node, isArticulation, numReticulations, reticulationID, splitting, siteID);
                if (node.isNetworkNode()) {
                    reticulationID++;
                }
            }

            return getProbabilityOfNetwork(speciesNetwork, Q, siteID);
        }
    }

    private static double processNode(QParameters Q, Map<String, R> nucleotideIndexMap, NetNode<SNAPPData[]> node, boolean isArticulation, int numReticulations, int reticulationID, Splitting splitting, int siteID)
    {
        double prob = 1.0;
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
            processInternalTreeNode(node, data, isArticulation, splitting, siteID);
        } else if(node.isNetworkNode()){
            prob *= processNetworkNode(node, data, numReticulations, reticulationID, splitting, siteID);
        }
        if(PRINT_DETAILS){
            System.out.println("FBottoms:");
            data.printFMatrixs("b");
        }

        prob *= processTopBranchOfNode(Q, node, data, splitting);
        if(PRINT_DETAILS){
            System.out.println("FTop:");
            data.printFMatrixs("t");
        }

        return prob;
    }

    private static FMatrix getFBottomNormal(FMatrix fTop1, FMatrix fTop2, Splitting splitting)
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



        FMatrix FBottom;

        FBottom = new FMatrix(fTop1.mx + fTop2.mx, fTop1.ifHasEmptyR() && fTop2.ifHasEmptyR());
        R r1 = null;
        R r2 = null;
        double p1 = 0;
        double p2 = 0;
        for (int n = 1; n <= fTop1.mx; n++) {
            for (R r : R.loopOver(n)) {
                if(fTop1.get(r) > 0) {
                    r1 = new R(r);
                    p1 = fTop1.get(r);
                }
            }
        }

        for (int n = 1; n <= fTop2.mx; n++) {
            for (R r : R.loopOver(n)) {
                if(fTop2.get(r) > 0) {
                    r2 = new R(r);
                    p2 = fTop2.get(r);
                }
            }
        }

        int nn[] = new int[]{r1.getNum(0) + r2.getNum(0)};
        R rr = new R((r1.getN() + r2.getN()), nn);

        FBottom.set(rr, p1 * p2 * rr.getProbabilityOfSelecting(r1));
        //splitting.logWeight += Math.log(rr.getProbabilityOfSelecting(r1));

        return FBottom;
    }

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

    private static void processInternalTreeNode(NetNode<SNAPPData[]> node, SNAPPData data, boolean isArticulation, Splitting splitting, int siteID) {

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
        for (Tuple<FMatrix, int[]> fTop1 : childData1.getFTops(node)) {
            for (Tuple<FMatrix, int[]> fTop2 : childData2.getFTops(node)) {
                int[] newIndex = mergeTwoSplittingIndices(fTop1.Item2, fTop2.Item2);
                if (newIndex != null) {
                    data.addFBottom(parent, getFBottomNormal(fTop1.Item1, fTop2.Item1, splitting), newIndex);
                    added = true;
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

    private static double processTopBranchOfNode(QParameters Q, NetNode<SNAPPData[]> node, SNAPPData data, Splitting splitting)
    {
        double final_prob = 1.0;
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
                        if (SWITCH_EXP_APPROX)
                            fTop.setMatrix(matQ.expQTtx(t, fBot.Item1.getArr(), fBot.Item1.mx));
                        else
                            fTop.setMatrix(matQ.getProbabilityForColumn(fBot.Item1.getArr()));

                        int mx = fTop.mx;
                        R sampled_r = splitting.lineageConfigurations.getTop(node, parent);
                        double sampled_v = fTop.get(sampled_r);
                        for (int n = 1; n <= mx; n++) {
                            for (R r : R.loopOver(n)) {
                                double prob = fTop.get(r);
                                fTop.set(r, 0.0);
                            }
                        }
                        final_prob *= sampled_v;
                        fTop.set(sampled_r, sampled_v);

                    }
                    data.addFTop(parent, fTop, fBot.Item2);
                    //if(node.getName().equals("I1"))
                    //    System.out.println(index + "\t" + fTop.getSum());
                    index++;
                }
            }
        }

        return final_prob;
    }

    private static double processNetworkNode(NetNode<SNAPPData[]> node, SNAPPData data, int numReticulations, int reticulationID, Splitting splitting, int siteID)
    {

        double total = 1.0;
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

                    if(prob > 0) {

                        R[] splitRPair = new R[]{new R(splitting.lineageConfigurations.getBottom(node, parents[0])), new R(splitting.lineageConfigurations.getBottom(node, parents[1]))};

                        boolean probSet = false;
                        for(int i=0; i<2; i++){
                            FMatrix fm = new FMatrix(splitRPair[i].n, splitRPair[i].n==0);
                            if(fm.mx != 0){
                                if(probSet){
                                    fm.set(splitRPair[i],1);
                                }
                                else{
                                    double weight = r.getProbabilityOfSelecting(splitRPair[0]) / calculateRWeight(r, splitRPair[0]);
                                    //System.out.println(splitRPair[0].n + " vs. " + splitRPair[1].n + ": " + Math.pow(inheritanceProbs[0],splitRPair[0].n) + " * " + Math.pow(inheritanceProbs[1],splitRPair[1].n) + " * " + weight);
                                    double newProb = prob / weight * Math.pow(inheritanceProbs[0],splitRPair[0].n) * Math.pow(inheritanceProbs[1],splitRPair[1].n);
                                    //System.out.println(splitRPair[0].n + " vs. " + splitRPair[1].n + ": " + newProb);
                                    //System.out.println(index);
                                    total *= newProb;
                                    fm.set(splitRPair[i],newProb);
                                    probSet = true;
                                }
                            }
                            int[] splittingIndex = tuple.Item2.clone();
                            splittingIndex[reticulationID] = 0;
                            data.addFBottom(parents[i], fm, splittingIndex);
                        }

                        return total;
                    }

                }
            }
        }
        return total;

    }



}
