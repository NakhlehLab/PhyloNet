package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

import java.util.*;

/**
 * Created by yunyu on 7/19/14.
 */

/**
 * This class represents the data SNAPP associates with each branch.
 * Note that Tree nodes correspond to the branches above them.
 * IE: (A,B)C.
 * Node A would hold the data for branch A-C.
 * Node B would hold the data for branch B-C.
 * Node C would hold the data for branch C-infinity.
 */
public class SNAPPData {

    int maxMX;

    /**
     * The "nucleotideIndex" of a node.
     * This index is what is used in R.getNum(index).
     */
    int nucleotideIndex;

    /**
     * The F Bottom matrix as described in the article.
     */
    Map<NetNode, List<Tuple<FMatrix, int[]>>> FBottomsMap;

    /**
     * The F Top matrix as described in the article.
     */
    Map<NetNode, List<Tuple<FMatrix, int[]>>> FTopsMap;


    public SNAPPData(){
        FBottomsMap = new HashMap<NetNode, List<Tuple<FMatrix, int[]>>>();
        FTopsMap = new HashMap<NetNode, List<Tuple<FMatrix, int[]>>>();
        maxMX = 0;
    }

    public void printFMatrixs(String type){
        Map<NetNode, List<Tuple<FMatrix, int[]>>> matricesToPrint = null;
        if(type.toLowerCase().equals("t")){
            matricesToPrint = FTopsMap;
        }else if(type.toLowerCase().equals("b")){
            matricesToPrint = FBottomsMap;
        }

        for(Map.Entry<NetNode, List<Tuple<FMatrix, int[]>>> FBottoms: matricesToPrint.entrySet()){
            NetNode parent = FBottoms.getKey();
            if(parent!=null) {
                System.out.println("To parent node " + parent.getName());
            }
            for(Tuple<FMatrix, int[]> FBot: FBottoms.getValue()){
                System.out.println(Arrays.toString(FBot.Item2) + ":" + FBot.Item1.hasEmptyR);
                System.out.println(FBot.Item1.toString());
            }
        }
    }

    public FMatrix addFBottom(NetNode parent, int mx, boolean emptyR, int[] splittingIndex){
        FMatrix fBot = new FMatrix(mx, emptyR);
        Tuple<FMatrix, int[]> matrixTuple = new Tuple<FMatrix, int[]>(fBot, splittingIndex);
        List<Tuple<FMatrix, int[]>> fBotList = FBottomsMap.get(parent);
        if(fBotList==null){
            fBotList = new ArrayList<Tuple<FMatrix, int[]>>();
            FBottomsMap.put(parent, fBotList);
        }
        fBotList.add(matrixTuple);
        maxMX = Math.max(mx, maxMX);
        return fBot;
    }


    public void addFBottom(NetNode parent, FMatrix fBot, int[] splittingIndex){
        maxMX = Math.max(fBot.mx, maxMX);
        Tuple<FMatrix, int[]> matrixTuple = new Tuple<FMatrix, int[]>(fBot, splittingIndex);
        List<Tuple<FMatrix, int[]>> fBotList = FBottomsMap.get(parent);
        if(fBotList==null){
            fBotList = new ArrayList<Tuple<FMatrix, int[]>>();
            FBottomsMap.put(parent, fBotList);
        }
        fBotList.add(matrixTuple);
    }

    public void addFTop(NetNode parent, FMatrix fTop, int[] splittingIndex){
        Tuple<FMatrix, int[]> matrixTuple = new Tuple<FMatrix, int[]>(fTop, splittingIndex);
        List<Tuple<FMatrix, int[]>> fTopList = FTopsMap.get(parent);
        if(fTopList==null){
            fTopList = new ArrayList<Tuple<FMatrix, int[]>>();
            FTopsMap.put(parent, fTopList);
        }
        fTopList.add(matrixTuple);
    }

    public void cleanFTop(){
        FTopsMap.clear();
    }

    public void cleanFBottom(){
        FBottomsMap.clear();
    }


    public List<Tuple<FMatrix, int[]>> getFTops(NetNode parent){
        return FTopsMap.get(parent);
    }


    public List<Tuple<FMatrix, int[]>> getFBottoms(NetNode parent){
        return FBottomsMap.get(parent);
    }

    public void cleanFBottomSplittingIndices(){
        NetNode parent = FBottomsMap.keySet().iterator().next();
        double[] arr = new double[R.getMatrixSize(maxMX)];
        int splittingIndexDimension = -1;
        for(Tuple<FMatrix, int[]> fm: FBottomsMap.get(parent)){
            if(splittingIndexDimension == -1){
                splittingIndexDimension = fm.Item2.length;
            }
            int i = 0;
            for(double value: fm.Item1.arr){
                arr[i++] += value;
            }
        }
        int[] splittingIndex = new int[splittingIndexDimension];
        FMatrix fBot = new FMatrix(maxMX, arr, false);
        FBottomsMap.clear();
        addFBottom(parent, fBot, splittingIndex);
    }
}
