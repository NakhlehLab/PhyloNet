package edu.rice.cs.bioinfo.programs.phylonet.algos.treeAugment;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferNetworkMLFromGT;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;
/*
 *@ClassName: inferMPL
 *@Description
 *@Author: Zhen Cao
 *@Date:  2019-06-12 14:15
 *@Version: 1.0
 */
public class inferMPL {
    protected HashMap<String, List<String>> _taxonMap = null;
    protected double _bootstrap = 100;
    protected int _maxReticulations = 5;
    protected long _maxExaminations = -1;
    protected int _maxFailure = 100;
    protected int _moveDiameter = -1;
    protected int _reticulationDiameter = -1;
    protected int _returnNetworks = 10;
    protected int _maxRounds = 100;
    protected int _maxTryPerBranch = 100;
    protected double _maxBranchLength = 6;
    protected double _improvementThreshold = 0.01;
    protected double _Brent1 = 0.01;
    protected double _Brent2 = 0.001;
    protected boolean _dentroscropeOutput = false;
    protected int _parallel = 1;
    protected Set<String> _fixedHybrid = new HashSet<String>();
    private double[] _operationWeight = {0.1,0.05,0.15,0.15,0.15,0.15,2.8}; //p1
//    private double[] _operationWeight = {0.1,0.1,0.15,0.55,0.15,0.15,2.8}; //p2
//    private double[] _operationWeight = {0.5,0,0,0.5,0,0,0};

    //    protected double[] _operationWeight = {0.05,0.05,0.25,0.25,0.2,0.2,2};
    protected boolean _postOptimization = true;
    protected boolean _optimizeBL = false;
    protected int _numRuns = 20;
    protected Long _seed = null;
    protected HashMap<String, Network> _trueNetworks = new HashMap<>();


    protected List initGeneTrees(String gtPath, String outgroup){
        List<MutableTuple<Tree,Double>> geneTrees = new ArrayList<MutableTuple<Tree, Double>>();
        List<List<MutableTuple<Tree,Double>> > geneTreeList = new ArrayList<List<MutableTuple<Tree,Double>>> ();
        String[] gtStrings = treeRoot.FileToString(gtPath).split("\n");
        STITree gtree  = null;
        MutableTuple<Tree,Double> t = null;
        for(String gtString : gtStrings){
            try{
                gtree = treeRoot.removeOutgroup(gtString, outgroup);
                geneTrees.add(new MutableTuple<Tree,Double>(gtree, 1.0));
                geneTreeList.add(geneTrees);
                geneTrees = new ArrayList<MutableTuple<Tree, Double>>();
            }
            catch (Exception e){
                System.err.println("initGeneTrees():"+e);
            }
        }
        return geneTreeList;

    }

    protected List initGeneTrees(String gtPath, String outgroup1, String outgroup2){
//        String geneFolderPath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rootedgT";
//        ArrayList<String> filenames = treeRoot.getFiles(geneFolderPath);
        List<MutableTuple<Tree,Double>> geneTrees = new ArrayList<MutableTuple<Tree, Double>>();
        List<List<MutableTuple<Tree,Double>> > geneTreeList = new ArrayList<List<MutableTuple<Tree,Double>>> ();
        String[] gtStrings = treeRoot.FileToString(gtPath).split("\n");
        STITree gtree  = null;
//        String[] outNodeNames = new String[2] {"Z_0, Z_1"};
        MutableTuple<Tree,Double> t = null;
        for(String gtString : gtStrings){
            try{
                gtree.rerootTreeAtEdge(outgroup1);
                gtree = treeRoot.removeOutgroup(gtString, outgroup1);
                gtree = treeRoot.removeOutgroup(gtString, outgroup2);
                geneTrees.add(new MutableTuple<Tree,Double>(gtree, 1.0));
                geneTreeList.add(geneTrees);
                geneTrees = new ArrayList<MutableTuple<Tree, Double>>();
//                geneTrees = new ArrayList<MutableTuple<Tree, Double>>();
            }
            catch (Exception e){
                System.err.println("initGeneTrees():"+e);
            }
        }
        return geneTreeList;

    }

    protected Network initSpeciesTree(String stPath, String outgroup){
        Network speciesNetwork = null;
        try{
            String stString = treeRoot.FileToString(stPath);
            STITree st = treeRoot.rootTree(stString, outgroup);
//            speciesNetwork = Networks.readNetwork(stString);
            speciesNetwork = Networks.readNetwork(treeRoot.removeOutgroup(st.toString(), outgroup).toString());

//            treeRoot.removeOutgroup(speciesNetwork, "Z");
//            if(_printDetails){
                System.out.println(speciesNetwork.toString());
//            }



        }
        catch (Exception e){
            System.err.println("initSpeciesTree():"+e);
        }
        return speciesNetwork;
//        String networkspath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/networks.json";
//        Map<String, String> networks = treeCompare.GetNetworks(networkspath);
//        String treename = "Reti0_A";
//        _speciesNetwork = Networks.readNetwork(networks.get(treename));
//        _speciesNetwork = Networks.readNetwork("(Z:1.0,(((A:1.0,B:1.0):1.0,(K:1.0,(M:1.0,N:1.0):1.0):1.0):1.0,(((E:1.0,F:1.0):1.0,(G:1.0,H:1.0):1.0):1.0,(C:1.0,(D:1.0,((I:1.0,J:1.0):1.0,(L:1.0,(O:1.0,P:1.0):1.0):1.0):1.0):1.0):1.0):1.0):1.0);");
    }

    protected void initAlleleSpecies(String mapfile, String outgroup){

        String[] lines = treeRoot.FileToString(mapfile).split("\n");
//        _alleles2species = new HashMap<String, String>();
        _taxonMap = new HashMap<String,  List<String>>();
        ArrayList<String> alleles = new ArrayList<String>();
        for(String line : lines){
//            System.out.println(line);
            try {
                String[] items = line.split(":");
                String[] individuals = items[1].split(",");
                if(items[0].equals(outgroup)){
                    continue;
                }
                for(int i = 0; i < individuals.length; i++){
//                    _alleles2species.put(individuals[i], items[0]);
                    alleles.add(individuals[i]);
                }

                _taxonMap.put(items[0], alleles);
                alleles = new ArrayList<String>();
            }catch (Exception e){

            }
        }

    }

//    protected void reorganizeStBrl(Network speciesNetwork){
//        for(Object o: speciesNetwork.dfs()){
//            NetNode node = (NetNode) o;
//            for(Object p: node.getParents()){
//                NetNode parent = (NetNode) p;
//                node.setParentDistance(parent,1.0);
//                if(node.isNetworkNode()){
//                    node.setParentProbability(parent, 0.5);
//                }
//            }
//        }
//
//
//    }

    protected void reorganizeStBrl(Network speciesNetwork){
        for(Object o: speciesNetwork.dfs()){
            NetNode node = (NetNode) o;
            for(Object p: node.getParents()){
                NetNode parent = (NetNode) p;
                if(node.getParentDistance(parent) == NetNode.NO_DISTANCE && !parent.isRoot()){

                    node.setParentDistance(parent,1.0);
                }

                if(node.isNetworkNode()){
                    node.setParentProbability(parent, 0.5);
                }
            }
        }


    }
    public void inferMPL(String gtPath, String stPath, String mapfile){
    }

    public void setParallel(int parallel){
        _parallel = parallel;
    }

//    public String produceResult(String gtPath, String stPath, String mapfile, String outpath, String outgroup){
//        List gts = initGeneTrees(gtPath, outgroup);
//        Network speciesNetwork = initSpeciesTree(stPath);
//        reorganizeStBrl(speciesNetwork);
//        initAlleleSpecies(mapfile);
//        File intermediateResultFile = new File(outpath);
//        InferNetworkMLFromGT inference = new treeAugmentInferNetworkMPL(intermediateResultFile);
////        int [] acceptCount = new int[6];
//        inference.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength, _Brent1, _Brent2, _maxExaminations, _maxFailure, _moveDiameter, _reticulationDiameter, _parallel, speciesNetwork, _fixedHybrid, _operationWeight, _numRuns, _optimizeBL, _seed);
//        LinkedList<Tuple<Network, Double>> resultTuples = new LinkedList<>();
//        inference.inferNetwork(gts,_taxonMap, _maxReticulations, _returnNetworks, _postOptimization, resultTuples);
//        int [] acceptCnt = ((treeAugmentInferNetworkMPL) inference).getAcceptCount();
//        double [] acceptRt = calculateAcceptRate(acceptCnt);
//        int index = 1;
//        StringBuffer result = new StringBuffer();
//        for(Tuple<Network, Double> tuple: resultTuples){
//            result.append("Inferred Network #" + index++ + ":");
//
//            Network n = tuple.Item1;
//
//            for(Object node : n.bfs())
//            {
//                NetNode netNode = (NetNode)node;
//                if(!netNode.isLeaf())
//                {
//                    netNode.setName(NetNode.NO_NAME);
//                }
//            }
//
//            result.append("\n" + n.toString());
//            result.append("\n" + "Total log probability: " + tuple.Item2);
//            result.append("\n Acceptance Rate for each move:");
//            for (double rate:acceptRt) {
//                result.append(rate+",");
//            }
//            result.append("\n");
//
//            if(_dentroscropeOutput){
//                edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.removeAllParameters(n);
//                result.append("\nVisualize in Dendroscope : " + n.toString());
//            }
//        }
//
//        return result.toString();
//
//    }


    //with mode to penalize retinum
    public String produceResult(String gtPath, String stPath, String mapfile, String outgroupG, String outgroupS, int mode){
        long startingTime = System.currentTimeMillis();
        List gts = initGeneTrees(gtPath, outgroupG);
        Network speciesNetwork = initSpeciesTree(stPath, outgroupS);
        reorganizeStBrl(speciesNetwork);
        initAlleleSpecies(mapfile, outgroupS);

        InferNetworkMLFromGT inference = new treeAugmentInferNetworkMPL();
//        int [] acceptCount = new int[6];
        inference.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength, _Brent1, _Brent2, _maxExaminations, _maxFailure, _moveDiameter, _reticulationDiameter, _parallel, speciesNetwork, _fixedHybrid, _operationWeight, _numRuns, _optimizeBL, _seed);
        LinkedList<Tuple<Network, Double>> resultTuples = new LinkedList<>();
        ((treeAugmentInferNetworkMPL) inference).inferNetwork(gts,_taxonMap, _maxReticulations, _returnNetworks, _postOptimization, resultTuples, mode);
        int [] acceptCnt = ((treeAugmentInferNetworkMPL) inference).getAcceptCount();
        double [] acceptRt = calculateAcceptRate(acceptCnt);
        int index = 1;
        StringBuffer result = new StringBuffer();
        for(Tuple<Network, Double> tuple: resultTuples){
            result.append("Inferred Network #" + index++ + ":");

            Network n = tuple.Item1;

            for(Object node : n.bfs())
            {
                NetNode netNode = (NetNode)node;
                if(!netNode.isLeaf())
                {
                    netNode.setName(NetNode.NO_NAME);
                }
            }

            result.append("\n" + n.toString());
            result.append("\n" + "Total log probability: " + tuple.Item2);
            result.append("\n Acceptance Rate for each move:");
            for (double rate:acceptRt) {
                result.append(rate+",");
            }
            result.append("\n");

            if(_dentroscropeOutput){
                edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.removeAllParameters(n);
                result.append("\nVisualize in Dendroscope : " + n.toString());
            }
        }
        double time = (System.currentTimeMillis()-startingTime)/60000.0;
        result.append("Running time:"+time);
        return result.toString();


    }

    public double[] calculateAcceptRate(int[] acceptCnt){
        double [] acceptRt = new double[acceptCnt.length];
        int sum = 0;
        for(int i = 0; i < acceptCnt.length; i++){
            sum += acceptCnt[i];
        }
        for(int i = 0; i < acceptCnt.length; i++){
            acceptRt[i] = acceptCnt[i]*1.0/sum;
        }
        return acceptRt;
    }

    //Here we set the max reticulation number to be the one we know
    public void setMaxReticulations(String treeName){
        _maxReticulations = Integer.parseInt(treeName.substring(4, 5));
    }

    public void initTrueNetworks(String trueNetworksPath, String outgroup){
//        String trueNetworksPath = path + "/networks.json";
        Map<String, String> trueNetworkst = treeCompare.GetNetworks(trueNetworksPath);
//        Map<String, Network> trueNetworks = new HashMap<>();
        for(Map.Entry<String, String> entry : trueNetworkst.entrySet()) {
            Network net = Networks.readNetwork(entry.getValue());
            if(!outgroup.equals("")){
                treeRoot.removeOutgroup(net, outgroup);
            }

//            System.out.println(net.toString());
            _trueNetworks.put(entry.getKey(), net);
        }
//        return trueNetworks;
    }
    /** Here is for the subnetwork
     * */
    public Network getTrueSubNetwork(Network trueNetwork, String selectedLeaves){
        List<String> leafList = new ArrayList<>();
        for(int i = 0; i < selectedLeaves.length(); i++){
            leafList.add(selectedLeaves.substring(i, i+1));
        }
        SuperNetwork superNetwork = new SuperNetwork(new ArrayList<>());
        System.out.println("True network: " + trueNetwork.toString());

        superNetwork.setTrueNetwork(trueNetwork);
        return superNetwork.getSubNetwork(trueNetwork, leafList, true, false);
    }


    /* Set the maxreti of sub to be the true subnetwork reti
    * */
    public void setMaxReticulations(String truepath, String treeName, String trinetname){
        initTrueNetworks(truepath, "");
        Network trueNet = _trueNetworks.get(treeName);
        Network subNet = getTrueSubNetwork(trueNet, trinetname);
        _maxReticulations = subNet.getReticulationCount();
        System.out.println(_maxReticulations);
        System.out.println("true subnet is "+subNet.toString());
        System.out.println("_maxReticulations for this subnetwork "+trinetname+" is "+_maxReticulations);
    }

    public void setMaxReticulations(int num){
        _maxReticulations = num;
    }

//    public String produceResultRetiNum(String treeName, String gtpath, String stpath, String mapfile, String outpath){
//        setMaxReticulations(treeName);
//        return produceResult(gtpath, stpath, mapfile, outpath);
//    }






//    public static void RunOnceLocalMPL(){
//        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/IQTreeAndAstralTree0225/";
//        String treename = "Reti1_A";
////        Test t = new Test();
//
////        int maxreticulation = 5;
//        if(treename.startsWith("R")) {
//            String stpath = path+treename + "/rooted_astral_tre.tre";
//            String gtpath = path+treename  + "/rgts.tre";
//            String outpath = path+treename  + "/MPL.txt";
//            String mapfile = path + "/iqmap.txt";
//            inferMPL mplinfer = new inferMPL();
//            String resultString = mplinfer.produceResult(gtpath, stpath, mapfile, outpath);
//            System.out.println(resultString);
//            String filepath = path+treename  + "/MPLresult.txt";
//            treeRoot.StringToFile(resultString, filepath);
//        }
//    }

    public static void RunOnceLocalMPLMode(int mode){
        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/treeAugment/data/rainbow_skink_data/";
        String treename = "Reti1_A";

        if(treename.startsWith("R")) {
            String stpath = path+treename + "/rooted_astral_tre_0.tre";
            String gtpath = path+treename  + "/rgts_0.tre";
            String outpath = path+treename  + "/MPL.txt";
            String mapfile = path + "/iqmap.txt";
            inferMPL mplinfer = new inferMPL();
            String resultString = mplinfer.produceResult(gtpath, stpath, mapfile, "Z_0", "Z", mode);
            System.out.println(resultString);
            String filepath = path+treename  + "/MPLresult.txt";
            treeRoot.StringToFile(resultString, filepath);
        }
    }

//    public static void runOnceNots(String path, String treename){
////        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/IQTreeAndAstralTree/";
////        String treename = "Reti0_A";
////        Test t = new Test();
//
////        int maxreticulation = 5;
//        if(treename.startsWith("R")) {
//            String stpath = path + treename + "/rooted_astral_tre.tre";
//            String gtpath = path + treename + "/rgts.tre";
//            String outpath = path + treename + "/MPL.txt";
//            String mapfile = path + "/iqmap.txt";
//            inferMPL mplinfer = new inferMPL();
//            String resultString = mplinfer.produceResult(gtpath, stpath, mapfile, outpath);
////            resultString += ;
//            System.out.println(resultString);
//            String filepath = path + treename + "/MPLresult.txt";
//            treeRoot.StringToFile(resultString, filepath);
//        }
//    }

//    public static void runOnceNotsReti(String path, String treename){
//        if(treename.startsWith("R")) {
//            String stpath = path + treename + "/rooted_astral_tre.tre";
//            String gtpath = path + treename + "/rgts.tre";
//            String outpath = path + treename + "/MPL.txt";
//            String mapfile = path + "/iqmap.txt";
//            inferMPL mplinfer = new inferMPL();
//            String resultString = mplinfer.produceResultRetiNum(treename, gtpath, stpath, mapfile, outpath);
////            resultString += ;
//            System.out.println(resultString);
//            String filepath = path + treename + "/MPLresult.txt";
//            treeRoot.StringToFile(resultString, filepath);
//        }
//    }

    //mode=0: normal; mode=1: aic; mode=2: bic; mode=3: set maxreti be true
    public static void runOnceNotsMode(String path, String treename, int mode, int parallel, int num){

        if(treename.startsWith("R")) {
            String stpath = path + treename + "/rooted_astral_tre_"+num+".tre";
            String gtpath = path + treename + "/rgts_"+num+".tre";
            String mapfile = path + "/iqmap.txt";
            inferMPL mplinfer = new inferMPL();
            mplinfer.setParallel(parallel);
            String modename = "0";
            if(mode == 1){
                modename = "AIC";
            }
            else if(mode == 2){
                modename = "BIC";
            }
            else if(mode == 3){
                modename = "MAX";
                mplinfer.setMaxReticulations(treename);
            }
            else if(mode == 4){
                modename = "4";
                mplinfer.setMaxReticulations(6);
            }

            String resultString = mplinfer.produceResult(gtpath, stpath, mapfile, "Z_0", "Z", mode);
            System.out.println(resultString);
            String filepath = path + treename + "/MPL"+ modename +"_"+num+"_result.txt";
            treeRoot.StringToFile(resultString, filepath);
        }
    }

    public static void runOnceLocalModeEmpirical(){
        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/treeAugment/data/rainbow_skink_data/";
        String stpath = path  + "/lizard_concatRST.txt";
        String gtpath = path  + "/rgts.tre";
        String mapfile = path + "/map.txt";
        int mode = 3;
        inferMPL mplinfer = new inferMPL();
        mplinfer.setParallel(2);
        String modename = "MAX";
        mplinfer.setMaxReticulations(1);

        String resultString = mplinfer.produceResult(gtpath, stpath, mapfile, "SP07_indexing28_h0", "SP07indexing28", mode);
        System.out.println(resultString);
        String filepath = path + "/concatMPL"+ modename +"_result.txt";
        treeRoot.StringToFile(resultString, filepath);

    }

    public static void runOnceNotsSub(String path, String treename, String trinetname, int mode, int parallel, int max){

        if(treename.startsWith("R")) {
            String stpath = path + treename+"/subgt/" + trinetname + "/rooted_astral_tre.tre";
            String gtpath = path + treename + "/subgt/" + trinetname + "/rgts.txt";
            String truepath = path + "/networks.json";
            String mapfile = path + treename + "/subgt/" + trinetname + "/map.txt";
            inferMPL mplinfer = new inferMPL();
            mplinfer.setParallel(parallel);
            String modename = "0";
            if(mode == 1){
                modename = "AIC";
            }
            else if(mode == 2){
                modename = "BIC";
            }
            else if(mode == 3){
                modename = "MAX";
                mplinfer.setMaxReticulations(truepath, treename, trinetname);
                System.out.println(treename+" : "+trinetname);
            }
            else if(mode == 4){
                modename = "4";
                mplinfer.setMaxReticulations(max);
            }

            String resultString = mplinfer.produceResult(gtpath, stpath, mapfile, "Z0", "Z", mode);
            System.out.println(resultString);
            String filepath = path + treename+"/subgt/" + trinetname + "/MPL"+ modename +"_result.txt";
            treeRoot.StringToFile(resultString, filepath);
        }
    }

    public static void runOnceNotsSub(String path, String treename, String trinetname, int mode, int parallel, int max, int num){

        if(treename.startsWith("R")) {
            String stpath = path + treename+"/subgt/" + trinetname + "/st"+num+".txt";
            String gtpath = path + treename + "/subgt/" + trinetname + "/rgts.txt";
            String truepath = path + "/networks.json";
            String mapfile = path + treename + "/subgt/" + trinetname + "/map.txt";
            inferMPL mplinfer = new inferMPL();
            mplinfer.setParallel(parallel);
            String modename = "0";
            if(mode == 1){
                modename = "AIC";
            }
            else if(mode == 2){
                modename = "BIC";
            }
            else if(mode == 3){
                modename = "MAX";
                mplinfer.setMaxReticulations(truepath, treename, trinetname);
                System.out.println(treename+" : "+trinetname);
            }
            else if(mode == 4){
                modename = "4";
                mplinfer.setMaxReticulations(max);
            }

            String resultString = mplinfer.produceResult(gtpath, stpath, mapfile, "", "", mode);
            System.out.println(resultString);
            String filepath = path + treename+"/subgt/" + trinetname + "/MPL"+ modename+num +"_result.txt";
            treeRoot.StringToFile(resultString, filepath);
        }
    }

    public static void runOnceLocalSub(String trinetname){
        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/IQTreeSub3/";
        String treename = "Reti4_D";
//        String trinetname = "ABCZ";
        int mode = 3;
        int parallel = 1;
        int max = 4;
        if(treename.startsWith("R")) {

            String stpath = path + treename+"/subgt/" + trinetname + "/st1.txt";
            String gtpath = path + treename + "/subgt/" + trinetname + "/rgts.txt";
            String truepath = path + "/networks.json";
            String mapfile = path + treename + "/subgt/" + trinetname + "/map.txt";
            inferMPL mplinfer = new inferMPL();
            mplinfer.setParallel(parallel);
            String modename = "0";
            if(mode == 1){
                modename = "AIC";
            }
            else if(mode == 2){
                modename = "BIC";
            }
            else if(mode == 3){
                modename = "MAX";
                mplinfer.setMaxReticulations(truepath, treename, trinetname);
            }
            else if(mode == 4){
                modename = "4";
                mplinfer.setMaxReticulations(max);
            }

            String resultString = mplinfer.produceResult(gtpath, stpath, mapfile, "", "", mode);
            System.out.println(resultString);
            String filepath = path + treename+"/subgt/" + trinetname + "/MPL"+ modename +"_result.txt";
            treeRoot.StringToFile(resultString, filepath);
        }
    }


    public static void runOneNetLocalSub(){
        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/IQTreeSub/";
        String treename = "Reti4_D";
//        String trinetname = "AEPZ";
//        int mode = 1;
//        int parallel = 1;
        ArrayList<String> treenames = treeRoot.getFolers(path+treename+"/subgt/");
        Collections.sort(treenames);
        for(String subname: treenames){
            String trinetname = subname.substring(subname.length()-4);
            System.out.println(trinetname);
            runOnceLocalSub(trinetname);
        }

    }



    public static void main(String[] args) {

//        RunOnceLocalMPL();
//        RunOnceLocalMPLMaxReti();
//        RunOnceLocalMPLMode(0);
        runOnceLocalModeEmpirical();

//        runOnceNots(args[0], args[1]);
//        runOnceNotsReti(args[0], args[1]);
//        runOnceNotsMode(args[0], args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]));
//        runOnceNotsMode(args[0], args[1], 1);
//        runOnceNotsSub(args[0], args[1], args[2], Integer.parseInt(args[3]), Integer.parseInt(args[4]), Integer.parseInt(args[5]), Integer.parseInt(args[6]));
//        runOnceLocalSub("ABZ");
//        runOneNetLocalSub();
    }
}
