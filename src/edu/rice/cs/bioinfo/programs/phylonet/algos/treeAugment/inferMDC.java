package edu.rice.cs.bioinfo.programs.phylonet.algos.treeAugment;
/*
 *@ClassName: inferMDC
 *@Description
 *@Author: Zhen Cao
 *@Date:  2019-06-12 14:15
 *@Version: 1.0
 */

import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;



public class inferMDC {
    private HashMap<String, List<String>> _species2alleles = null;
    private HashMap<String, String> _allele2species = null;
    //    private List<List<NetworkNonEmpty>> _geneTrees;
    private double _bootstrap = 100;
    private NetworkNonEmpty _startSpeciesNetwork = null;
    private int _maxReticulations = 5;
    private long _maxExaminations = -1;
    private int _moveDiameter = -1;
    private int _reticulationDiameter = -1;
    private int _maxFailure = 100;
    private int _returnNetworks = 1;
    private int _numProcessors = 1;
    private boolean _dentroscropeOutput = false;
    private int _numRuns = 20;
    private double[] _operationWeight = {0.1,0.1,0.15,0.55,0.15,0.15, 2.8};
    private Long _seed = null;

    private Set<String> _fixedHybrid = new HashSet<String>();
    protected HashMap<String, Network> _trueNetworks = new HashMap<>();


    protected List initGeneTrees(String gtPath, String outgroup){
//        String geneFolderPath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rootedgT";
//        ArrayList<String> filenames = treeRoot.getFiles(geneFolderPath);
        List<MutableTuple<Tree,Double>> geneTrees = new ArrayList<MutableTuple<Tree, Double>>();
//        List<List<MutableTuple<Tree,Double>> > geneTreeList = new ArrayList<List<MutableTuple<Tree,Double>>> ();
        String[] gtStrings = treeRoot.FileToString(gtPath).split("\n");
        STITree gtree  = null;
        MutableTuple<Tree,Double> t = null;
        for(String gtString : gtStrings){
            try{
//                gtree = new STITree(gtString);

                gtree = treeRoot.removeOutgroup(gtString, outgroup);
                geneTrees.add(new MutableTuple<Tree,Double>(gtree, 1.0));
//                geneTreeList.add(geneTrees);
//                geneTrees = new ArrayList<MutableTuple<Tree, Double>>();
//                geneTrees = new ArrayList<MutableTuple<Tree, Double>>();
            }
            catch (Exception e){
                System.err.println("initGeneTrees():"+e);
            }
        }
//        return geneTreeList;
        return geneTrees;
    }

    protected Network initSpeciesTree(String stPath, String outgroup){
        Network speciesNetwork = null;
        try{
            String stString = treeRoot.FileToString(stPath);
//            speciesNetwork = Networks.readNetwork(stString);
//            treeRoot.removeOutgroup(speciesNetwork, "Z")
            speciesNetwork = Networks.readNetwork(treeRoot.removeOutgroup(stString, outgroup).toString());
//            if(_printDetails){
//                System.out.println(_speciesNetwork.toString());
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
        _allele2species = new HashMap<String, String>();
        _species2alleles = new HashMap<String, List<String>>();
//        _taxonMap = new HashMap<String,  List<String>>();
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
                    _allele2species.put(individuals[i], items[0]);
                    alleles.add(individuals[i]);
                }

                _species2alleles.put(items[0], alleles);
                alleles = new ArrayList<String>();
            }catch (Exception e){

            }
        }

    }
    public void setParallel(int parallel){
        _numProcessors = parallel;
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
            treeRoot.removeOutgroup(net, outgroup);
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
        System.out.println("true subnet is "+subNet.toString());
        System.out.println("_maxReticulations for this subnetwork "+trinetname+" is "+_maxReticulations);
    }

    public String produceResult(String gtPath, String stPath, String mapfile, String outgroupG, String outgroupS){
        long startingTime = System.currentTimeMillis();
        List gts = initGeneTrees(gtPath, outgroupG);
        Network speciesNetwork = initSpeciesTree(stPath, outgroupS);
        reorganizeStBrl(speciesNetwork);
        initAlleleSpecies(mapfile, outgroupS);
//        File intermediateResultFile = new File(outpath);

        treeAugmentInferNetworkMDC inference = new treeAugmentInferNetworkMDC();
        inference.setSearchParameter(_maxExaminations, _maxFailure, _moveDiameter, _reticulationDiameter, speciesNetwork, _fixedHybrid, _numProcessors, _operationWeight, _numRuns, _seed);
//        inference.setIntermediateResultFile(intermediateResultFile);
        List<Tuple<Network, Double>> resultTuples = inference.inferNetwork(gts,_species2alleles,_maxReticulations, _returnNetworks);
        //System.out.print(System.currentTimeMillis()-start);
        StringBuffer result = new StringBuffer();
        int index = 1;
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
            result.append("\n" + "Total number of extra lineages: " + tuple.Item2);
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

    public static void RunOnceLocalMDC(){
        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/IQTreeAndAstralTree0313/";
        String treename = "Reti4_A";
//        Test t = new Test();
//        int maxreticulation = 5;
        if(treename.startsWith("R")) {
            String stpath = path+treename + "/rooted_astral_tre.tre";
            String gtpath = path+treename  + "/rgts.tre";
            String outpath = path+treename  + "/MDC.txt";
            String mapfile = path + "/iqmap.txt";
            inferMDC mdcinfer = new inferMDC();
            mdcinfer.setMaxReticulations(treename);
            String resultString = mdcinfer.produceResult(gtpath, stpath, mapfile, "Z_0", "Z");
            System.out.println(resultString);
            String filepath = path+treename  + "/MDCresult.txt";
            treeRoot.StringToFile(resultString, filepath);
        }
    }

    public static void runOnceNotsReti(String path, String treename){
        if(treename.startsWith("R")) {
            String stpath = path + treename + "/rooted_astral_tre.tre";
            String gtpath = path + treename + "/rgts.tre";
            String outpath = path + treename + "/MDC.txt";
            String mapfile = path + "/iqmap.txt";
            inferMDC mdcinfer = new inferMDC();
            mdcinfer.setMaxReticulations(treename);
            String resultString = mdcinfer.produceResult(gtpath, stpath, mapfile, "Z_0", "Z");
            System.out.println(resultString);
            String filepath = path + treename + "/MDCresult.txt";
            treeRoot.StringToFile(resultString, filepath);
        }
    }

    //mode=1, max reticulation number=true number;  mode =0, no a prior
    public static void runOnceNotsM(String path, String treename, int mode, int parallel){
        if(treename.startsWith("R")) {
            String stpath = path + treename + "/rooted_astral_tre.tre";
            String gtpath = path + treename + "/rgts.tre";
            String outpath = path + treename + "/MDC.txt";
            String mapfile = path + "/iqmap.txt";
            inferMDC mdcinfer = new inferMDC();
            mdcinfer.setParallel(parallel);
            if(mode == 1){
                mdcinfer.setMaxReticulations(treename);
            }
            String resultString = mdcinfer.produceResult(gtpath, stpath, mapfile, "Z_0", "Z");
            System.out.println(resultString);
            String filepath = path + treename + "/MDC"+mode+"result.txt";
            treeRoot.StringToFile(resultString, filepath);
        }
    }

    public static void runOnceNotsSub(String path, String treename, String trinetname, int mode, int parallel, int num){

        if(treename.startsWith("R")) {
            String stpath = path + treename+"/subgt/" + trinetname + "/st"+num+".txt";
            String gtpath = path + treename + "/subgt/" + trinetname + "/rgts.txt";
            String truepath = path + "/networks.json";
            String mapfile = path + treename + "/subgt/" + trinetname + "/map.txt";
            inferMDC mdcinfer = new inferMDC();
            mdcinfer.setParallel(parallel);
            String modename = "0";
            if(mode == 1){
                modename = "MAX";
                mdcinfer.setMaxReticulations(truepath, treename, trinetname);
            }

            String resultString = mdcinfer.produceResult(gtpath, stpath, mapfile, "", "");
            System.out.println(resultString);
            String filepath = path + treename+"/subgt/" + trinetname + "/MDC"+ modename + num +"_result.txt";
            treeRoot.StringToFile(resultString, filepath);
        }
    }

    public static void runOnceLocalSub(String trinetname){
        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/IQTreeSub3/";
        String treename = "Reti5_D";
//        String trinetname = "AEPZ";
        int mode = 1;
        int parallel = 1;
        ArrayList<String> treenames = treeRoot.getFolers(path);
        Collections.sort(treenames);
        if(treename.startsWith("R")) {

            String stpath = path + treename+"/subgt/" + trinetname + "/st1.txt";
            String gtpath = path + treename + "/subgt/" + trinetname + "/rgts.txt";
            String truepath = path + "/networks.json";
            String mapfile = path + treename + "/subgt/" + trinetname + "/map.txt";
            inferMDC mdcinfer = new inferMDC();
            mdcinfer.setParallel(parallel);
            String modename = "0";
            if(mode == 1){
                modename = "MAX";
                mdcinfer.setMaxReticulations(truepath, treename, trinetname);
            }


            String resultString = mdcinfer.produceResult(gtpath, stpath, mapfile, "", "");
            System.out.println(resultString);
            String filepath = path + treename+"/subgt/" + trinetname + "/MDC"+ modename +"_result.txt";
            treeRoot.StringToFile(resultString, filepath);
        }
    }

    public static void runOneNetLocalSub(){
        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/IQTreeSub/";
        String treename = "Reti4_D";
//        String trinetname = "AEPZ";
        int mode = 1;
        int parallel = 1;
        ArrayList<String> treenames = treeRoot.getFolers(path+treename+"/subgt/");
        Collections.sort(treenames);
        for(String subname: treenames){
            String trinetname = subname.substring(subname.length()-4);
            System.out.println(trinetname);
            runOnceLocalSub(trinetname);
        }

    }

    public static void main(String[] args) {
//        RunOnceLocalMDC();
//        runOnceNotsReti(args[0], args[1], args[3]);
//        runOnceNotsM(args[0], args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]));
//        runOnceLocalSub("ABZ");
        runOnceNotsSub(args[0], args[1], args[2], Integer.parseInt(args[3]), Integer.parseInt(args[4]), Integer.parseInt(args[5]));
//        runOneNetLocalSub();
    }

}
