package edu.rice.cs.bioinfo.programs.phylonet.algos.treeAugment;


import com.google.gson.Gson;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;
/*
 * @author: Zhen Cao
 * @Date: 1/19/2019
 * */
public class treeCompare {

    public static double getRootedRobinsonFouldsDistance(Tree tree1, Tree tree2) {
        SymmetricDifference symmetricDifference = new SymmetricDifference();
        symmetricDifference.computeDifference(tree1, tree2, true);
        double diff = symmetricDifference.getWeightedAverage();
        double undiff = symmetricDifference.getUnweightedAverage();
        return undiff;
    }

    public static double getRFDistance(String tree, String net){
//        System.out.println(tree+"\n"+net+"\n*****");
        double RF = 20;
        String minN="";

        try{
            Network<String> network = Networks.readNetwork(net);
            for (NetworkTree<String> nt: Networks.getTrees(network)) {
                Tree t = Trees.readTree(tree);
                Tree n = Trees.readTree(nt.makeTree().toString());
                double tmpRF = getRootedRobinsonFouldsDistance(t, n);
                if(RF > tmpRF){
                    RF = tmpRF;
                    minN = n.toString();
                }
            }
            if(RF > 0){
                System.out.println(RF);
                System.out.println("tree:\t"+tree);
                System.out.println("network:\t"+net);
                System.out.println("induced:\t"+minN);
            }


        }
        catch (Exception e){
            System.err.println("rootTree():"+e);
            System.err.println(RF+":"+tree);
        }
//        System.out.println(RF);
        return RF;


    }



    public static class jsonNets{

        private List<jsonNet> networks;

        public jsonNets(){

        }
        public List<jsonNet> getNetworks(){
            return networks;
        }

        public void setNetworks(List<jsonNet> networks){
            this.networks = networks;
        }

        public static class jsonNet{
            private String tag;
            private String netstring;

            public jsonNet(){

            }
            public String getTag(){
                return tag;
            }

            public void setTag(String tag){
                this.tag = tag;
            }


            public String getNetstring(){
                return netstring;
            }

            public void setNetstring(String netstring){
                this.netstring = netstring;
            }

        }
    }

    public static Map<String, String> GetNetworks(String networkspath){
        String json = treeRoot.FileToString(networkspath);
        Gson gson = new Gson();
        jsonNets j = gson.fromJson(json,jsonNets.class);
        List<jsonNets.jsonNet> networks=j.getNetworks();
        Map m = new HashMap<String, String>();
        for(jsonNets.jsonNet jnet : networks){
            m.put(jnet.getTag(), jnet.getNetstring());
        }
        return m;
    }

    public static void writeNetworkToFile(Map<String, String> networks, String dir){
        for(Map.Entry<String, String> entry : networks.entrySet()){
            treeRoot.StringToFile(entry.getValue(), dir+entry.getKey()+".txt");
        }
    }

    public static ArrayList<String> runASTRAL(){
        String treefolderpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rootedSpeciesTree_Astral";
        String networkspath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/networks.json";

        Map<String, String> networks = GetNetworks(networkspath);


        String networkOutPath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/networks/";
        writeNetworkToFile(networks, networkOutPath);


        ArrayList<String> treefilepaths = treeRoot.getFiles(treefolderpath);
        int countNotInducedTree = 0;
        int countInduced = 0;
        ArrayList<String> RightTreeNameList = new ArrayList<String>();
        for(String treefilename : treefilepaths){

            String[] arr = treefilename.split("/");
            String treefilepath = arr[arr.length - 1];
            double dist = 0;
            try {
                String treeName = treefilepath.substring(0, 5) + "_" + treefilepath.substring(5,6);
                if(!treeName.startsWith("R")){
//                    System.out.println(treeName);
                    continue;
                }
                String tree = treeRoot.FileToString(treefilename);
                dist = getRFDistance(tree, networks.get(treeName));
                if(dist > 0){
                    System.out.println("treeName:\t"+treeName);
                    System.out.println("*****************");
                    countNotInducedTree++;
                    RightTreeNameList.add(treeName);
                }
                else{
                    countInduced++;

                }

            }catch (Exception e) {
                System.err.println(treefilepath);
            }

        }
//        System.out.println("trees:"+treefilepaths.size());
        System.out.println("trees not included in network:"+countNotInducedTree);
        System.out.println("trees included in network:"+countInduced);
        Collections.sort(RightTreeNameList);
        return RightTreeNameList;
    }

    public static ArrayList<String> runRAxML(){
        String treefolderpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rootedSpeciesTree_RAxML";
        String networkspath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/networks.json";

        Map<String, String> networks = GetNetworks(networkspath);
        ArrayList<String> RightTreeNameList = new ArrayList<String>();

//        String networkOutPath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/networks/";
//        writeNetworkToFile(networks, networkOutPath);


        ArrayList<String> treefilepaths = treeRoot.getFiles(treefolderpath);
        int countNotInducedTree = 0;
        int countInduced = 0;
        double distsum = 0;
        for(String treefilename : treefilepaths){
            String[] arr = treefilename.split("/");
            String treefilepath = arr[arr.length - 1];
            double dist = 0;
            try {
                String treeName = treefilepath.substring(0,7);
                if(!treeName.startsWith("R")){
//                    System.out.println(treeName);
                    continue;
                }
                String tree = treeRoot.FileToString(treefilename);
                dist = getRFDistance(tree, networks.get(treeName));
                if(dist > 0){
                    System.out.println("treeName:\t"+treeName+"distance="+dist);

                    System.out.println("*****************");
                    countNotInducedTree++;
                    RightTreeNameList.add(treeName);


                }
                else{
                    countInduced++;

                }

            }catch (Exception e) {
                System.err.println(treefilepath);
            }
            distsum += dist;
        }
        System.out.println("trees not included in network:"+countNotInducedTree);
        System.out.println("trees included in network:"+countInduced);
        Collections.sort(RightTreeNameList);
        System.out.println("distsum="+distsum);
        return RightTreeNameList;
    }



    //todo:root, delete outgroup, map individual to species
    public static ArrayList<String> concatedIQST(){
        ArrayList<String> RightTreeNameList = new ArrayList<>();
        String treefolderpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rootedSpeciesTree_IQTree/";
        String networkspath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/networks.json";

        Map<String, String> networks = GetNetworks(networkspath);

        ArrayList<String> treefilepaths = treeRoot.getFiles(treefolderpath);
        int countNotInducedTree = 0;
        int countInduced = 0;
        double distsum = 0;
        for(String treefilename : treefilepaths){
            String[] arr = treefilename.split("/");
            String treefilepath = arr[arr.length - 1];
            double dist = 0;
            try {
                String treeName = treefilepath.substring(0,7);
                if(!treeName.startsWith("R")){
//                    System.out.println(treeName);
                    continue;
                }
                String tree = treeRoot.FileToString(treefilename);
                dist = getRFDistance(tree, networks.get(treeName));
                if(dist > 0){
                    System.out.println("treeName:\t"+treeName+"distance="+dist);

                    System.out.println("*****************");
                    countNotInducedTree++;
                    RightTreeNameList.add(treeName);


                }
                else{
                    countInduced++;

                }

            }catch (Exception e) {
                System.err.println(treefilepath);
            }
            distsum += dist;
        }
        System.out.println("trees not included in network:"+countNotInducedTree);
        System.out.println("trees included in network:"+countInduced);
        Collections.sort(RightTreeNameList);
        System.out.println("distsum="+distsum);

        return RightTreeNameList;
    }

    public static ArrayList<String> IQASTRAL(){
        String treefolderpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/IQTreeAndAstralTree";
//        String networkspath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/networks.json";
//        String treefolderpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/proj/code/data1000gt/";
        String networkspath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/networks.json";
        boolean rootgenetree = false;
        boolean writespeciestree = false;
        Map<String, String> networks = GetNetworks(networkspath);


        ArrayList<String> treefilepaths = treeRoot.getFolers(treefolderpath);
        System.out.println(treefilepaths);
        int countNotInducedTree = 0;
        int countInduced = 0;
        ArrayList<String> RightTreeNameList = new ArrayList<String>();



        for(String treefilename : treefilepaths){
            String[] arr = treefilename.split("/");
            String treeName = arr[arr.length - 1];
//            System.out.println(treeName);
            double dist = 0;
            String infile = treefilename + "/astral_tree.tre";
            try {
                if(!treeName.startsWith("R")){
                    continue;
                }
                String tree = treeRoot.FileToString(infile);
                String rootedTree = treeRoot.rootTree(tree, "Z", "nodelete");

                if(writespeciestree){
                    String outfilename = treefilename+"/rooted_astral_tre_0.tre";
                    treeRoot.StringToFile(rootedTree, outfilename);
                }


                if (rootgenetree) {
                    String gtfile = treefilename + "/genetrees_0.txt";
                    String[] gt = treeRoot.FileToString(gtfile).split("\n");
                    String gtoutname = treefilename+"/rgts_0.tre";
                    try{
                        treeRoot.StringToFile(treeRoot.rootTrees(gt, "Z_0", "Z_1"), gtoutname);
                    } catch (Exception e){
                        System.err.println(e);
                    }
                }
                dist = getRFDistance(rootedTree, networks.get(treeName));
                if(dist > 0){
                    System.out.println("treeName:\t"+treeName);
                    System.out.println("*****************");
                    countNotInducedTree++;
                    RightTreeNameList.add(treeName);
                }
                else{
                    countInduced++;
                }

            }catch (Exception e) {
                System.err.println();
            }
        }
        System.out.println("trees not included in network:"+countNotInducedTree);
        System.out.println("trees included in network:"+countInduced);
        Collections.sort(RightTreeNameList);
        return RightTreeNameList;
    }

    public static List<Tree> getGeneTrees(String gtPath){
        List<Tree> geneTrees = new ArrayList<Tree>();
        String[] gtStrings = treeRoot.FileToString(gtPath).split("\n");
        STITree gtree  = null;
        MutableTuple<Tree,Double> t = null;
        for(String gtString : gtStrings){
            try{
                gtree = treeRoot.removeOutgroup(gtString, "Z_0");
                geneTrees.add(gtree);
            }
            catch (Exception e){
                System.err.println("initGeneTrees():"+e);
            }
        }
        return geneTrees;
    }

    public static Tree getSpeciesTree(String stPath){
        STITree stree = null;
        try{
            String stString = treeRoot.FileToString(stPath);
            stree = treeRoot.removeOutgroup(stString, "Z");
        }
        catch (Exception e){
            System.err.println("initSpeciesTree():"+e);
        }
        return stree;
    }

    public static HashMap<Tree, Integer> getDifferentGTs(List<Tree> gts){
        HashMap<Tree, Integer> distinctGeneTrees = new HashMap();

        for(Tree gt : gts){
            boolean inDistinctList = false;
            for (Tree ref:distinctGeneTrees.keySet()){
                if(getRootedRobinsonFouldsDistance(gt, ref) == 0){
                    inDistinctList = true;
                    int count  = distinctGeneTrees.get(ref);
                    distinctGeneTrees.remove(ref);
                    distinctGeneTrees.put(ref, count+1);
                    break;
                }
            }
            if(!inDistinctList){
                distinctGeneTrees.put(gt, 1);
            }
//            else{
//                inDistinctList = false;
//            }
        }
        System.out.println("the number of different gene trees: "+ distinctGeneTrees.size());
        return distinctGeneTrees;
    }

    public static void getGTs_ST_Diff(List<Tree> gts, Tree st){
        double aveDiff = 0;
        int sameTopo = 0;
        for(Tree gt: gts){
            double diff = getRootedRobinsonFouldsDistance(gt, st);
            if(diff == 0){
                sameTopo ++;
            }
            aveDiff += diff;
        }
        aveDiff /= gts.size();

        System.out.println("the number of gene trees==species trees: "+sameTopo);
        System.out.println("the average difference with species trees is: "+ aveDiff);
    }

    public static void checkTreeDifference(String path, String treename, int num){
        if(treename.startsWith("R")) {
            String stpath = "";
            String gtpath = "";
            if(num == -1){
                stpath = path + treename + "/rooted_astral_tre.tre";
                gtpath = path + treename + "/rgts.tre";
            }
            else{
                stpath = path + treename + "/rooted_astral_tre_" + num + ".tre";
                gtpath = path + treename + "/rgts_" + num + ".tre";
            }

            List<Tree> gts = getGeneTrees(gtpath);

            HashMap<Tree, Integer> diffGts = getDifferentGTs(gts);
//            Tree speciesTree = getSpeciesTree(stpath);
//            getGTs_ST_Diff(gts, speciesTree);

        }
    }

    public static void TestgetDifferentGTs(){
        String gtString1 = "(Z:6.5032896712070904,((B:6.5032896712070904,A:6.5032896712070904)i20:4.711530201978998,(((M:1.21493846897087,N:1.2015614615369377)i15:0.5266215533418312,K:1.6769779365754285)i16:6.5032896712070904,(((H:2.4173133586554716,G:2.4515047234037515)i10:1.8649265178632008,(F:2.632088660299164,E:3.2451931331855737)i7:1.2856402077264741)i11:4.502775388116889,(C:6.5032896712070904,(D:4.557379522151745,((J:2.0259528567288494,I:1.9642592877235094)i28:1.572419345579662,(L:1.4249957286369859,(P:1.1098894234123262,O)i31:0.0)i30:0.387910133713086)i29:2.114256390507651)i25:1.5444794449753867)i24:1.1415281321024846)i23:2.4399772182326847)i22:1.760845106315608))i21:1.9484562993805583;";
        String gtString2 = "(Z:6.5032896712070904,((A:6.5032896712070904,M:6.5032896712070904)i20:4.711530201978998,(((B:1.21493846897087,N:1.2015614615369377)i15:0.5266215533418312,K:1.6769779365754285)i16:6.5032896712070904,(((H:2.4173133586554716,G:2.4515047234037515)i10:1.8649265178632008,(F:2.632088660299164,E:3.2451931331855737)i7:1.2856402077264741)i11:4.502775388116889,(C:6.5032896712070904,(D:4.557379522151745,((J:2.0259528567288494,I:1.9642592877235094)i28:1.572419345579662,(L:1.4249957286369859,(P:1.1098894234123262,O)i31:0.0)i30:0.387910133713086)i29:2.114256390507651)i25:1.5444794449753867)i24:1.1415281321024846)i23:2.4399772182326847)i22:1.760845106315608))i21:1.9484562993805583;";
        List<Tree> gts = new ArrayList<>();
        STITree gtree = treeRoot.removeOutgroup(gtString1, "Z");
        STITree gtree2 = treeRoot.removeOutgroup(gtString2, "Z");
        gts.add(gtree);
        gts.add(gtree2);
        gts.add(gtree);
        gts.add(gtree2);
        HashMap<Tree, Integer> diffGts = getDifferentGTs(gts);
    }

    public static void checkAllDiff(){
//        String treefolderpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/proj/code/data1000gt/";
        String treefolderpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/IQTreeAndAstralTree/";
        ArrayList<String> treefilepaths = treeRoot.getFolers(treefolderpath);
        Collections.sort(treefilepaths);
        for(String treefilename : treefilepaths){
            String[] arr = treefilename.split("/");
            String treeName = arr[arr.length - 1];
//            checkTreeDifference(treefolderpath, treeName, 0);
            checkTreeDifference(treefolderpath, treeName, -1);
        }
    }
    public static void compareRAxMLIQTree(){
        ArrayList<String> RightRAxML = new ArrayList<String>();
        RightRAxML = runRAxML();
        ArrayList<String> RightIQ = new ArrayList<>();
        RightIQ = concatedIQST();
        Collections.sort(RightRAxML);
        Collections.sort(RightIQ);
        for(String s: RightRAxML){
            System.out.print(s+",");
        }
        System.out.println();
        for(String s: RightIQ){
            System.out.print(s+",");
        }
    }

    public static void compareSkinks(){
        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rainbow_skink_data/iqtree";
        List<String> files = treeRoot.getFiles(path);
        List<Tree> gts = new ArrayList<>();
        for(String file: files){
            try {
                String treestring = treeRoot.FileToString(file);
                gts.add(new STITree<>(treestring));
            }catch (Exception e){
                System.err.println("[error compareSkinks] "+e);
            }
        }
        System.out.println(gts.size());
        HashMap<Tree, Integer> diffGts = getDifferentGTs(gts);
    }

    public static void main(String[] args) {
//        ArrayList<String> RightRAxML = new ArrayList<String>();
//        ArrayList<String> RightASTRAL = new ArrayList<String>();
//        RightRAxML=runRAxML();
//        System.out.println("*******************************************************");
//        RightASTRAL=runASTRAL();
//        System.out.println(RightRAxML);
//        System.out.println("*******************************************************");
//
//        System.out.println(RightASTRAL);

//        run("/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rootedSpeciesTree_RAxML");
//        IQASTRAL();
//        checkTreeDifference("/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/proj/code/data1000gt/", "Reti1_A", 0);
//        TestgetDifferentGTs();
//        checkAllDiff();
        compareSkinks();
    }
}
