package edu.rice.cs.bioinfo.programs.phylonet.algos.treeAugment;
/*
 * @author: Zhen Cao
 * @Date: 1/19/2019
 * This class is for utility
 * */
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class treeRoot {
    public treeRoot(){

    }

    public static String FileToString(String Filename)
    {
        String str="";
//        File file=new File(Filename);
        try {
            FileInputStream in=new FileInputStream(new File(Filename));
            int size=in.available();
            byte[] buffer=new byte[size];
            in.read(buffer);
            in.close();
            str=new String(buffer);

        } catch (IOException e) {
            // TODO Auto-generated catch block
            StringWriter sw = new StringWriter();
            PrintWriter pw = new PrintWriter(sw);
            e.printStackTrace(pw);
        }
        return str;
    }

    public static void StringToFile(String s, String Filename)
    {

        try {
            File writename = new File(Filename);
            writename.createNewFile();
            BufferedWriter out = new BufferedWriter(new FileWriter(writename));
            out.write(s);
            out.flush();
            out.close();


        } catch (IOException e) {
            // TODO Auto-generated catch block
            StringWriter sw = new StringWriter();
            PrintWriter pw = new PrintWriter(sw);
            e.printStackTrace(pw);
        }
    }

    public static ArrayList<String> getFiles(String path) {
        ArrayList<String> files = new ArrayList<String>();
        File file = new File(path);
        File[] tempList = file.listFiles();

        for (int i = 0; i < tempList.length; i++) {
            if (tempList[i].isFile()) {
                files.add(tempList[i].toString());
            }

        }
        return files;
    }


    public static ArrayList<String> getFolers(String path) {
        ArrayList<String> files = new ArrayList<String>();
        File file = new File(path);
        File[] tempList = file.listFiles();

        for (int i = 0; i < tempList.length; i++) {
            if (tempList[i].isDirectory()) {
                files.add(tempList[i].toString());
            }
        }
        return files;
    }
    public static STITree removeOutgroup(String tree, String nodename){
//        System.out.println(tree);
        STITree gtree  = null;
        try{
            gtree = new STITree(tree);
        }
        catch (Exception e){
            System.err.println("removeOutgroup Tree: "+e);
        }
        STINode outNode = gtree.getNode(nodename);
        if(outNode == null){
            return gtree;
        }
        STINode parent = outNode.getParent();
        for(Object o : parent.getChildren()) {

            STINode sibling = (STINode) o;
            if(!sibling.equals(outNode)){
                gtree.rerootTreeAtEdge(sibling);
                gtree.removeNode(nodename);
                Trees.removeBinaryNodes(gtree);
                break;
//                Trees.rootAndRemoveOutgroup(gtree, nodename);
            }
        }

        return gtree;
    }


    public static Network removeOutgroup(Network net, String outNodeName){
        NetNode outNode = net.findNode(outNodeName);
        for(Object o : outNode.getParents()) {
            NetNode parent = (NetNode) o;
            parent.removeChild(outNode);
//            for(Object oc: parent.getChildren()) {
//                NetNode sibling = (NetNode) oc;
//                sibling.setParentDistance(parent, NetNode.NO_DISTANCE);
//            }
        }
        Networks.removeBinaryNodes(net);
        return net;
    }

    public static STITree rootTree(String goodTree, String outgroup){
        Network network = Networks.readNetwork(goodTree);
        int index = 0;
        for (Object o : Networks.postTraversal(network)) {
            NetNode node = (NetNode) o;
            if (!node.isLeaf()) {
                node.setName("i" + index);
            }

            index++;
        }
        STITree stitree = null;
        String gt = network.toString();
//        STITree gtree  = null;
        try{
            stitree = new STITree(gt);
            stitree.rerootTreeAtEdge(outgroup);
        }
        catch (Exception e){
            System.err.println("rootTree():"+e);
        }

        return stitree;
    }

    public static String rootTree(String goodTree, String outgroup, String deleteindividual) throws Exception{
        Network network = Networks.readNetwork(goodTree);
        int index = 0;
        for (Object o : Networks.postTraversal(network)) {
            NetNode node = (NetNode) o;
            if (!node.isLeaf()) {
                node.setName("i" + index);
            }

            index++;
        }
        String gt = network.toString();
        STITree gtree  = null;
        try{
            gtree = new STITree(gt);
        }
        catch (Exception e){
            System.err.println("rootTree():"+e);
        }

//        STITree gtree = new STITree(gt);
        STINode outgroupN = gtree.getNode(outgroup);
        gtree.rerootTreeAtEdge(outgroupN);
        if(gtree.getNode(deleteindividual) != null){
            gtree.removeNode(deleteindividual);
            network = Networks.readNetwork(gtree.toString());
            index = 0;
            for (Object o : gtree.postTraverse()) {
                STINode node = (STINode) o;
                if (!node.isLeaf()) {
                    node.setName("in" + index);
                    index++;
                }
            }

        }

        return gtree.toString();
    }

    public static String rootTrees(String [] goodTrees, String outgroup, String deleteindividual){

//        ArrayList<String> trees = new ArrayList<String>();
        String s="";
        for(String goodTree: goodTrees) {
//            System.out.println(goodTree);
             try {
                 s += rootTree(goodTree, outgroup, deleteindividual);
                 s += "\r\n";

//                       trees.add(s);
             } catch (Exception e) {
                 System.err.println("rootTrees"+e);

             }
        }
//        System.out.println("**************");
//        System.out.println(s);

        return s;
    }
    public static void rootSpeciesTree(String folderpath, String outpath, String outgroup){

//        String outgroup = "Z";

        ArrayList<String> filenames = getFiles(folderpath);

        for(String filename : filenames) {
            String fileString = FileToString(filename);
            String[] paths = filename.split("/");
            String outfilename = outpath + paths[paths.length - 1];
//            System.out.println(outfilename);
            try {
                StringToFile(rootTree(fileString, outgroup, "nodelete"), outfilename);
//                System.out.println(rootTree(fileString, outgroup));
            } catch (Exception e) {
                System.err.println(e);
            }
        }
    }

    public void rootGeneTrees(){

        String folderpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rainbow_skink_data";
        String outpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rainbow_skink_data/rgts.tre/";
        String outgroup = "Z_0";

        ArrayList<String> filenames = getFiles(folderpath);

        for(String filename : filenames){


            String fileString = FileToString(filename);
            String[] trees=fileString.split("\n");
            String[] paths = filename.split("/");
            if(paths[paths.length-1].startsWith("R")){
                String outfilename = outpath+paths[paths.length-1];
                try{
//                    STITree rootedgt = rootTrees(trees, outgroup);
                    StringToFile(rootTrees(trees, outgroup, ""), outfilename);
//                    System.out.println(rootTrees(trees, outgroup));
                } catch (Exception e){
                    System.err.println(e);
                }
            }

        }
    }

    public void rootGeneTreesAndDelete(){

        String folderpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/bestgTOneFILE";
        String outpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rootedgT/";
        String outgroup = "Z_0";

        ArrayList<String> filenames = getFiles(folderpath);

        for(String filename : filenames){


            String fileString = FileToString(filename);
            String[] trees=fileString.split("\n");
            String[] paths = filename.split("/");
            if(paths[paths.length-1].startsWith("R")){
                String outfilename = outpath+paths[paths.length-1];
                try{
//                    STITree rootedgt = rootTrees(trees, outgroup);
//                    STITree
                    StringToFile(rootTrees(trees, outgroup, ""), outfilename);
//                    System.out.println(rootTrees(trees, outgroup));
                } catch (Exception e){
                    System.err.println(e);
                }
            }

        }
    }
//    public void rootIQtree(){
//        String folderpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/IQTreeAndAndAstralTree/";
////        String outpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rootedSpeciesTree_Astral/";
//
//        ArrayList<String> foldernames = getFiles(folderpath);
//
//        for(String fname : foldernames) {
//
//
//        }
//
//        rootSpeciesTree(folderpath, outpath);
//
//
//
//
//
//    }


    public void writemap(){
        String gtpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rainbow_skink_data/rgts.tre";
        String gtstring = treeRoot.FileToString(gtpath);
        STITree gt = null;
        HashMap<String, List<String>> taxonmap = new HashMap<>();
        try {
            gt = new STITree(gtstring);
            String [] leaves = gt.getLeaves();
            for(String leaf:leaves){
                String [] arr = leaf.split("_");
                String taxon = arr[0]+arr[1];
                List<String> individuals = new ArrayList<>();
                if(taxonmap.containsKey(taxon)){
                    individuals = taxonmap.get(taxon);
                    individuals.add(leaf);
                    taxonmap.remove(taxon);
                    taxonmap.put(taxon, individuals);
                }
                else{
                    individuals.add(leaf);
                    taxonmap.put(taxon, individuals);
                }
            }
        }catch (Exception e){
            System.out.println("writemap wrong");
        }
        String result = "";
        for(Map.Entry<String, List<String>> entry: taxonmap.entrySet()){
            result += entry.getKey();
            result += ":";
            result += String.join(",", entry.getValue());
            result += "\n";
        }
        StringToFile(result, "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rainbow_skink_data/map.txt");
        System.out.println(result);
    }

    public static void rerootAstral(){

        treeRoot t = new treeRoot();
        String stpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rainbow_skink_data/lizard_st.tre";
        String outgroup = "SP07indexing28";
        String ststring = FileToString(stpath);
        STITree stree = treeRoot.rootTree(ststring, outgroup);
        System.out.println(stree.toString());
        StringToFile(stree.toString(), "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rainbow_skink_data/rst.tre");
    }
    public static void mapleaves(){
        HashMap<String, String> map = new HashMap<>();
//        String net = "((((SP07indexing11:1.0,(((SP07indexing30:1.0,SP09indexing15:1.0):0.17039775050774653,SP07indexing29:1.0):2.442667272328067)#H1:0.0011774181844964955::0.5453077749386916):0.901504939784219,(#H1:5.938859691040344::0.45469222506130835,((SP09indexing3:1.0,SP03indexing25:1.0):0.4088109758343756,(SP10indexing20:1.0,SP07indexing9:1.0):0.0858406183246179):1.10603492051881):0.7466203937102353):2.7346197411319366,SP09indexing21:1.0):1.8541019662496847,SP02Aindexing4:1.0);";
//        String net = "(SP07indexing28:4.209655408733094,(SP02Aindexing4:3.1110431200649846,(SP09indexing21:4.209655408733094,((((SP07indexing30:4.209655408733094,SP09indexing15:2.000160738805061)i9:0.18765755009637727,SP07indexing29:4.209655408733094)i10:2.1302138670532598,SP07indexing11:4.209655408733094)i11:0.2122195254568186,((SP09indexing3:4.209655408733094,SP03indexing25:3.573666642013102)i4:0.4352537811809968,(SP10indexing20:3.6342912638295357,SP07indexing9)i19:0.0)i18:0.07391758662244717)i17:1.5448897224076301)i16:2.948003492820481))i14:2.062263205144751;";
        String net = "(SP07indexing28:0.0238584516,(SP02Aindexing4:0.018927614199999998,(SP09indexing21:0.0227495387,((SP07indexing11:0.014226511099999999,(SP07indexing30:0.0103455713,(SP07indexing29:0.0101497146,SP09indexing15:0.0069544131)i8:9.740611E-4)i10:0.0065751779)i12:0.0015760351,(SP10indexing20:0.009432411799999998,(SP07indexing9:0.014994415,(SP09indexing3:0.0124550144,SP03indexing25)i19:0.0167588245)i17:0.0015003071)i15:6.431855E-4)i13:0.00511702)i5:0.0100453832))i3:0.0089631301;";
        map.put("SP07indexing28", "Lampropholis guichenoti");
        map.put("SP02Aindexing4", "Lampropholis coggeri");
        map.put("SP09indexing21", "Pygmaeascincus timlowi");
        map.put("SP03indexing25", "Carlia amax");
        map.put("SP09indexing3", "Carlia vivax");
        map.put("SP07indexing9", "Carlia longipes");
        map.put("SP10indexing20", "Carlia rhomboidalis");
        map.put("SP07indexing11", "Liburnascincus mundivensis");
        map.put("SP07indexing29", "Lygisaurus foliorum");
        map.put("SP07indexing30", "Lygisaurus cf. macfarlani");
        map.put("SP09indexing15", "Lygisaurus sesbrauna");
        Network n = Networks.readNetwork(net);
        for(Object o: n.getLeaves()){
            NetNode node = (NetNode) o;
            if(map.containsKey(node.getName())){
                node.setName(map.get(node.getName()));
            }
        }
        System.out.println(n.toString());
    }

    static public void main(String[] args){
        mapleaves();
//        treeRoot t = new treeRoot();
//        String gtpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rainbow_skink_data/gts.tre";
//        String[] gts = FileToString(gtpath).split("\n");
//        String result = rootTrees(gts, "SP07_indexing28_h0", "SP07_indexing28_h1");
//        StringToFile(result, "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rainbow_skink_data/rgts.tre");
//        Network net = Networks.readNetwork("((((((((L:1.728,(M:1.44)#H1:0.28800000000000003::0.8000000000000002)S8:20.458111067404356,(P:1.0,O:1.0)S7:21.186111067404358)S6:4.437222213480872,((C:8.916100448255998,B:8.916100448255998)S10:9.572325441247633,(E:6.191736422399999)#H2:12.296689467103633::0.2499999999999999)S9:8.1349073913816)S5:19.381786628484445,((G:4.299816959999999,(N:1.2)#H3:3.099816959999999::0.65)S12:8.539367685488635,(((K:2.0736,#H1:0.6335999999999999::0.19999999999999984)S15:0.41472,J:2.48832)S14:2.6714603519999995,F:5.159780351999999)S13:7.679404293488635)S11:33.16593526388104)S4:9.20102398187393,A:55.206143891243606)S3:11.041228778248716,(#H2:25.756263514662276::0.7500000000000001)#H4:34.299372732430044::0.7000000000000001)S2:13.249474533898464,(((((#H3:2.3831808::0.35)#H5:3.846902906879999::0.4,D:7.430083706879999)S20:3.2692368310271975,(H:2.9859839999999997,I:2.9859839999999997)S19:7.713336537907196)S18:4.707701036679165,#H5:11.823840774586362::0.6)S17:22.93057834988837,#H4:6.389599987412456::0.29999999999999993)S16:41.159247278916055)S1:20.503152796609214,Z:100.0);");
////        String tree = "((((((((L:1.728,(M:1.44)#H1:0.28800000000000003::0.8000000000000002)S8:20.458111067404356,(P:1.0,O:1.0)S7:21.186111067404358)S6:4.437222213480872,((C:8.916100448255998,B:8.916100448255998)S10:9.572325441247633,(E:6.191736422399999)#H2:12.296689467103633::0.2499999999999999)S9:8.1349073913816)S5:19.381786628484445,((G:4.299816959999999,(N:1.2)#H3:3.099816959999999::0.65)S12:8.539367685488635,(((K:2.0736,#H1:0.6335999999999999::0.19999999999999984)S15:0.41472,J:2.48832)S14:2.6714603519999995,F:5.159780351999999)S13:7.679404293488635)S11:33.16593526388104)S4:9.20102398187393,A:55.206143891243606)S3:11.041228778248716,(#H2:25.756263514662276::0.7500000000000001)#H4:34.299372732430044::0.7000000000000001)S2:13.249474533898464,(((((#H3:2.3831808::0.35)#H5:3.846902906879999::0.4,D:7.430083706879999)S20:3.2692368310271975,(H:2.9859839999999997,I:2.9859839999999997)S19:7.713336537907196)S18:4.707701036679165,#H5:11.823840774586362::0.6)S17:22.93057834988837,#H4:6.389599987412456::0.29999999999999993)S16:41.159247278916055)S1:20.503152796609214,Z:100.0);";
//        String tree = "(Z:4.209655408733094,((((F:2.586972269548974,E:2.8233610476132034)i24:3.436465520499609,((N:0.9709769565687146,M:1.2918846766488157)i20:0.4660510549299122,J:1.7247487589450947)i21:3.6988297849671126)i25:1.533866069849107,(((L:0.8620932429035711,K:1.5309195952996604)i14:1.1186129553747794,G:2.4511819171385074)i15:2.3854421946553037,(H:2.1302138670532598,I:2.121915064238565)i10:3.1110431200649846)i16:2.788269727801935)i26:1.0569193863694388,(A:4.209655408733094,((D:2.1302138670532598,C:3.1110431200649846)i5:2.743318339939666,(B:3.4814169083618802,(P:0.5853144757567299,O)i31:0.9325106757409185)i30:3.0992086624179858)i29:1.49606502725527)i28:3.0148150879704287))i27:0.668696084695781;";
//        STITree stitree = treeRoot.removeOutgroup(tree, "Z");
//        System.out.println(stitree.getRoot().getName());
//        System.out.println( stitree.toString());
//
//        treeRoot.removeOutgroup(net, "Z");
//        System.out.println(net.toString());
//        t.writemap();
        //        String treeout = t.removeOutgroup(net, "Z");
//        t.removeOutgroup(net, "Z");

//        System.out.println(net.toString());
//        t.rootGeneTrees();


//        String folderpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/speciesTree_Astral";
//        String outpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rootedSpeciesTree_Astral/";
//        rootSpeciesTree(folderpath, outpath);
//
//        t.rootGeneTrees();



//        String STRAxML = "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/speciesTree_";
//        String rSTRAxML= "/Users/zhen/Desktop/Zhen/research/phylogenetics/scalablePhylogeneticNetworkInference/data/rootedSpeciesTree/rooted";
//        t.rootSpeciesTree(STRAxML, rSTRAxML);

        //        String s="((((L_0:0.00501707538645165818,L_1:0.00000100000050002909):0.00941055181733906759,(P_1:0.00000100000050002909,P_0:0.00100067719709055021):0.01182005505507253973):0.01276548200760693579,(((((((B_0:0.00737530323592503270,B_1:0.01181953322142060767):0.05390659093152030262,(A_1:0.01776624308817553732,A_0:0.02312569256254220357):0.07156671895159436025):0.06537911151269693022,((K_1:0.00585318703355132246,K_0:0.00319042244355325104):0.06882767864767656840,(((M_1:0.00201515161775757011,M_0:0.00199306551570837404):0.01541113621641108933,(N_1:0.00198649305466938814,N_0:0.00202170815849337607):0.01719286740184226836):0.04770703333734686857,(Z_0:0.00000100000050002909,Z_1:0.00401013805913900522):1.19380895112187146445):0.00574304177256377071):0.05870695609395511910):0.03000135692794294456,(((G_0:0.01071291759534525635,G_1:0.01175948736434921302):0.02007136690414206698,(H_1:0.00000100000050002909,H_0:0.00000100000050002909):0.02754612191114961042):0.00840566063506667689,((E_1:0.02410342362752090559,E_0:0.02264813709010397164):0.00484308814207149579,(F_1:0.02510544034454245677,F_0:0.02097506868739815872):0.01803472997600509900):0.00406669589923412417):0.05031690887827611897):0.03943585870205408572,(C_0:0.00199958075264799638,C_1:0.00200873250677360360):0.07292966879138924885):0.00550789756414620142,(D_0:0.00607037189001380546,D_1:0.00398272908021218754):0.06194193000479023686):0.01374969581570348147,((I_0:0.00000100000050002909,I_1:0.00000100000050002909):0.03001541130931727891,(J_1:0.01214772436426610723,J_0:0.01317949661937398501):0.01936250440772952289):0.01115564845063425625):0.01444205435482900229):0.01400102317616337168,O_0:0.00199586200552785871,O_1:0.00402337370525488597):0.0;";
//        String outcome="";
//        try{
//            outcome=t.rootTree(s, "Z");
////            System.out.println("root"+rootTree(s, "Z"));
//        }catch (Exception e){
//            System.err.println(e);
//        }
//
//        System.out.println(outcome);

    }
}
