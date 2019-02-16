package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import com.google.gson.Gson;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.util.CombinatoricsUtils;

import java.io.*;
import java.nio.file.Paths;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 10/19/18
 * Time: 10:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimTest {
    public static Network PrepareSimulatedNetwork(String input) {
        String netstring = Pipeline.convertBeastNetworkString(input);
        Network<NetNodeInfo> net = Networks.readNetwork(netstring);

        NetNode rootChild = net.getRoot().getChildren().iterator().next();
        NetNode root = net.getRoot();
        // Reset height
        rootChild.setParentDistance(root, rootChild.getParentDistance(root) + 0.25);

        // Scale branch lengths
        for(NetNode<NetNodeInfo> node : net.dfs()) {
            for(NetNode parent : node.getParents()) {
                node.setParentDistance(parent, node.getParentDistance(parent) * 200.0);
            }
        }

        // Reset inheritance probabilities
        double prob = 0.6;
        for(NetNode<NetNodeInfo> node : net.dfs()) {
            if(node.isNetworkNode()) {
                Iterator<NetNode<NetNodeInfo>> it = node.getParents().iterator();
                NetNode parent1 = it.next();
                NetNode parent2 = it.next();
                node.setParentProbability(parent1, prob);
                node.setParentProbability(parent2, 1.0 - prob);
                prob += 0.05;
            }
        }

        Pipeline.initNetHeights(net);

        double height = net.getRoot().getData().getHeight();

        System.out.println("Height: " + height);

        BniNetNode<NetNodeInfo> outgroup = new BniNetNode<>();
        outgroup.setData(new NetNodeInfo(0.0));
        outgroup.setName("Z");
        net.getRoot().adoptChild(outgroup, height);

        // Reset heights
        List<NetNode<NetNodeInfo>> internalNodes = new ArrayList<>();
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(net)) {
            if(!node.isLeaf()) {
                internalNodes.add(node);
            }
        }
        Collections.sort(internalNodes, (NetNode<NetNodeInfo> a, NetNode<NetNodeInfo> b)->Double.compare(a.getData().getHeight(), b.getData().getHeight()));

        double curheight = 1.0;
        double increment = 1.0;
        for(NetNode<NetNodeInfo> node : internalNodes) {
            if(!node.isLeaf()) {
                node.getData().setHeight(curheight);
                curheight *= 1.2;//+= increment;
                increment += 0.2;
            }
        }

        net.getRoot().getData().setHeight(100.0);

        // Reset branch lengths
        for(NetNode<NetNodeInfo> node : net.dfs()) {
            for(NetNode<NetNodeInfo> parent : node.getParents()) {
                node.setParentDistance(parent, parent.getData().getHeight() - node.getData().getHeight());
            }
        }

        System.out.println("Full string: " + Networks.getFullString(net));
        System.out.println("Dendroscope: " + Networks.getDendroscopeCompatibleString(net));

        return net;
    }

    static void WalkThroughNetworks() {
        String filename = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/no_parallel_edges.trees";
        int count = 0;
        int exceptions = 0;
        int readfailed = 0;
        int d3net = 0;
        List<String> problematic = new ArrayList<>();
        Map<Integer, List<Network>> retiNumToNet = new TreeMap<>();
        Map<Network, Double> scores = new HashMap<>();

        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String s;
            int index = 0;
            while((s = in.readLine()) != null) {
                //if(index < 7182) {index++;continue;}
                //if(index == 682) {index++;continue;}

                try {
                    Network trueNetwork = PrepareSimulatedNetwork(s);
                    if(trueNetwork == null) {index++;readfailed++;continue;}


                    int retiNum = trueNetwork.getReticulationCount();
                    if(!retiNumToNet.containsKey(retiNum)) {
                        retiNumToNet.put(retiNum, new ArrayList<>());
                    }

                    retiNumToNet.get(retiNum).add(trueNetwork);

                    scores.put(trueNetwork, Pipeline.ComputeNetworkDifficulty(trueNetwork));
                } catch (Exception e) {
                    //e.printStackTrace();
                    problematic.add(s);
                    exceptions++;
                    //break;
                }
                index++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        Map<String, String> keyToNetwork = new TreeMap<>();

        for(Integer retiNum : retiNumToNet.keySet()) {
            System.out.println(retiNum + "\t" + retiNumToNet.get(retiNum).size());

            Collections.sort(retiNumToNet.get(retiNum), (Network a, Network b)->Double.compare(scores.get(a), scores.get(b)));
            List<Double> s = new ArrayList<>();
            for(Network net : retiNumToNet.get(retiNum)) {
                s.add(scores.get(net));
            }

            keyToNetwork.put("Reti" + retiNum + "_A", retiNumToNet.get(retiNum).get(0).toString());
            keyToNetwork.put("Reti" + retiNum + "_B", retiNumToNet.get(retiNum).get((int)(retiNumToNet.get(retiNum).size() / 3.0)).toString());
            keyToNetwork.put("Reti" + retiNum + "_C", retiNumToNet.get(retiNum).get((int)(retiNumToNet.get(retiNum).size() * 2.0 / 3.0)).toString());
            keyToNetwork.put("Reti" + retiNum + "_D", retiNumToNet.get(retiNum).get(retiNumToNet.get(retiNum).size() - 1).toString());
        }

        for(String key : keyToNetwork.keySet()) {
            System.out.println(key);
            System.out.println(keyToNetwork.get(key));
        }

        String outfilename = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/networks.json";
        try {
            PrintWriter out = new PrintWriter(outfilename);
            out.println("{");
            out.println("\"networks\":[");
            boolean first = true;
            for(String key : keyToNetwork.keySet()) {
                if(first) first = false;
                else out.println(",");
                out.println("\t{");
                out.println("\t\"tag\": \"" + key + "\",");
                out.println("\t\"netstring\": \"" + keyToNetwork.get(key) + "\"");
                out.print("\t}");
            }
            out.println("]");
            out.println("}");
            out.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static Map<String, String> readIQTree(String resultFolder) {
        File path = new File(resultFolder);

        List<String> filenames = new ArrayList<>();
        Map<String, String> filename2locusname = new HashMap<>();
        Map<String, String> locus2treestring = new TreeMap<>();

        File[] files = path.listFiles();

        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){
                if(files[i].getName().endsWith(".treefile")) {
                    String ss[] = files[i].getName().split("\\.");
                    String sss[] = ss[0].split("_");
                    filename2locusname.put(files[i].toString(), sss[1]);
                    filenames.add(files[i].toString());
                }
            }
        }

        Collections.sort(filenames);
        System.out.println("IQTREE read: " + filenames.size());
        for(String filename : filenames) {
            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                s = in.readLine();
                STITree tree = null;
                try{ tree = new STITree(s); } catch (Exception e) {e.printStackTrace();}

                for(Object nodeObj : tree.postTraverse()) {
                    STINode node = (STINode) nodeObj;
                    if(node.isLeaf()) {
                        String name = node.getName();
                        // convert name of leaves from iqtree
                        name = name.charAt(0) + "_" + name.substring(1);
                        node.setName(name);
                    }
                }

                locus2treestring.put(filename2locusname.get(filename), tree.toNewick());

            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        return locus2treestring;
    }

    static void prepare(String[] args) {
        int nchar = 1000000;
        String root = args.length > 0 ? args[0] : ".";
        String outgroup = null;
        String gtoutgroup = "Z_1";
        boolean sgt = true;


        Map<String, String> allele2species = new HashMap<>();
        Map<String, List<String>> species2alleles = new HashMap<>();
        List<String> species_list = new ArrayList<>();
        for(int i = 0 ; i < 16 ; i++) {
            char ch = (char)((int)'A' + i);
            String species = "";
            species += ch;
            species_list.add(species);
            species2alleles.put(species, new ArrayList<>());
        }
        species_list.add("Z");
        species2alleles.put("Z", new ArrayList<>());

        for(String species : species_list) {
            for(int j = 0 ; j < 2 ; j++) {
                String allele = "";
                allele += species;
                allele += "_";
                allele += j;
                allele2species.put(allele, species);
                species2alleles.get(species).add(allele);
            }
        }


        int N = species2alleles.size();
        int total = N;
        if (outgroup != null) total--;

        total = (int) CombinatoricsUtils.binomialCoefficient(total, 3);

        File path = new File(root);

        List<String> dirnames = new ArrayList<>();

        File [] files = path.listFiles();
        for (int i = 0; i < files.length; i++){
            if (files[i].isDirectory()){
                if(files[i].getName().startsWith("Reti")) {
                    System.out.println("Detected: " + files[i]);
                    dirnames.add(files[i].getName());
                }
            }
        }

        Collections.sort(dirnames);

        for(String dirname : dirnames) {
            Map<String, String> locus2tree = null;
            if(sgt) {
                locus2tree = readIQTree(Paths.get(root, "iqtree", dirname).toString());
            }

            for(int nnn = 0 ; nnn < total ; nnn++) {
                String nexusPath = Paths.get(root, dirname,  "run_" + nnn + ".nex").toString();
                System.out.println(nexusPath);

                try {
                    PrintWriter out = new PrintWriter(nexusPath);

                    out.print(String.format("#NEXUS \nBegin data;\nDimensions ntax=%d nchar=%d;\nFormat datatype=dna symbols=\"012\" missing=? gap=-;\nMatrix\n\n", species2alleles.size(), nchar));

                    BufferedReader in = new BufferedReader(new FileReader(Paths.get(root, "markers", dirname, "markers_0.txt").toFile()));
                    String line = null;
                    while ((line = in.readLine()) != null) {
                        out.println(line);
                    }
                    in.close();

                    out.println();
                    out.println(";End;\n");

                    if(sgt) {
                        out.println("BEGIN TREES;");
                        for(String locus : locus2tree.keySet()) {
                            out.println("Tree " + locus + " = " + locus2tree.get(locus));
                        }
                        out.println("End;");
                        out.println();
                    }

                    out.println("BEGIN PHYLONET; \n");
                    out.println("SN_SEQ -cl 6000000 -bl 3000000 -sf 5000 -ee -sd 123456 -pl 2 -gtburnin -pre 10 ");
                    if(outgroup != null) {
                        out.println("-outgroup \"" + outgroup +"\"");
                    }

                    if(gtoutgroup != null) {
                        out.println("-gtoutgroup \"" + gtoutgroup +"\"");
                    }

                    if(sgt) {
                        out.print("-sgt (");
                        boolean first = true;
                        for(String locus : locus2tree.keySet()) {
                            if(first) first = false;
                            else out.print(",");
                            out.print("" + locus);
                        }
                        out.println(") ");
                    }

                    out.print("-tm <");
                    int index = 0;
                    for (String species : species2alleles.keySet()) {
                        if (index > 0) out.print(";");
                        out.print(species + ":");
                        for (int i = 0; i < species2alleles.get(species).size(); i++) {
                            if (i > 0) out.print(",");
                            out.print(species2alleles.get(species).get(i));
                        }
                        index++;
                    }

                    out.println("> ");
                    out.println(" -nnn " + nnn + " ;");
                    out.println("END;");

                    out.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

    }

//    static void CheckSimulationResult(String tag, String beastString, String resultFolder) {
//        Network trueNetwork = PrepareSimulatedNetwork(beastString);
//        File path = new File(resultFolder);
//
//        List<String> filenames = new ArrayList<>();
//
//        File [] files = path.listFiles();
//        for (int i = 0; i < files.length; i++){
//            if (files[i].isFile()){ //this line weeds out other directories/folders
//                if(files[i].toString().endsWith(".out")) {
//                    System.out.println(files[i]);
//                    filenames.add(files[i].toString());
//                }
//            }
//        }
//
//        Collections.sort(filenames);
//        SuperNetwork3 sn = Pipeline.stage2(filenames, 400000, 80000, 1000);
//        Network inferred = sn.compute();
//
//        if(Networks.hasTheSameTopology(trueNetwork, inferred)) {
//            System.out.println(tag + " : " + "Good");
//        } else {
//            System.out.println(tag + " : " + "Not Good");
//        }
//    }

    static class ResultInfo {
        Tuple3<Network, Network, Double> closest;
    }

    static void CheckIdealResult(String tag, String netString, String reduceFolder) {
        Map<String, String> allele2species = new HashMap<>();
        List<String> species_list = new ArrayList<>();
        for(int i = 0 ; i < 16 ; i++) {
            char ch = (char)((int)'A' + i);
            String species = "";
            species += ch;
            species_list.add(species);
        }
        species_list.add("Z");

        for(String species : species_list) {
            for(int j = 0 ; j < 2 ; j++) {
                String allele = "";
                allele += species;
                allele += j;
                allele2species.put(allele, species);
            }
        }



        System.out.println("Tag: " + tag);

        Network trueNetwork = Networks.readNetwork(netString);

        SNOptions options = new SNOptions();
        if(reduceFolder != null)
            options.tripletFilename = reduceFolder + "/" + tag + "/triplets.txt";
        options.trueNetwork = trueNetwork.clone();
        options.reconcileHeights = true;
        options.allele2species = allele2species;
        SuperNetwork3.printDetails_ = true;

        SNSummary summary = Pipeline.stage2_1(null, 0, 0, 1, options);
        Network inferred = summary.inferredNetwork;

        System.out.println("True: " + Networks.getDendroscopeCompatibleString(trueNetwork));
        System.out.println("Inferred: " + inferred);
        System.out.println("Inferred # Reti: " + inferred.getReticulationCount());

        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(inferred, trueNetwork);
        System.out.println("Closest true # Reti: " + closest.Item1.getReticulationCount());
        System.out.println(closest.Item2);
        System.out.println("Closest inferred # Reti: " + closest.Item2.getReticulationCount());
        System.out.println(closest.Item1);
        System.out.println("Distance: " + closest.Item3);
    }

    static ResultInfo CheckSimulationResult2(String tag, String netString, String resultFolder, String reduceFolder) {
        Map<String, String> allele2species = new HashMap<>();
        List<String> species_list = new ArrayList<>();
        for(int i = 0 ; i < 16 ; i++) {
            char ch = (char)((int)'A' + i);
            String species = "";
            species += ch;
            species_list.add(species);
        }
        species_list.add("Z");

        for(String species : species_list) {
            for(int j = 0 ; j < 2 ; j++) {
                String allele = "";
                allele += species;
                allele += j;
                allele2species.put(allele, species);
            }
        }



        System.out.println("Tag: " + tag);

        Network trueNetwork = Networks.readNetwork(netString);
        File path = new File(resultFolder);

        if(!path.exists()) return null;

        List<String> filenames = new ArrayList<>();

        File [] files = path.listFiles();
        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){ //this line weeds out other directories/folders
                if(files[i].toString().endsWith(".out")) {
                    //System.out.println(files[i]);
                    filenames.add(files[i].toString());
                }
            }
        }

        Collections.sort(filenames);
        SNOptions options = new SNOptions();
        if(reduceFolder != null)
            options.tripletFilename = reduceFolder + "/" + tag + "/triplets.txt";
        options.trueNetwork = trueNetwork.clone();
        options.reconcileHeights = true;
        options.allele2species = allele2species;
        SuperNetwork3.printDetails_ = true;
        //SNSummary summary = Pipeline.stage2_1(filenames, 6000000, 2000000, 5000, options); // good
        SNSummary summary = Pipeline.stage2_1(filenames, 2000000, 1000000, 5000, options);


        int trinetCorrect = 0;
        int trinetCorrectBackbone = 0;
        Mean error = new Mean();
        for(SuperNetwork3.NetworkWithInfo netinfo : summary.netinfos) {
            if(netinfo.dirty) {
                netinfo.network = Networks.readNetwork(netinfo.backup);
            }
            Tuple<Network, Map<NetNode, NetNode>> tuple = SuperNetwork3.getSubNetwork(trueNetwork, netinfo.taxa, true);
            Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueBackbone(netinfo.network, tuple.Item1);
            if(closest.Item3 == 0.0) {
                if(tuple.Item1.getReticulationCount() == closest.Item2.getReticulationCount()) {
                    trinetCorrect++;
                } else {
                    trinetCorrectBackbone++;
                }

                SuperNetwork3.initNetHeights(tuple.Item1);
                Map<NetNode, NetNode> mapTwo = SuperNetwork3.mapTwoNetworks(netinfo.network, tuple.Item1);
                if(mapTwo != null) {
                    for(NetNode node1 : mapTwo.keySet()) {
                        if(node1.isLeaf()) continue;
                        NetNode node2 = mapTwo.get(node1);
                        double height1 = 100.0 * ((NetNodeInfo) node1.getData()).getHeight();
                        double height2 = ((NetNodeInfo) node2.getData()).getHeight();
                        error.increment(Math.abs(height1 - height2) / height2);
                    }
                }
            }
        }

        System.out.println("# Trinets: " + summary.netinfos.size());
        System.out.println("# Trinets correct: " + trinetCorrect);
        System.out.println("# Trinets correct backbone: " + trinetCorrectBackbone);
        System.out.println("Avg time error: " + error.getResult());

        Network inferred = summary.inferredNetwork;


        ResultInfo result = new ResultInfo();

        System.out.println("True: " + Networks.getDendroscopeCompatibleString(trueNetwork));
        System.out.println("Inferred: " + inferred);
        System.out.println("Inferred # Reti: " + inferred.getReticulationCount());

        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(inferred, trueNetwork);
        result.closest = closest;
        System.out.println("Closest true # Reti: " + closest.Item2.getReticulationCount());
        System.out.println(closest.Item2);
        System.out.println("Closest inferred # Reti: " + closest.Item1.getReticulationCount());
        System.out.println(closest.Item1);
        System.out.println("Distance: " + closest.Item3);


        Set<String> trueHybridizations = new HashSet<>(Networks.getTaxaNamesUnderReticulation(trueNetwork));
        Set<String> inferredHybridizations = new HashSet<>(Networks.getTaxaNamesUnderReticulation(inferred));
        Set<String> truePositive = new HashSet<>(trueHybridizations);
        truePositive.retainAll(inferredHybridizations);
        Set<String> falsePositive = new HashSet<>(inferredHybridizations);
        falsePositive.removeAll(trueHybridizations);
        Set<String> falseNegative = new HashSet<>(trueHybridizations);
        falseNegative.removeAll(inferredHybridizations);
        Set<String> trueNegative = new HashSet<>(summary.taxaNames);
        trueNegative.removeAll(truePositive);
        trueNegative.removeAll(falseNegative);
        trueNegative.removeAll(falsePositive);

        System.out.println("True positive: " + truePositive.size() + "/" + trueHybridizations.size());
        System.out.println("False positive: " + falsePositive.size() + "/" + (trueNegative.size() + falsePositive.size()));
        System.out.println("False negative: " + falseNegative.size() + "/" + trueHybridizations.size());


        System.out.println();

        return result;
    }

    public static class NetworkListJson {
        public static class NetworkJson {
            String tag;
            String netstring;
        }

        List<NetworkJson> networks = new ArrayList();

        public List<NetworkJson> getNetworks() {
            return networks;
        }

        public void setNetworks(List<NetworkJson> networks) {
            this.networks = networks;
        }
    }

    public static void check24(String[] args) {
        String networkfilename = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/networks.json";
        String resultRoot = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run61/"; // run61 run47
        String reduceRoot = null;//"/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/Reduce/";
        int requiredIndex = 19;
        if(args.length > 0) {
            networkfilename = args[0];
            resultRoot = args[1].equals("null") ? null : args[1];
            reduceRoot = args[2].equals("null") ? null : args[2];
            if(!args[3].equals("def"))
                SuperNetwork3.eps = Double.parseDouble(args[3]);
            requiredIndex = Integer.parseInt(args[4]);
        }

        NetworkListJson networkListjson = null;
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(networkfilename));

            Gson gson = new Gson();
            networkListjson = gson.fromJson(bufferedReader, NetworkListJson.class);

        } catch(Exception e) {
            e.printStackTrace();
        }

        int index = 0;
        for(NetworkListJson.NetworkJson networkJson : networkListjson.networks) {
            if(index++ < requiredIndex) {continue;}
//            Network net = Networks.readNetworkWithRootPop(networkJson.netstring);
//            for(Object nodeObj : net.dfs()) {
//                NetNode node = (NetNode) nodeObj;
//                if(!node.isLeaf()) node.setName("");
//                for(Object parentObj : node.getParents()) {
//                    NetNode parent = (NetNode) parentObj;
//                    node.setParentDistance(parent, node.getParentDistance(parent) * 0.01);
//                }
//            }
//            System.out.println(Networks.getDendroscopeCompatibleString(net));
            if(resultRoot != null) {
                String resultPath = resultRoot + "/" + networkJson.tag;
                CheckSimulationResult2(networkJson.tag, networkJson.netstring, resultPath, reduceRoot);
            } else {
                CheckIdealResult(networkJson.tag, networkJson.netstring, reduceRoot);
            }
            break;
        }
    }

    public static void measureCorrectness() {
        Network trueNetwork = Networks.readNetwork("(Z:100.0,(((((E:4.299816959999999)#H1:3.1302667468799994,(D:5.159780351999999,#H1:0.859963392)S16:2.2703033548799993)S15:3.2692368310271975,(F:3.5831807999999996,(H:2.9859839999999997,G:2.9859839999999997)S14:0.5971967999999999)S13:7.116139737907197)S12:2.139864107581438,A:12.839184645488634)S11:5.649241244014997,(((M:1.2,N:1.2)S4:0.24,L:1.44)S3:13.967021574586362,(((((O:1.0,P:1.0)S10:0.728,K:1.728)S9:0.3455999999999999,J:2.0736)S8:0.41472,I:2.48832)S7:6.427780448255998,(C:6.191736422399999,B:6.191736422399999)S6:2.7243640258559987)S5:6.490921126330363)S2:3.08140431491727)S1:81.51157411049637);");
        Network inferredNetwork = Networks.readNetwork("(Z:0.7255560612824605,(((L:0.013874224671003466,(N:0.012007666301733881,M:0.012007666301733881)I7:0.0018665583692695814)I4:0.1318268794169739,((I:0.024485711596326017,(J:0.019893413894892884,((O:0.009847745975532486,P:0.009847745975532486)I15:0.007031621727710621,K:0.016879367703243112)I14:0.0030140461916497697)I12:0.004592297701433128)I8:0.0607196709051126,(B:0.06070577476604321,C:0.06070577476604321)I9:0.024499607735395418)I5:0.06049572158653875)I2:0.026629114814542613,(((D:0.053658715959090976,E:0.053658715959090976)I11:0.04810311687417632,(F:0.033686787717168064,(G:0.029213972928706228,H:0.029213972928706228)I13:0.004472814788461852)I10:0.06807504511609927)I6:0.01798915790458915,A:0.11975099073785656)I3:0.05257922816466351)I1:0.5532258423799403)I0;");
        System.out.println(Networks.hasTheSameTopology(trueNetwork, inferredNetwork));
    }

//    public static void check15(String[] args) {
//        Map<String, String> beastStrings = new HashMap<>();
//        beastStrings.put("#14", "(((((((J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.017178654766260792)#H1[&gamma=0.9437783623704922,gamma_95%HPD={0.9437783623704922,0.9437783623704922},gamma_mean=0.9437783623704922,gamma_median=0.9437783623704922,gamma_range={0.9437783623704922,0.9437783623704922},height_95%HPD={0.017178654766260792,0.017178654766260792},height_mean=0.017178654766260792,height_median=0.017178654766260792,height_range={0.017178654766260792,0.017178654766260792}]:0.0017490974609733478,I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01892775222723414)S5[&height_95%HPD={0.01892775222723414,0.01892775222723414},height_mean=0.01892775222723414,height_median=0.01892775222723414,height_range={0.01892775222723414,0.01892775222723414}]:0.009397644206168954,(K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.014259823085634504,L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.014259823085634504)S6[&height_95%HPD={0.014259823085634504,0.014259823085634504},height_mean=0.014259823085634504,height_median=0.014259823085634504,height_range={0.014259823085634504,0.014259823085634504}]:0.01406557334776859)S4[&height_95%HPD={0.028325396433403094,0.028325396433403094},height_mean=0.028325396433403094,height_median=0.028325396433403094,height_range={0.028325396433403094,0.028325396433403094}]:0.04526369353431495,(F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02613005187787895,G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02613005187787895)S7[&height_95%HPD={0.02613005187787895,0.02613005187787895},height_mean=0.02613005187787895,height_median=0.02613005187787895,height_range={0.02613005187787895,0.02613005187787895}]:0.04745903808983909)S3[&height_95%HPD={0.07358908996771804,0.07358908996771804},height_mean=0.07358908996771804,height_median=0.07358908996771804,height_range={0.07358908996771804,0.07358908996771804}]:0.12857368487312532,((((P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0017728763570904749,O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0017728763570904749)S11[&height_95%HPD={0.0017728763570904749,0.0017728763570904749},height_mean=0.0017728763570904749,height_median=0.0017728763570904749,height_range={0.0017728763570904749,0.0017728763570904749}]:0.044573871489949984,(N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007180867541767444)#H3[&gamma=0.06790226188489634,gamma_95%HPD={0.06790226188489634,0.06790226188489634},gamma_mean=0.06790226188489634,gamma_median=0.06790226188489634,gamma_range={0.06790226188489634,0.06790226188489634},height_95%HPD={0.007180867541767444,0.007180867541767444},height_mean=0.007180867541767444,height_median=0.007180867541767444,height_range={0.007180867541767444,0.007180867541767444}]:0.039165880305273015)S10[&height_95%HPD={0.04634674784704046,0.04634674784704046},height_mean=0.04634674784704046,height_median=0.04634674784704046,height_range={0.04634674784704046,0.04634674784704046}]:0.014216548271922619,(#H1[&height_95%HPD={0.017178654766260792,0.017178654766260792},height_mean=0.017178654766260792,height_median=0.017178654766260792,height_range={0.017178654766260792,0.017178654766260792}]:0.03054142430926321,((M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0105777613246191)#H4[&gamma=0.9030942968747303,gamma_95%HPD={0.9030942968747303,0.9030942968747303},gamma_mean=0.9030942968747303,gamma_median=0.9030942968747303,gamma_range={0.9030942968747303,0.9030942968747303},height_95%HPD={0.0105777613246191,0.0105777613246191},height_mean=0.0105777613246191,height_median=0.0105777613246191,height_range={0.0105777613246191,0.0105777613246191}]:0.013993245564295348,H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02457100688891445)S13[&height_95%HPD={0.02457100688891445,0.02457100688891445},height_mean=0.02457100688891445,height_median=0.02457100688891445,height_range={0.02457100688891445,0.02457100688891445}]:0.023149072186609554)S12[&height_95%HPD={0.047720079075524,0.047720079075524},height_mean=0.047720079075524,height_median=0.047720079075524,height_range={0.047720079075524,0.047720079075524}]:0.012843217043439076)S9[&height_95%HPD={0.06056329611896308,0.06056329611896308},height_mean=0.06056329611896308,height_median=0.06056329611896308,height_range={0.06056329611896308,0.06056329611896308}]:0.05670591702727815,(D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.04622383233943289)#H2[&gamma=0.562406767659462,gamma_95%HPD={0.562406767659462,0.562406767659462},gamma_mean=0.562406767659462,gamma_median=0.562406767659462,gamma_range={0.562406767659462,0.562406767659462},height_95%HPD={0.04622383233943289,0.04622383233943289},height_mean=0.04622383233943289,height_median=0.04622383233943289,height_range={0.04622383233943289,0.04622383233943289}]:0.07104538080680833)S8[&height_95%HPD={0.11726921314624122,0.11726921314624122},height_mean=0.11726921314624122,height_median=0.11726921314624122,height_range={0.11726921314624122,0.11726921314624122}]:0.08489356169460213)S2[&height_95%HPD={0.20216277484084336,0.20216277484084336},height_mean=0.20216277484084336,height_median=0.20216277484084336,height_range={0.20216277484084336,0.20216277484084336}]:0.03352162617529264,(((#H2[&height_95%HPD={0.04622383233943289,0.04622383233943289},height_mean=0.04622383233943289,height_median=0.04622383233943289,height_range={0.04622383233943289,0.04622383233943289}]:0.025966589669454115,((E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.028071716137723163,#H3[&height_95%HPD={0.007180867541767444,0.007180867541767444},height_mean=0.007180867541767444,height_median=0.007180867541767444,height_range={0.007180867541767444,0.007180867541767444}]:0.02089084859595572)S18[&height_95%HPD={0.028071716137723163,0.028071716137723163},height_mean=0.028071716137723163,height_median=0.028071716137723163,height_range={0.028071716137723163,0.028071716137723163}]:0.04206951902203401,(#H4[&height_95%HPD={0.0105777613246191,0.0105777613246191},height_mean=0.0105777613246191,height_median=0.0105777613246191,height_range={0.0105777613246191,0.0105777613246191}]:0.041855261911728586,C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.05243302323634769)S19[&height_95%HPD={0.05243302323634769,0.05243302323634769},height_mean=0.05243302323634769,height_median=0.05243302323634769,height_range={0.05243302323634769,0.05243302323634769}]:0.01770821192340949)S17[&height_95%HPD={0.07014123515975718,0.07014123515975718},height_mean=0.07014123515975718,height_median=0.07014123515975718,height_range={0.07014123515975718,0.07014123515975718}]:0.0020491868491298304)S16[&height_95%HPD={0.072190422008887,0.072190422008887},height_mean=0.072190422008887,height_median=0.072190422008887,height_range={0.072190422008887,0.072190422008887}]:0.011064718658908207,B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.08325514066779521)S15[&height_95%HPD={0.08325514066779521,0.08325514066779521},height_mean=0.08325514066779521,height_median=0.08325514066779521,height_range={0.08325514066779521,0.08325514066779521}]:0.09553959109631066,A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.17879473176410587)S14[&height_95%HPD={0.17879473176410587,0.17879473176410587},height_mean=0.17879473176410587,height_median=0.17879473176410587,height_range={0.17879473176410587,0.17879473176410587}]:0.05688966925203012)S1[&height_95%HPD={0.235684401016136,0.235684401016136},height_mean=0.235684401016136,height_median=0.235684401016136,height_range={0.235684401016136,0.235684401016136}]:0.014315598983864003)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#13", "((((((J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.014426863620686431)#H1[&gamma=0.07376259016304232,gamma_95%HPD={0.07376259016304232,0.07376259016304232},gamma_mean=0.07376259016304232,gamma_median=0.07376259016304232,gamma_range={0.07376259016304232,0.07376259016304232},height_95%HPD={0.014426863620686431,0.014426863620686431},height_mean=0.014426863620686431,height_median=0.014426863620686431,height_range={0.014426863620686431,0.014426863620686431}]:0.016407910272920018,((P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0012554609009508555,O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0012554609009508555)S6[&height_95%HPD={0.0012554609009508555,0.0012554609009508555},height_mean=0.0012554609009508555,height_median=0.0012554609009508555,height_range={0.0012554609009508555,0.0012554609009508555}]:0.006311167053107941,L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007566627954058797)S5[&height_95%HPD={0.007566627954058797,0.007566627954058797},height_mean=0.007566627954058797,height_median=0.007566627954058797,height_range={0.007566627954058797,0.007566627954058797}]:0.023268145939547652)S4[&height_95%HPD={0.03083477389360645,0.03083477389360645},height_mean=0.03083477389360645,height_median=0.03083477389360645,height_range={0.03083477389360645,0.03083477389360645}]:0.08645162496055112,(B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.07367087920214882)#H2[&gamma=0.35466326532669745,gamma_95%HPD={0.35466326532669745,0.35466326532669745},gamma_mean=0.35466326532669745,gamma_median=0.35466326532669745,gamma_range={0.35466326532669745,0.35466326532669745},height_95%HPD={0.07367087920214882,0.07367087920214882},height_mean=0.07367087920214882,height_median=0.07367087920214882,height_range={0.07367087920214882,0.07367087920214882}]:0.043615519652008744)S3[&height_95%HPD={0.11728639885415756,0.11728639885415756},height_mean=0.11728639885415756,height_median=0.11728639885415756,height_range={0.11728639885415756,0.11728639885415756}]:0.009443847082525464,(((N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007205156252452688,M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007205156252452688)S9[&height_95%HPD={0.007205156252452688,0.007205156252452688},height_mean=0.007205156252452688,height_median=0.007205156252452688,height_range={0.007205156252452688,0.007205156252452688}]:0.009748888176151521,(I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.015152913270436802)#H3[&gamma=0.27352536739602384,gamma_95%HPD={0.27352536739602384,0.27352536739602384},gamma_mean=0.27352536739602384,gamma_median=0.27352536739602384,gamma_range={0.27352536739602384,0.27352536739602384},height_95%HPD={0.015152913270436802,0.015152913270436802},height_mean=0.015152913270436802,height_median=0.015152913270436802,height_range={0.015152913270436802,0.015152913270436802}]:0.0018011311581674072)S8[&height_95%HPD={0.01695404442860421,0.01695404442860421},height_mean=0.01695404442860421,height_median=0.01695404442860421,height_range={0.01695404442860421,0.01695404442860421}]:0.10886129081053952,((#H3[&height_95%HPD={0.015152913270436802,0.015152913270436802},height_mean=0.015152913270436802,height_median=0.015152913270436802,height_range={0.015152913270436802,0.015152913270436802}]:0.06663209179387383)#H4[&gamma=0.8512880466954317,gamma_95%HPD={0.8512880466954317,0.8512880466954317},gamma_mean=0.8512880466954317,gamma_median=0.8512880466954317,gamma_range={0.8512880466954317,0.8512880466954317},height_95%HPD={0.08178500506431063,0.08178500506431063},height_mean=0.08178500506431063,height_median=0.08178500506431063,height_range={0.08178500506431063,0.08178500506431063}]:0.040273765506160186,((A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.08660203619052231,((D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0439485114180293,C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0439485114180293)S14[&height_95%HPD={0.0439485114180293,0.0439485114180293},height_mean=0.0439485114180293,height_median=0.0439485114180293,height_range={0.0439485114180293,0.0439485114180293}]:0.03237774964410378,#H2[&height_95%HPD={0.07367087920214882,0.07367087920214882},height_mean=0.07367087920214882,height_median=0.07367087920214882,height_range={0.07367087920214882,0.07367087920214882}]:0.0026553818599842627)S13[&height_95%HPD={0.07632626106213308,0.07632626106213308},height_mean=0.07632626106213308,height_median=0.07632626106213308,height_range={0.07632626106213308,0.07632626106213308}]:0.01027577512838923)S12[&height_95%HPD={0.08660203619052231,0.08660203619052231},height_mean=0.08660203619052231,height_median=0.08660203619052231,height_range={0.08660203619052231,0.08660203619052231}]:0.031070456260878482,(#H4[&height_95%HPD={0.08178500506431063,0.08178500506431063},height_mean=0.08178500506431063,height_median=0.08178500506431063,height_range={0.08178500506431063,0.08178500506431063}]:0.017939250589755956,(#H1[&height_95%HPD={0.014426863620686431,0.014426863620686431},height_mean=0.014426863620686431,height_median=0.014426863620686431,height_range={0.014426863620686431,0.014426863620686431}]:0.04078472644090045,(E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03010110873702812,F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03010110873702812)S17[&height_95%HPD={0.03010110873702812,0.03010110873702812},height_mean=0.03010110873702812,height_median=0.03010110873702812,height_range={0.03010110873702812,0.03010110873702812}]:0.02511048132455876)S16[&height_95%HPD={0.05521159006158688,0.05521159006158688},height_mean=0.05521159006158688,height_median=0.05521159006158688,height_range={0.05521159006158688,0.05521159006158688}]:0.04451266559247971)S15[&height_95%HPD={0.09972425565406659,0.09972425565406659},height_mean=0.09972425565406659,height_median=0.09972425565406659,height_range={0.09972425565406659,0.09972425565406659}]:0.017948236797334205)S11[&height_95%HPD={0.1176724924514008,0.1176724924514008},height_mean=0.1176724924514008,height_median=0.1176724924514008,height_range={0.1176724924514008,0.1176724924514008}]:0.004386278119070025)S10[&height_95%HPD={0.12205877057047082,0.12205877057047082},height_mean=0.12205877057047082,height_median=0.12205877057047082,height_range={0.12205877057047082,0.12205877057047082}]:0.0037565646686729126)S7[&height_95%HPD={0.12581533523914373,0.12581533523914373},height_mean=0.12581533523914373,height_median=0.12581533523914373,height_range={0.12581533523914373,0.12581533523914373}]:9.149106975392951E-4)S2[&height_95%HPD={0.12673024593668303,0.12673024593668303},height_mean=0.12673024593668303,height_median=0.12673024593668303,height_range={0.12673024593668303,0.12673024593668303}]:0.060685228483116804,((H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01555672299259242,G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01555672299259242)S19[&height_95%HPD={0.01555672299259242,0.01555672299259242},height_mean=0.01555672299259242,height_median=0.01555672299259242,height_range={0.01555672299259242,0.01555672299259242}]:0.08718575593129829,K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.10274247892389071)S18[&height_95%HPD={0.10274247892389071,0.10274247892389071},height_mean=0.10274247892389071,height_median=0.10274247892389071,height_range={0.10274247892389071,0.10274247892389071}]:0.08467299549590912)S1[&height_95%HPD={0.18741547441979983,0.18741547441979983},height_mean=0.18741547441979983,height_median=0.18741547441979983,height_range={0.18741547441979983,0.18741547441979983}]:0.06258452558020017)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#12", "((((((O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.004579762788431058,P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.004579762788431058)S15[&height_95%HPD={0.004579762788431058,0.004579762788431058},height_mean=0.004579762788431058,height_median=0.004579762788431058,height_range={0.004579762788431058,0.004579762788431058}]:0.05148779330917733)#H3[&gamma=0.972101314480719,gamma_95%HPD={0.972101314480719,0.972101314480719},gamma_mean=0.972101314480719,gamma_median=0.972101314480719,gamma_range={0.972101314480719,0.972101314480719},height_95%HPD={0.056067556097608384,0.056067556097608384},height_mean=0.056067556097608384,height_median=0.056067556097608384,height_range={0.056067556097608384,0.056067556097608384}]:0.021077489822582784,((H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.020661152728918497,I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.020661152728918497)S12[&height_95%HPD={0.020661152728918497,0.020661152728918497},height_mean=0.020661152728918497,height_median=0.020661152728918497,height_range={0.020661152728918497,0.020661152728918497}]:0.041399215626152425)#H4[&gamma=0.9823864589615545,gamma_95%HPD={0.9823864589615545,0.9823864589615545},gamma_mean=0.9823864589615545,gamma_median=0.9823864589615545,gamma_range={0.9823864589615545,0.9823864589615545},height_95%HPD={0.06206036835507092,0.06206036835507092},height_mean=0.06206036835507092,height_median=0.06206036835507092,height_range={0.06206036835507092,0.06206036835507092}]:0.015084677565120247)S6[&height_95%HPD={0.07714504592019117,0.07714504592019117},height_mean=0.07714504592019117,height_median=0.07714504592019117,height_range={0.07714504592019117,0.07714504592019117}]:0.06200129312206876,((#H4[&height_95%HPD={0.06206036835507092,0.06206036835507092},height_mean=0.06206036835507092,height_median=0.06206036835507092,height_range={0.06206036835507092,0.06206036835507092}]:0.04325782720638707,(((E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.028824316149344137)#H1[&gamma=0.7349188608241479,gamma_95%HPD={0.7349188608241479,0.7349188608241479},gamma_mean=0.7349188608241479,gamma_median=0.7349188608241479,gamma_range={0.7349188608241479,0.7349188608241479},height_95%HPD={0.028824316149344137,0.028824316149344137},height_mean=0.028824316149344137,height_median=0.028824316149344137,height_range={0.028824316149344137,0.028824316149344137}]:0.014155826136100946)#H2[&gamma=0.7949347517734078,gamma_95%HPD={0.7949347517734078,0.7949347517734078},gamma_mean=0.7949347517734078,gamma_median=0.7949347517734078,gamma_range={0.7949347517734078,0.7949347517734078},height_95%HPD={0.04298014228544508,0.04298014228544508},height_mean=0.04298014228544508,height_median=0.04298014228544508,height_range={0.04298014228544508,0.04298014228544508}]:0.005654781222549171,((N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007000256731683002,M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007000256731683002)S19[&height_95%HPD={0.007000256731683002,0.007000256731683002},height_mean=0.007000256731683002,height_median=0.007000256731683002,height_range={0.007000256731683002,0.007000256731683002}]:0.003694085510706574,L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.010694342242389576)S18[&height_95%HPD={0.010694342242389576,0.010694342242389576},height_mean=0.010694342242389576,height_median=0.010694342242389576,height_range={0.010694342242389576,0.010694342242389576}]:0.03794058126560468)S9[&height_95%HPD={0.048634923507994254,0.048634923507994254},height_mean=0.048634923507994254,height_median=0.048634923507994254,height_range={0.048634923507994254,0.048634923507994254}]:0.05668327205346374)S7[&height_95%HPD={0.105318195561458,0.105318195561458},height_mean=0.105318195561458,height_median=0.105318195561458,height_range={0.105318195561458,0.105318195561458}]:0.015281188954862246,A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.12059938451632024)S5[&height_95%HPD={0.12059938451632024,0.12059938451632024},height_mean=0.12059938451632024,height_median=0.12059938451632024,height_range={0.12059938451632024,0.12059938451632024}]:0.018546954525939685)S2[&height_95%HPD={0.13914633904225993,0.13914633904225993},height_mean=0.13914633904225993,height_median=0.13914633904225993,height_range={0.13914633904225993,0.13914633904225993}]:0.030662467294945484,(#H3[&height_95%HPD={0.056067556097608384,0.056067556097608384},height_mean=0.056067556097608384,height_median=0.056067556097608384,height_range={0.056067556097608384,0.056067556097608384}]:0.10538513203525263,((((F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.025684455899571995,G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.025684455899571995)S13[&height_95%HPD={0.025684455899571995,0.025684455899571995},height_mean=0.025684455899571995,height_median=0.025684455899571995,height_range={0.025684455899571995,0.025684455899571995}]:0.035360705370648976,(#H2[&height_95%HPD={0.04298014228544508,0.04298014228544508},height_mean=0.04298014228544508,height_median=0.04298014228544508,height_range={0.04298014228544508,0.04298014228544508}]:0.012441253872856967,((K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.016082994246625737,J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.016082994246625737)S17[&height_95%HPD={0.016082994246625737,0.016082994246625737},height_mean=0.016082994246625737,height_median=0.016082994246625737,height_range={0.016082994246625737,0.016082994246625737}]:0.03791054972353833,C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.053993543970164065)S16[&height_95%HPD={0.053993543970164065,0.053993543970164065},height_mean=0.053993543970164065,height_median=0.053993543970164065,height_range={0.053993543970164065,0.053993543970164065}]:0.0014278521881379846)S14[&height_95%HPD={0.05542139615830205,0.05542139615830205},height_mean=0.05542139615830205,height_median=0.05542139615830205,height_range={0.05542139615830205,0.05542139615830205}]:0.005623765111918921)S10[&height_95%HPD={0.06104516127022097,0.06104516127022097},height_mean=0.06104516127022097,height_median=0.06104516127022097,height_range={0.06104516127022097,0.06104516127022097}]:0.0020576322851856843,(#H1[&height_95%HPD={0.028824316149344137,0.028824316149344137},height_mean=0.028824316149344137,height_median=0.028824316149344137,height_range={0.028824316149344137,0.028824316149344137}]:0.0016644894060449111,D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.030488805555389048)S11[&height_95%HPD={0.030488805555389048,0.030488805555389048},height_mean=0.030488805555389048,height_median=0.030488805555389048,height_range={0.030488805555389048,0.030488805555389048}]:0.03261398800001761)S8[&height_95%HPD={0.06310279355540666,0.06310279355540666},height_mean=0.06310279355540666,height_median=0.06310279355540666,height_range={0.06310279355540666,0.06310279355540666}]:0.04856071939096657,B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.11166351294637322)S4[&height_95%HPD={0.11166351294637322,0.11166351294637322},height_mean=0.11166351294637322,height_median=0.11166351294637322,height_range={0.11166351294637322,0.11166351294637322}]:0.04978917518648779)S3[&height_95%HPD={0.161452688132861,0.161452688132861},height_mean=0.161452688132861,height_median=0.161452688132861,height_range={0.161452688132861,0.161452688132861}]:0.008356118204344398)S1[&height_95%HPD={0.1698088063372054,0.1698088063372054},height_mean=0.1698088063372054,height_median=0.1698088063372054,height_range={0.1698088063372054,0.1698088063372054}]:0.08019119366279459)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#11", "((((G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03965453551189635,(N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.011321231783757157,M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.011321231783757157)S17[&height_95%HPD={0.011321231783757157,0.011321231783757157},height_mean=0.011321231783757157,height_median=0.011321231783757157,height_range={0.011321231783757157,0.011321231783757157}]:0.028333303728139192)S10[&height_95%HPD={0.03965453551189635,0.03965453551189635},height_mean=0.03965453551189635,height_median=0.03965453551189635,height_range={0.03965453551189635,0.03965453551189635}]:0.02844292249568864)#H3[&gamma=0.557755754090305,gamma_95%HPD={0.557755754090305,0.557755754090305},gamma_mean=0.557755754090305,gamma_median=0.557755754090305,gamma_range={0.557755754090305,0.557755754090305},height_95%HPD={0.06809745800758499,0.06809745800758499},height_mean=0.06809745800758499,height_median=0.06809745800758499,height_range={0.06809745800758499,0.06809745800758499}]:0.08044476924407662,((#H3[&height_95%HPD={0.06809745800758499,0.06809745800758499},height_mean=0.06809745800758499,height_median=0.06809745800758499,height_range={0.06809745800758499,0.06809745800758499}]:0.03860081953401345,(((H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.025212301486450206,(J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01630540468626554)#H1[&gamma=0.7968608495940704,gamma_95%HPD={0.7968608495940704,0.7968608495940704},gamma_mean=0.7968608495940704,gamma_median=0.7968608495940704,gamma_range={0.7968608495940704,0.7968608495940704},height_95%HPD={0.01630540468626554,0.01630540468626554},height_mean=0.01630540468626554,height_median=0.01630540468626554,height_range={0.01630540468626554,0.01630540468626554}]:0.008906896800184666)S11[&height_95%HPD={0.025212301486450206,0.025212301486450206},height_mean=0.025212301486450206,height_median=0.025212301486450206,height_range={0.025212301486450206,0.025212301486450206}]:0.04016925405810973,B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.06538155554455993)S8[&height_95%HPD={0.06538155554455993,0.06538155554455993},height_mean=0.06538155554455993,height_median=0.06538155554455993,height_range={0.06538155554455993,0.06538155554455993}]:0.0061076979461861525,(D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.052467426571753184,(I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.016375376490705196,(L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.012908188176587149,K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.012908188176587149)S18[&height_95%HPD={0.012908188176587149,0.012908188176587149},height_mean=0.012908188176587149,height_median=0.012908188176587149,height_range={0.012908188176587149,0.012908188176587149}]:0.003467188314118047)S14[&height_95%HPD={0.016375376490705196,0.016375376490705196},height_mean=0.016375376490705196,height_median=0.016375376490705196,height_range={0.016375376490705196,0.016375376490705196}]:0.03609205008104799)S9[&height_95%HPD={0.052467426571753184,0.052467426571753184},height_mean=0.052467426571753184,height_median=0.052467426571753184,height_range={0.052467426571753184,0.052467426571753184}]:0.0190218269189929)S5[&height_95%HPD={0.07148925349074609,0.07148925349074609},height_mean=0.07148925349074609,height_median=0.07148925349074609,height_range={0.07148925349074609,0.07148925349074609}]:0.035209024050852356)S3[&height_95%HPD={0.10669827754159844,0.10669827754159844},height_mean=0.10669827754159844,height_median=0.10669827754159844,height_range={0.10669827754159844,0.10669827754159844}]:0.022609205693167406,((((O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007223362730762406,P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007223362730762406)S16[&height_95%HPD={0.007223362730762406,0.007223362730762406},height_mean=0.007223362730762406,height_median=0.007223362730762406,height_range={0.007223362730762406,0.007223362730762406}]:0.03415196489371597,F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.041375327624478375)S15[&height_95%HPD={0.041375327624478375,0.041375327624478375},height_mean=0.041375327624478375,height_median=0.041375327624478375,height_range={0.041375327624478375,0.041375327624478375}]:0.0054430854072518575)#H2[&gamma=0.8670717803958126,gamma_95%HPD={0.8670717803958126,0.8670717803958126},gamma_mean=0.8670717803958126,gamma_median=0.8670717803958126,gamma_range={0.8670717803958126,0.8670717803958126},height_95%HPD={0.04681841303173023,0.04681841303173023},height_mean=0.04681841303173023,height_median=0.04681841303173023,height_range={0.04681841303173023,0.04681841303173023}]:0.03602145431783005,(((C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.05322296995707024,(E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.05169191517677535,#H2[&height_95%HPD={0.04681841303173023,0.04681841303173023},height_mean=0.04681841303173023,height_median=0.04681841303173023,height_range={0.04681841303173023,0.04681841303173023}]:0.004873502145045117)S13[&height_95%HPD={0.05169191517677535,0.05169191517677535},height_mean=0.05169191517677535,height_median=0.05169191517677535,height_range={0.05169191517677535,0.05169191517677535}]:0.0015310547802948882)S12[&height_95%HPD={0.05322296995707024,0.05322296995707024},height_mean=0.05322296995707024,height_median=0.05322296995707024,height_range={0.05322296995707024,0.05322296995707024}]:0.00954517820449971,#H1[&height_95%HPD={0.01630540468626554,0.01630540468626554},height_mean=0.01630540468626554,height_median=0.01630540468626554,height_range={0.01630540468626554,0.01630540468626554}]:0.04646274347530441)S7[&height_95%HPD={0.06276814816156995,0.06276814816156995},height_mean=0.06276814816156995,height_median=0.06276814816156995,height_range={0.06276814816156995,0.06276814816156995}]:0.01599373973185403,A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.07876188789342398)S6[&height_95%HPD={0.07876188789342398,0.07876188789342398},height_mean=0.07876188789342398,height_median=0.07876188789342398,height_range={0.07876188789342398,0.07876188789342398}]:0.004077979456136305)S4[&height_95%HPD={0.08283986734956028,0.08283986734956028},height_mean=0.08283986734956028,height_median=0.08283986734956028,height_range={0.08283986734956028,0.08283986734956028}]:0.046467615885205565)S2[&height_95%HPD={0.12930748323476585,0.12930748323476585},height_mean=0.12930748323476585,height_median=0.12930748323476585,height_range={0.12930748323476585,0.12930748323476585}]:0.01923474401689576)S1[&height_95%HPD={0.1485422272516616,0.1485422272516616},height_mean=0.1485422272516616,height_median=0.1485422272516616,height_range={0.1485422272516616,0.1485422272516616}]:0.10145777274833839)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#10", "((((K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0184486472566881,L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0184486472566881)S5[&height_95%HPD={0.0184486472566881,0.0184486472566881},height_mean=0.0184486472566881,height_median=0.0184486472566881,height_range={0.0184486472566881,0.0184486472566881}]:0.10206208660187102,(B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.10709290651270595,(((C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.07391224459918935,(H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.023413037330219544,G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.023413037330219544)S11[&height_95%HPD={0.023413037330219544,0.023413037330219544},height_mean=0.023413037330219544,height_median=0.023413037330219544,height_range={0.023413037330219544,0.023413037330219544}]:0.05049920726896981)S10[&height_95%HPD={0.07391224459918935,0.07391224459918935},height_mean=0.07391224459918935,height_median=0.07391224459918935,height_range={0.07391224459918935,0.07391224459918935}]:0.00438433740875388,(((I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.018462554508548168,J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.018462554508548168)S17[&height_95%HPD={0.018462554508548168,0.018462554508548168},height_mean=0.018462554508548168,height_median=0.018462554508548168,height_range={0.018462554508548168,0.018462554508548168}]:0.008005040517736817)#H2[&gamma=0.28278096579039347,gamma_95%HPD={0.28278096579039347,0.28278096579039347},gamma_mean=0.28278096579039347,gamma_median=0.28278096579039347,gamma_range={0.28278096579039347,0.28278096579039347},height_95%HPD={0.026467595026284985,0.026467595026284985},height_mean=0.026467595026284985,height_median=0.026467595026284985,height_range={0.026467595026284985,0.026467595026284985}]:0.03850849215637342,((((P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.010237725762182737,O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.010237725762182737)S18[&height_95%HPD={0.010237725762182737,0.010237725762182737},height_mean=0.010237725762182737,height_median=0.010237725762182737,height_range={0.010237725762182737,0.010237725762182737}]:0.006026926983616954,N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01626465274579969)S16[&height_95%HPD={0.01626465274579969,0.01626465274579969},height_mean=0.01626465274579969,height_median=0.01626465274579969,height_range={0.01626465274579969,0.01626465274579969}]:0.011080325302090743,#H2[&height_95%HPD={0.026467595026284985,0.026467595026284985},height_mean=0.026467595026284985,height_median=0.026467595026284985,height_range={0.026467595026284985,0.026467595026284985}]:8.77383021605449E-4)S15[&height_95%HPD={0.027344978047890434,0.027344978047890434},height_mean=0.027344978047890434,height_median=0.027344978047890434,height_range={0.027344978047890434,0.027344978047890434}]:0.007298920938904191,E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.034643898986794625)S12[&height_95%HPD={0.034643898986794625,0.034643898986794625},height_mean=0.034643898986794625,height_median=0.034643898986794625,height_range={0.034643898986794625,0.034643898986794625}]:0.03033218819586378)S9[&height_95%HPD={0.0649760871826584,0.0649760871826584},height_mean=0.0649760871826584,height_median=0.0649760871826584,height_range={0.0649760871826584,0.0649760871826584}]:0.013320494825284829)S8[&height_95%HPD={0.07829658200794323,0.07829658200794323},height_mean=0.07829658200794323,height_median=0.07829658200794323,height_range={0.07829658200794323,0.07829658200794323}]:0.0062520570313351065)#H3[&gamma=0.7229921048849862,gamma_95%HPD={0.7229921048849862,0.7229921048849862},gamma_mean=0.7229921048849862,gamma_median=0.7229921048849862,gamma_range={0.7229921048849862,0.7229921048849862},height_95%HPD={0.08454863903927834,0.08454863903927834},height_mean=0.08454863903927834,height_median=0.08454863903927834,height_range={0.08454863903927834,0.08454863903927834}]:0.02254426747342761)S6[&height_95%HPD={0.10709290651270595,0.10709290651270595},height_mean=0.10709290651270595,height_median=0.10709290651270595,height_range={0.10709290651270595,0.10709290651270595}]:0.013417827345853173)S2[&height_95%HPD={0.12051073385855912,0.12051073385855912},height_mean=0.12051073385855912,height_median=0.12051073385855912,height_range={0.12051073385855912,0.12051073385855912}]:0.1254897890755196,((((D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.037031303338413146,((M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.017421068235199066)#H1[&gamma=0.209498304961934,gamma_95%HPD={0.209498304961934,0.209498304961934},gamma_mean=0.209498304961934,gamma_median=0.209498304961934,gamma_range={0.209498304961934,0.209498304961934},height_95%HPD={0.017421068235199066,0.017421068235199066},height_mean=0.017421068235199066,height_median=0.017421068235199066,height_range={0.017421068235199066,0.017421068235199066}]:0.01564442374735328,F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03306549198255235)S14[&height_95%HPD={0.03306549198255235,0.03306549198255235},height_mean=0.03306549198255235,height_median=0.03306549198255235,height_range={0.03306549198255235,0.03306549198255235}]:0.003965811355860799)S13[&height_95%HPD={0.037031303338413146,0.037031303338413146},height_mean=0.037031303338413146,height_median=0.037031303338413146,height_range={0.037031303338413146,0.037031303338413146}]:0.02114818353933981,#H1[&height_95%HPD={0.017421068235199066,0.017421068235199066},height_mean=0.017421068235199066,height_median=0.017421068235199066,height_range={0.017421068235199066,0.017421068235199066}]:0.04075841864255389)S7[&height_95%HPD={0.058179486877752956,0.058179486877752956},height_mean=0.058179486877752956,height_median=0.058179486877752956,height_range={0.058179486877752956,0.058179486877752956}]:0.030579248178853297,#H3[&height_95%HPD={0.08454863903927834,0.08454863903927834},height_mean=0.08454863903927834,height_median=0.08454863903927834,height_range={0.08454863903927834,0.08454863903927834}]:0.004210096017327913)S4[&height_95%HPD={0.08875873505660625,0.08875873505660625},height_mean=0.08875873505660625,height_median=0.08875873505660625,height_range={0.08875873505660625,0.08875873505660625}]:0.041873680233158767,A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.13063241528976502)S3[&height_95%HPD={0.13063241528976502,0.13063241528976502},height_mean=0.13063241528976502,height_median=0.13063241528976502,height_range={0.13063241528976502,0.13063241528976502}]:0.11536810764431371)S1[&height_95%HPD={0.24600052293407873,0.24600052293407873},height_mean=0.24600052293407873,height_median=0.24600052293407873,height_range={0.24600052293407873,0.24600052293407873}]:0.00399947706592127)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#9", "(((((((G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.027465330951573663)#H3[&gamma=0.366705248579973,gamma_95%HPD={0.366705248579973,0.366705248579973},gamma_mean=0.366705248579973,gamma_median=0.366705248579973,gamma_range={0.366705248579973,0.366705248579973},height_95%HPD={0.027465330951573663,0.027465330951573663},height_mean=0.027465330951573663,height_median=0.027465330951573663,height_range={0.027465330951573663,0.027465330951573663}]:0.05019128020513522,(O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0018641472034237327,P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0018641472034237327)S11[&height_95%HPD={0.0018641472034237327,0.0018641472034237327},height_mean=0.0018641472034237327,height_median=0.0018641472034237327,height_range={0.0018641472034237327,0.0018641472034237327}]:0.07579246395328515)S9[&height_95%HPD={0.07765661115670888,0.07765661115670888},height_mean=0.07765661115670888,height_median=0.07765661115670888,height_range={0.07765661115670888,0.07765661115670888}]:0.00438947410461793,((H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.022988627042481724,(I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02125952954997723,((L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.008451719369610922,K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.008451719369610922)S17[&height_95%HPD={0.008451719369610922,0.008451719369610922},height_mean=0.008451719369610922,height_median=0.008451719369610922,height_range={0.008451719369610922,0.008451719369610922}]:0.009849077322158784)#H2[&gamma=0.8260753572175957,gamma_95%HPD={0.8260753572175957,0.8260753572175957},gamma_mean=0.8260753572175957,gamma_median=0.8260753572175957,gamma_range={0.8260753572175957,0.8260753572175957},height_95%HPD={0.018300796691769705,0.018300796691769705},height_mean=0.018300796691769705,height_median=0.018300796691769705,height_range={0.018300796691769705,0.018300796691769705}]:0.0029587328582075245)S16[&height_95%HPD={0.02125952954997723,0.02125952954997723},height_mean=0.02125952954997723,height_median=0.02125952954997723,height_range={0.02125952954997723,0.02125952954997723}]:0.0017290974925044944)S10[&height_95%HPD={0.022988627042481724,0.022988627042481724},height_mean=0.022988627042481724,height_median=0.022988627042481724,height_range={0.022988627042481724,0.022988627042481724}]:0.056254307339719783,B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.07924293438220151)S8[&height_95%HPD={0.07924293438220151,0.07924293438220151},height_mean=0.07924293438220151,height_median=0.07924293438220151,height_range={0.07924293438220151,0.07924293438220151}]:0.0028031508791253046)S7[&height_95%HPD={0.08204608526132681,0.08204608526132681},height_mean=0.08204608526132681,height_median=0.08204608526132681,height_range={0.08204608526132681,0.08204608526132681}]:0.025540201264476647,#H2[&height_95%HPD={0.018300796691769705,0.018300796691769705},height_mean=0.018300796691769705,height_median=0.018300796691769705,height_range={0.018300796691769705,0.018300796691769705}]:0.08928548983403375)S4[&height_95%HPD={0.10758628652580346,0.10758628652580346},height_mean=0.10758628652580346,height_median=0.10758628652580346,height_range={0.10758628652580346,0.10758628652580346}]:0.029645750090363826,(((E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02986178300988787,F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02986178300988787)S15[&height_95%HPD={0.02986178300988787,0.02986178300988787},height_mean=0.02986178300988787,height_median=0.02986178300988787,height_range={0.02986178300988787,0.02986178300988787}]:0.012056948678513452,C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.04191873168840132)S6[&height_95%HPD={0.04191873168840132,0.04191873168840132},height_mean=0.04191873168840132,height_median=0.04191873168840132,height_range={0.04191873168840132,0.04191873168840132}]:0.07717753265080574,A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.11909626433920706)S5[&height_95%HPD={0.11909626433920706,0.11909626433920706},height_mean=0.11909626433920706,height_median=0.11909626433920706,height_range={0.11909626433920706,0.11909626433920706}]:0.018135772276960227)S3[&height_95%HPD={0.13723203661616729,0.13723203661616729},height_mean=0.13723203661616729,height_median=0.13723203661616729,height_range={0.13723203661616729,0.13723203661616729}]:0.06991292816644656,((N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0035316640245187936)#H1[&gamma=0.03648192657551197,gamma_95%HPD={0.03648192657551197,0.03648192657551197},gamma_mean=0.03648192657551197,gamma_median=0.03648192657551197,gamma_range={0.03648192657551197,0.03648192657551197},height_95%HPD={0.0035316640245187936,0.0035316640245187936},height_mean=0.0035316640245187936,height_median=0.0035316640245187936,height_range={0.0035316640245187936,0.0035316640245187936}]:0.06562660548385671,(((#H1[&height_95%HPD={0.0035316640245187936,0.0035316640245187936},height_mean=0.0035316640245187936,height_median=0.0035316640245187936,height_range={0.0035316640245187936,0.0035316640245187936}]:7.985883312095488E-4,M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.004330252355728342)S18[&height_95%HPD={0.004330252355728342,0.004330252355728342},height_mean=0.004330252355728342,height_median=0.004330252355728342,height_range={0.004330252355728342,0.004330252355728342}]:0.011038974804094964,J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.015369227159823307)S13[&height_95%HPD={0.015369227159823307,0.015369227159823307},height_mean=0.015369227159823307,height_median=0.015369227159823307,height_range={0.015369227159823307,0.015369227159823307}]:0.03815510267226424,(#H3[&height_95%HPD={0.027465330951573663,0.027465330951573663},height_mean=0.027465330951573663,height_median=0.027465330951573663,height_range={0.027465330951573663,0.027465330951573663}]:0.0024225786514855963,D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02988790960305926)S14[&height_95%HPD={0.02988790960305926,0.02988790960305926},height_mean=0.02988790960305926,height_median=0.02988790960305926,height_range={0.02988790960305926,0.02988790960305926}]:0.023636420229028288)S12[&height_95%HPD={0.05352432983208755,0.05352432983208755},height_mean=0.05352432983208755,height_median=0.05352432983208755,height_range={0.05352432983208755,0.05352432983208755}]:0.015633939676287956)S2[&height_95%HPD={0.0691582695083755,0.0691582695083755},height_mean=0.0691582695083755,height_median=0.0691582695083755,height_range={0.0691582695083755,0.0691582695083755}]:0.13798669527423835)S1[&height_95%HPD={0.20714496478261385,0.20714496478261385},height_mean=0.20714496478261385,height_median=0.20714496478261385,height_range={0.20714496478261385,0.20714496478261385}]:0.042855035217386206)[&height_95%HPD={0.25000000000000006,0.25000000000000006},height_mean=0.25000000000000006,height_median=0.25000000000000006,height_range={0.25000000000000006,0.25000000000000006},topologySupport=1.0E-4];");
//        beastStrings.put("#8", "((((((P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.005324571168481895)#H1[&gamma=0.0633977492675657,gamma_95%HPD={0.0633977492675657,0.0633977492675657},gamma_mean=0.0633977492675657,gamma_median=0.0633977492675657,gamma_range={0.0633977492675657,0.0633977492675657},height_95%HPD={0.005324571168481895,0.005324571168481895},height_mean=0.005324571168481895,height_median=0.005324571168481895,height_range={0.005324571168481895,0.005324571168481895}]:0.032068566233910106,(O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007909852394431177,N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007909852394431177)S17[&height_95%HPD={0.007909852394431177,0.007909852394431177},height_mean=0.007909852394431177,height_median=0.007909852394431177,height_range={0.007909852394431177,0.007909852394431177}]:0.029483285007960824)S13[&height_95%HPD={0.037393137402392,0.037393137402392},height_mean=0.037393137402392,height_median=0.037393137402392,height_range={0.037393137402392,0.037393137402392}]:0.030801382321658133,B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.06819451972405013)S4[&height_95%HPD={0.06819451972405013,0.06819451972405013},height_mean=0.06819451972405013,height_median=0.06819451972405013,height_range={0.06819451972405013,0.06819451972405013}]:0.06394096272185265,((A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.06950255655749554,((D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.04075228697270514,(I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.019593050347006646)#H2[&gamma=0.0854092089981272,gamma_95%HPD={0.0854092089981272,0.0854092089981272},gamma_mean=0.0854092089981272,gamma_median=0.0854092089981272,gamma_range={0.0854092089981272,0.0854092089981272},height_95%HPD={0.019593050347006646,0.019593050347006646},height_mean=0.019593050347006646,height_median=0.019593050347006646,height_range={0.019593050347006646,0.019593050347006646}]:0.02115923662569849)S14[&height_95%HPD={0.04075228697270514,0.04075228697270514},height_mean=0.04075228697270514,height_median=0.04075228697270514,height_range={0.04075228697270514,0.04075228697270514}]:0.02533443512506986,#H1[&height_95%HPD={0.005324571168481895,0.005324571168481895},height_mean=0.005324571168481895,height_median=0.005324571168481895,height_range={0.005324571168481895,0.005324571168481895}]:0.0607621509292931)S10[&height_95%HPD={0.066086722097775,0.066086722097775},height_mean=0.066086722097775,height_median=0.066086722097775,height_range={0.066086722097775,0.066086722097775}]:0.003415834459720546)S7[&height_95%HPD={0.06950255655749554,0.06950255655749554},height_mean=0.06950255655749554,height_median=0.06950255655749554,height_range={0.06950255655749554,0.06950255655749554}]:0.015502649356816861,((#H2[&height_95%HPD={0.019593050347006646,0.019593050347006646},height_mean=0.019593050347006646,height_median=0.019593050347006646,height_range={0.019593050347006646,0.019593050347006646}]:0.018700382866204485,(H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.026243453916657083,G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.026243453916657083)S16[&height_95%HPD={0.026243453916657083,0.026243453916657083},height_mean=0.026243453916657083,height_median=0.026243453916657083,height_range={0.026243453916657083,0.026243453916657083}]:0.012049979296554048)S9[&height_95%HPD={0.03829343321321113,0.03829343321321113},height_mean=0.03829343321321113,height_median=0.03829343321321113,height_range={0.03829343321321113,0.03829343321321113}]:0.04431178974779307,((C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0463467074552916,(J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01636923524114528,K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01636923524114528)S15[&height_95%HPD={0.01636923524114528,0.01636923524114528},height_mean=0.01636923524114528,height_median=0.01636923524114528,height_range={0.01636923524114528,0.01636923524114528}]:0.029977472214146317)S12[&height_95%HPD={0.0463467074552916,0.0463467074552916},height_mean=0.0463467074552916,height_median=0.0463467074552916,height_range={0.0463467074552916,0.0463467074552916}]:0.02195685149306989,(E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03252142351467635,F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03252142351467635)S11[&height_95%HPD={0.03252142351467635,0.03252142351467635},height_mean=0.03252142351467635,height_median=0.03252142351467635,height_range={0.03252142351467635,0.03252142351467635}]:0.035782135433685136)S8[&height_95%HPD={0.06830355894836149,0.06830355894836149},height_mean=0.06830355894836149,height_median=0.06830355894836149,height_range={0.06830355894836149,0.06830355894836149}]:0.014301664012642712)S6[&height_95%HPD={0.0826052229610042,0.0826052229610042},height_mean=0.0826052229610042,height_median=0.0826052229610042,height_range={0.0826052229610042,0.0826052229610042}]:0.002399982953308205)S5[&height_95%HPD={0.0850052059143124,0.0850052059143124},height_mean=0.0850052059143124,height_median=0.0850052059143124,height_range={0.0850052059143124,0.0850052059143124}]:0.04713027653159038)S3[&height_95%HPD={0.13213548244590279,0.13213548244590279},height_mean=0.13213548244590279,height_median=0.13213548244590279,height_range={0.13213548244590279,0.13213548244590279}]:0.054582382646049116,(M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.009208708241443542,L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.009208708241443542)S2[&height_95%HPD={0.009208708241443542,0.009208708241443542},height_mean=0.009208708241443542,height_median=0.009208708241443542,height_range={0.009208708241443542,0.009208708241443542}]:0.17750915685050836)S1[&height_95%HPD={0.1867178650919519,0.1867178650919519},height_mean=0.1867178650919519,height_median=0.1867178650919519,height_range={0.1867178650919519,0.1867178650919519}]:0.0632821349080481)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#7", "(((((P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:1.4569942828179805E-4,O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:1.4569942828179805E-4)S17[&height_95%HPD={1.4569942828179805E-4,1.4569942828179805E-4},height_mean=1.4569942828179805E-4,height_median=1.4569942828179805E-4,height_range={1.4569942828179805E-4,1.4569942828179805E-4}]:0.014390705996872255,J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.014536405425154053)S8[&height_95%HPD={0.014536405425154053,0.014536405425154053},height_mean=0.014536405425154053,height_median=0.014536405425154053,height_range={0.014536405425154053,0.014536405425154053}]:0.0966955323990436,(B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.05653620742881535,((G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03358447303437434)#H1[&gamma=0.8255783017242199,gamma_95%HPD={0.8255783017242199,0.8255783017242199},gamma_mean=0.8255783017242199,gamma_median=0.8255783017242199,gamma_range={0.8255783017242199,0.8255783017242199},height_95%HPD={0.03358447303437434,0.03358447303437434},height_mean=0.03358447303437434,height_median=0.03358447303437434,height_range={0.03358447303437434,0.03358447303437434}]:0.015868052400132532,((H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.029475698245895132,I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.029475698245895132)S16[&height_95%HPD={0.029475698245895132,0.029475698245895132},height_mean=0.029475698245895132,height_median=0.029475698245895132,height_range={0.029475698245895132,0.029475698245895132}]:0.012010917801460824,(L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.009103444547421902,K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.009103444547421902)S15[&height_95%HPD={0.009103444547421902,0.009103444547421902},height_mean=0.009103444547421902,height_median=0.009103444547421902,height_range={0.009103444547421902,0.009103444547421902}]:0.032383171499934055)S14[&height_95%HPD={0.041486616047355956,0.041486616047355956},height_mean=0.041486616047355956,height_median=0.041486616047355956,height_range={0.041486616047355956,0.041486616047355956}]:0.007965909387150916)S13[&height_95%HPD={0.04945252543450687,0.04945252543450687},height_mean=0.04945252543450687,height_median=0.04945252543450687,height_range={0.04945252543450687,0.04945252543450687}]:0.007083681994308477)S7[&height_95%HPD={0.05653620742881535,0.05653620742881535},height_mean=0.05653620742881535,height_median=0.05653620742881535,height_range={0.05653620742881535,0.05653620742881535}]:0.05469573039538231)S2[&height_95%HPD={0.11123193782419766,0.11123193782419766},height_mean=0.11123193782419766,height_median=0.11123193782419766,height_range={0.11123193782419766,0.11123193782419766}]:0.09842758494214099,(((D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.05183517456437342,C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.05183517456437342)S12[&height_95%HPD={0.05183517456437342,0.05183517456437342},height_mean=0.05183517456437342,height_median=0.05183517456437342,height_range={0.05183517456437342,0.05183517456437342}]:0.007704798081034553,(F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.034771350355517344)#H2[&gamma=0.554831845491623,gamma_95%HPD={0.554831845491623,0.554831845491623},gamma_mean=0.554831845491623,gamma_median=0.554831845491623,gamma_range={0.554831845491623,0.554831845491623},height_95%HPD={0.034771350355517344,0.034771350355517344},height_mean=0.034771350355517344,height_median=0.034771350355517344,height_range={0.034771350355517344,0.034771350355517344}]:0.02476862228989063)S5[&height_95%HPD={0.059539972645407974,0.059539972645407974},height_mean=0.059539972645407974,height_median=0.059539972645407974,height_range={0.059539972645407974,0.059539972645407974}]:0.11285028903661318,(#H2[&height_95%HPD={0.034771350355517344,0.034771350355517344},height_mean=0.034771350355517344,height_median=0.034771350355517344,height_range={0.034771350355517344,0.034771350355517344}]:0.10763810177414512,((M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.006968833074350089,N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.006968833074350089)S10[&height_95%HPD={0.006968833074350089,0.006968833074350089},height_mean=0.006968833074350089,height_median=0.006968833074350089,height_range={0.006968833074350089,0.006968833074350089}]:0.08028710948181786,((#H1[&height_95%HPD={0.03358447303437434,0.03358447303437434},height_mean=0.03358447303437434,height_median=0.03358447303437434,height_range={0.03358447303437434,0.03358447303437434}]:0.012514205423397484,E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.046098678457771824)S11[&height_95%HPD={0.046098678457771824,0.046098678457771824},height_mean=0.046098678457771824,height_median=0.046098678457771824,height_range={0.046098678457771824,0.046098678457771824}]:0.03296698062506892,A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.07906565908284074)S9[&height_95%HPD={0.07906565908284074,0.07906565908284074},height_mean=0.07906565908284074,height_median=0.07906565908284074,height_range={0.07906565908284074,0.07906565908284074}]:0.008190283473327203)S6[&height_95%HPD={0.08725594255616795,0.08725594255616795},height_mean=0.08725594255616795,height_median=0.08725594255616795,height_range={0.08725594255616795,0.08725594255616795}]:0.05515350957349452)S4[&height_95%HPD={0.14240945212966247,0.14240945212966247},height_mean=0.14240945212966247,height_median=0.14240945212966247,height_range={0.14240945212966247,0.14240945212966247}]:0.029980809552358684)S3[&height_95%HPD={0.17239026168202115,0.17239026168202115},height_mean=0.17239026168202115,height_median=0.17239026168202115,height_range={0.17239026168202115,0.17239026168202115}]:0.037269261084317495)S1[&height_95%HPD={0.20965952276633865,0.20965952276633865},height_mean=0.20965952276633865,height_median=0.20965952276633865,height_range={0.20965952276633865,0.20965952276633865}]:0.040340477233661354)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#6", "((((A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.07428991324506573,(((J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01251803332741927,I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01251803332741927)S14[&height_95%HPD={0.01251803332741927,0.01251803332741927},height_mean=0.01251803332741927,height_median=0.01251803332741927,height_range={0.01251803332741927,0.01251803332741927}]:0.03105118236441451,(D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03128520698616055,E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03128520698616055)S13[&height_95%HPD={0.03128520698616055,0.03128520698616055},height_mean=0.03128520698616055,height_median=0.03128520698616055,height_range={0.03128520698616055,0.03128520698616055}]:0.012284008705673227)S10[&height_95%HPD={0.04356921569183378,0.04356921569183378},height_mean=0.04356921569183378,height_median=0.04356921569183378,height_range={0.04356921569183378,0.04356921569183378}]:0.007221611986876042,((N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0047427194904210435,M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0047427194904210435)S17[&height_95%HPD={0.0047427194904210435,0.0047427194904210435},height_mean=0.0047427194904210435,height_median=0.0047427194904210435,height_range={0.0047427194904210435,0.0047427194904210435}]:0.012019815416266011,H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.016762534906687054)S11[&height_95%HPD={0.016762534906687054,0.016762534906687054},height_mean=0.016762534906687054,height_median=0.016762534906687054,height_range={0.016762534906687054,0.016762534906687054}]:0.03402829277202277)S7[&height_95%HPD={0.05079082767870982,0.05079082767870982},height_mean=0.05079082767870982,height_median=0.05079082767870982,height_range={0.05079082767870982,0.05079082767870982}]:0.02349908556635591)S4[&height_95%HPD={0.07428991324506573,0.07428991324506573},height_mean=0.07428991324506573,height_median=0.07428991324506573,height_range={0.07428991324506573,0.07428991324506573}]:0.044451341436778674,(((O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0042660216792135275,P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0042660216792135275)S16[&height_95%HPD={0.0042660216792135275,0.0042660216792135275},height_mean=0.0042660216792135275,height_median=0.0042660216792135275,height_range={0.0042660216792135275,0.0042660216792135275}]:0.019731281821523783,G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02399730350073731)S15[&height_95%HPD={0.02399730350073731,0.02399730350073731},height_mean=0.02399730350073731,height_median=0.02399730350073731,height_range={0.02399730350073731,0.02399730350073731}]:0.016654832201972558)#H2[&gamma=0.6251630957959864,gamma_95%HPD={0.6251630957959864,0.6251630957959864},gamma_mean=0.6251630957959864,gamma_median=0.6251630957959864,gamma_range={0.6251630957959864,0.6251630957959864},height_95%HPD={0.04065213570270987,0.04065213570270987},height_mean=0.04065213570270987,height_median=0.04065213570270987,height_range={0.04065213570270987,0.04065213570270987}]:0.07808911897913454)S2[&height_95%HPD={0.1187412546818444,0.1187412546818444},height_mean=0.1187412546818444,height_median=0.1187412546818444,height_range={0.1187412546818444,0.1187412546818444}]:0.04086417307620652,(((F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.029125365035769624)#H1[&gamma=0.7452128570146695,gamma_95%HPD={0.7452128570146695,0.7452128570146695},gamma_mean=0.7452128570146695,gamma_median=0.7452128570146695,gamma_range={0.7452128570146695,0.7452128570146695},height_95%HPD={0.029125365035769624,0.029125365035769624},height_mean=0.029125365035769624,height_median=0.029125365035769624,height_range={0.029125365035769624,0.029125365035769624}]:0.04275379024623671,((#H1[&height_95%HPD={0.029125365035769624,0.029125365035769624},height_mean=0.029125365035769624,height_median=0.029125365035769624,height_range={0.029125365035769624,0.029125365035769624}]:0.008123973864285039,C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03724933890005466)S12[&height_95%HPD={0.03724933890005466,0.03724933890005466},height_mean=0.03724933890005466,height_median=0.03724933890005466,height_range={0.03724933890005466,0.03724933890005466}]:0.00788299500226744,#H2[&height_95%HPD={0.04065213570270987,0.04065213570270987},height_mean=0.04065213570270987,height_median=0.04065213570270987,height_range={0.04065213570270987,0.04065213570270987}]:0.004480198199612234)S8[&height_95%HPD={0.0451323339023221,0.0451323339023221},height_mean=0.0451323339023221,height_median=0.0451323339023221,height_range={0.0451323339023221,0.0451323339023221}]:0.026746821379684232)S5[&height_95%HPD={0.07187915528200634,0.07187915528200634},height_mean=0.07187915528200634,height_median=0.07187915528200634,height_range={0.07187915528200634,0.07187915528200634}]:0.03612713802682141,((K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.008725207530716733,L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.008725207530716733)S9[&height_95%HPD={0.008725207530716733,0.008725207530716733},height_mean=0.008725207530716733,height_median=0.008725207530716733,height_range={0.008725207530716733,0.008725207530716733}]:0.056485553231029634,B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.06521076076174637)S6[&height_95%HPD={0.06521076076174637,0.06521076076174637},height_mean=0.06521076076174637,height_median=0.06521076076174637,height_range={0.06521076076174637,0.06521076076174637}]:0.04279553254708138)S3[&height_95%HPD={0.10800629330882774,0.10800629330882774},height_mean=0.10800629330882774,height_median=0.10800629330882774,height_range={0.10800629330882774,0.10800629330882774}]:0.05159913444922318)S1[&height_95%HPD={0.15960542775805092,0.15960542775805092},height_mean=0.15960542775805092,height_median=0.15960542775805092,height_range={0.15960542775805092,0.15960542775805092}]:0.09039457224194908)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#5", "((((M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.006511976626668359,N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.006511976626668359)S10[&height_95%HPD={0.006511976626668359,0.006511976626668359},height_mean=0.006511976626668359,height_median=0.006511976626668359,height_range={0.006511976626668359,0.006511976626668359}]:0.05796408854169238,(((G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.012716650222989873,H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.012716650222989873)S15[&height_95%HPD={0.012716650222989873,0.012716650222989873},height_mean=0.012716650222989873,height_median=0.012716650222989873,height_range={0.012716650222989873,0.012716650222989873}]:0.003987312570966739,E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01670396279395661)S13[&height_95%HPD={0.01670396279395661,0.01670396279395661},height_mean=0.01670396279395661,height_median=0.01670396279395661,height_range={0.01670396279395661,0.01670396279395661}]:0.01646783245455466,((I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.011328477641036921,J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.011328477641036921)S14[&height_95%HPD={0.011328477641036921,0.011328477641036921},height_mean=0.011328477641036921,height_median=0.011328477641036921,height_range={0.011328477641036921,0.011328477641036921}]:0.00818671545179847)#H1[&gamma=0.9209349514232135,gamma_95%HPD={0.9209349514232135,0.9209349514232135},gamma_mean=0.9209349514232135,gamma_median=0.9209349514232135,gamma_range={0.9209349514232135,0.9209349514232135},height_95%HPD={0.019515193092835392,0.019515193092835392},height_mean=0.019515193092835392,height_median=0.019515193092835392,height_range={0.019515193092835392,0.019515193092835392}]:0.013656602155675879)S11[&height_95%HPD={0.03317179524851127,0.03317179524851127},height_mean=0.03317179524851127,height_median=0.03317179524851127,height_range={0.03317179524851127,0.03317179524851127}]:0.03130426991984947)S3[&height_95%HPD={0.06447606516836074,0.06447606516836074},height_mean=0.06447606516836074,height_median=0.06447606516836074,height_range={0.06447606516836074,0.06447606516836074}]:0.027356680082464113,((((K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.00884644591736572,L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.00884644591736572)S16[&height_95%HPD={0.00884644591736572,0.00884644591736572},height_mean=0.00884644591736572,height_median=0.00884644591736572,height_range={0.00884644591736572,0.00884644591736572}]:0.005153904595242564,F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.014000350512608284)S8[&height_95%HPD={0.014000350512608284,0.014000350512608284},height_mean=0.014000350512608284,height_median=0.014000350512608284,height_range={0.014000350512608284,0.014000350512608284}]:0.05554490942588686,((D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01893055403227445,C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01893055403227445)S12[&height_95%HPD={0.01893055403227445,0.01893055403227445},height_mean=0.01893055403227445,height_median=0.01893055403227445,height_range={0.01893055403227445,0.01893055403227445}]:0.041376451479717746,#H1[&height_95%HPD={0.019515193092835392,0.019515193092835392},height_mean=0.019515193092835392,height_median=0.019515193092835392,height_range={0.019515193092835392,0.019515193092835392}]:0.0407918124191568)S7[&height_95%HPD={0.060307005511992196,0.060307005511992196},height_mean=0.060307005511992196,height_median=0.060307005511992196,height_range={0.060307005511992196,0.060307005511992196}]:0.009238254426502945)S5[&height_95%HPD={0.06954525993849514,0.06954525993849514},height_mean=0.06954525993849514,height_median=0.06954525993849514,height_range={0.06954525993849514,0.06954525993849514}]:0.014923291984040576,(A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.08040204700490594,((P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.005080017845952273,O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.005080017845952273)S9[&height_95%HPD={0.005080017845952273,0.005080017845952273},height_mean=0.005080017845952273,height_median=0.005080017845952273,height_range={0.005080017845952273,0.005080017845952273}]:0.06342963277725971,B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.06850965062321199)S6[&height_95%HPD={0.06850965062321199,0.06850965062321199},height_mean=0.06850965062321199,height_median=0.06850965062321199,height_range={0.06850965062321199,0.06850965062321199}]:0.011892396381693954)S4[&height_95%HPD={0.08040204700490594,0.08040204700490594},height_mean=0.08040204700490594,height_median=0.08040204700490594,height_range={0.08040204700490594,0.08040204700490594}]:0.004066504917629776)S2[&height_95%HPD={0.08446855192253572,0.08446855192253572},height_mean=0.08446855192253572,height_median=0.08446855192253572,height_range={0.08446855192253572,0.08446855192253572}]:0.007364193328289137)S1[&height_95%HPD={0.09183274525082485,0.09183274525082485},height_mean=0.09183274525082485,height_median=0.09183274525082485,height_range={0.09183274525082485,0.09183274525082485}]:0.15816725474917515)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#4", "(((((D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.01253296527202008,((K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.00331938251871558,L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.00331938251871558)S16[&height_95%HPD={0.00331938251871558,0.00331938251871558},height_mean=0.00331938251871558,height_median=0.00331938251871558,height_range={0.00331938251871558,0.00331938251871558}]:0.0063090087857836374,(H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.006288265367084406)#H1[&gamma=0.40145069761339935,gamma_95%HPD={0.40145069761339935,0.40145069761339935},gamma_mean=0.40145069761339935,gamma_median=0.40145069761339935,gamma_range={0.40145069761339935,0.40145069761339935},height_95%HPD={0.006288265367084406,0.006288265367084406},height_mean=0.006288265367084406,height_median=0.006288265367084406,height_range={0.006288265367084406,0.006288265367084406}]:0.003340125937414812)S14[&height_95%HPD={0.009628391304499218,0.009628391304499218},height_mean=0.009628391304499218,height_median=0.009628391304499218,height_range={0.009628391304499218,0.009628391304499218}]:0.0029045739675208626)S11[&height_95%HPD={0.01253296527202008,0.01253296527202008},height_mean=0.01253296527202008,height_median=0.01253296527202008,height_range={0.01253296527202008,0.01253296527202008}]:0.00746514322268621,(E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.011194052268471033,(M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.002244139591575628,N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.002244139591575628)S15[&height_95%HPD={0.002244139591575628,0.002244139591575628},height_mean=0.002244139591575628,height_median=0.002244139591575628,height_range={0.002244139591575628,0.002244139591575628}]:0.008949912676895405)S12[&height_95%HPD={0.011194052268471033,0.011194052268471033},height_mean=0.011194052268471033,height_median=0.011194052268471033,height_range={0.011194052268471033,0.011194052268471033}]:0.008804056226235257)S4[&height_95%HPD={0.01999810849470629,0.01999810849470629},height_mean=0.01999810849470629,height_median=0.01999810849470629,height_range={0.01999810849470629,0.01999810849470629}]:0.10950563141171984,(F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.008681184007997583,G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.008681184007997583)S5[&height_95%HPD={0.008681184007997583,0.008681184007997583},height_mean=0.008681184007997583,height_median=0.008681184007997583,height_range={0.008681184007997583,0.008681184007997583}]:0.12082255589842855)S2[&height_95%HPD={0.12950373990642614,0.12950373990642614},height_mean=0.12950373990642614,height_median=0.12950373990642614,height_range={0.12950373990642614,0.12950373990642614}]:0.012673908455837157,((((O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0021956568641774143,P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0021956568641774143)S13[&height_95%HPD={0.0021956568641774143,0.0021956568641774143},height_mean=0.0021956568641774143,height_median=0.0021956568641774143,height_range={0.0021956568641774143,0.0021956568641774143}]:0.011088951800291763,C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.013284608664469177)S8[&height_95%HPD={0.013284608664469177,0.013284608664469177},height_mean=0.013284608664469177,height_median=0.013284608664469177,height_range={0.013284608664469177,0.013284608664469177}]:0.023408878765174695,A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03669348742964387)S7[&height_95%HPD={0.03669348742964387,0.03669348742964387},height_mean=0.03669348742964387,height_median=0.03669348742964387,height_range={0.03669348742964387,0.03669348742964387}]:0.04415260703949869,(#H1[&height_95%HPD={0.006288265367084406,0.006288265367084406},height_mean=0.006288265367084406,height_median=0.006288265367084406,height_range={0.006288265367084406,0.006288265367084406}]:0.02275449688897188,(B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.024142143669066307,(J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.005812851389995682,I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.005812851389995682)S10[&height_95%HPD={0.005812851389995682,0.005812851389995682},height_mean=0.005812851389995682,height_median=0.005812851389995682,height_range={0.005812851389995682,0.005812851389995682}]:0.018329292279070625)S9[&height_95%HPD={0.024142143669066307,0.024142143669066307},height_mean=0.024142143669066307,height_median=0.024142143669066307,height_range={0.024142143669066307,0.024142143669066307}]:0.00490061858698998)S6[&height_95%HPD={0.029042762256056287,0.029042762256056287},height_mean=0.029042762256056287,height_median=0.029042762256056287,height_range={0.029042762256056287,0.029042762256056287}]:0.05180333221308628)S3[&height_95%HPD={0.08084609446914257,0.08084609446914257},height_mean=0.08084609446914257,height_median=0.08084609446914257,height_range={0.08084609446914257,0.08084609446914257}]:0.06133155389312073)S1[&height_95%HPD={0.1421776483622633,0.1421776483622633},height_mean=0.1421776483622633,height_median=0.1421776483622633,height_range={0.1421776483622633,0.1421776483622633}]:0.10782235163773671)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#3", "(((((C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.024023511698857775,B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.024023511698857775)S14[&height_95%HPD={0.024023511698857775,0.024023511698857775},height_mean=0.024023511698857775,height_median=0.024023511698857775,height_range={0.024023511698857775,0.024023511698857775}]:0.013511004194263793,(J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.012115584149811204)#H1[&gamma=0.6932078964591376,gamma_95%HPD={0.6932078964591376,0.6932078964591376},gamma_mean=0.6932078964591376,gamma_median=0.6932078964591376,gamma_range={0.6932078964591376,0.6932078964591376},height_95%HPD={0.012115584149811204,0.012115584149811204},height_mean=0.012115584149811204,height_median=0.012115584149811204,height_range={0.012115584149811204,0.012115584149811204}]:0.025418931743310363)S4[&height_95%HPD={0.03753451589312157,0.03753451589312157},height_mean=0.03753451589312157,height_median=0.03753451589312157,height_range={0.03753451589312157,0.03753451589312157}]:0.03182800772051861,(#H1[&height_95%HPD={0.012115584149811204,0.012115584149811204},height_mean=0.012115584149811204,height_median=0.012115584149811204,height_range={0.012115584149811204,0.012115584149811204}]:0.043952310746603024,(((G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.019617367950908182,H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.019617367950908182)S10[&height_95%HPD={0.019617367950908182,0.019617367950908182},height_mean=0.019617367950908182,height_median=0.019617367950908182,height_range={0.019617367950908182,0.019617367950908182}]:0.023011227064418327,(P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.001943210945015983,O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.001943210945015983)S11[&height_95%HPD={0.001943210945015983,0.001943210945015983},height_mean=0.001943210945015983,height_median=0.001943210945015983,height_range={0.001943210945015983,0.001943210945015983}]:0.040685384070310526)S7[&height_95%HPD={0.04262859501532651,0.04262859501532651},height_mean=0.04262859501532651,height_median=0.04262859501532651,height_range={0.04262859501532651,0.04262859501532651}]:0.005489914745591212,((D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02286216414812478,(I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.013609790096632685,(N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.00738694921682484,M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.00738694921682484)S16[&height_95%HPD={0.00738694921682484,0.00738694921682484},height_mean=0.00738694921682484,height_median=0.00738694921682484,height_range={0.00738694921682484,0.00738694921682484}]:0.006222840879807845)S15[&height_95%HPD={0.013609790096632685,0.013609790096632685},height_mean=0.013609790096632685,height_median=0.013609790096632685,height_range={0.013609790096632685,0.013609790096632685}]:0.009252374051492096)S9[&height_95%HPD={0.02286216414812478,0.02286216414812478},height_mean=0.02286216414812478,height_median=0.02286216414812478,height_range={0.02286216414812478,0.02286216414812478}]:0.022052153205261216,A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.044914317353385996)S8[&height_95%HPD={0.044914317353385996,0.044914317353385996},height_mean=0.044914317353385996,height_median=0.044914317353385996,height_range={0.044914317353385996,0.044914317353385996}]:0.003204192407531725)S6[&height_95%HPD={0.04811850976091772,0.04811850976091772},height_mean=0.04811850976091772,height_median=0.04811850976091772,height_range={0.04811850976091772,0.04811850976091772}]:0.007949385135496506)S5[&height_95%HPD={0.05606789489641423,0.05606789489641423},height_mean=0.05606789489641423,height_median=0.05606789489641423,height_range={0.05606789489641423,0.05606789489641423}]:0.013294628717225948)S2[&height_95%HPD={0.06936252361364018,0.06936252361364018},height_mean=0.06936252361364018,height_median=0.06936252361364018,height_range={0.06936252361364018,0.06936252361364018}]:0.07582700508537846,((F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.022825266704155633,E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.022825266704155633)S13[&height_95%HPD={0.022825266704155633,0.022825266704155633},height_mean=0.022825266704155633,height_median=0.022825266704155633,height_range={0.022825266704155633,0.022825266704155633}]:0.01897532804510585,(L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.008007054268468095,K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.008007054268468095)S12[&height_95%HPD={0.008007054268468095,0.008007054268468095},height_mean=0.008007054268468095,height_median=0.008007054268468095,height_range={0.008007054268468095,0.008007054268468095}]:0.03379354048079339)S3[&height_95%HPD={0.04180059474926148,0.04180059474926148},height_mean=0.04180059474926148,height_median=0.04180059474926148,height_range={0.04180059474926148,0.04180059474926148}]:0.10338893394975715)S1[&height_95%HPD={0.14518952869901863,0.14518952869901863},height_mean=0.14518952869901863,height_median=0.14518952869901863,height_range={0.14518952869901863,0.14518952869901863}]:0.10481047130098137)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#2", "((((D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.05771036965726989,E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.05771036965726989)S5[&height_95%HPD={0.05771036965726989,0.05771036965726989},height_mean=0.05771036965726989,height_median=0.05771036965726989,height_range={0.05771036965726989,0.05771036965726989}]:0.04626213739593327,(A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.09658144531234442,((P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.004862532945891923,O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.004862532945891923)S8[&height_95%HPD={0.004862532945891923,0.004862532945891923},height_mean=0.004862532945891923,height_median=0.004862532945891923,height_range={0.004862532945891923,0.004862532945891923}]:0.08446001396410843,((((N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.009663551955136324,M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.009663551955136324)S15[&height_95%HPD={0.009663551955136324,0.009663551955136324},height_mean=0.009663551955136324,height_median=0.009663551955136324,height_range={0.009663551955136324,0.009663551955136324}]:0.0018502956684125216,L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.011513847623548845)S12[&height_95%HPD={0.011513847623548845,0.011513847623548845},height_mean=0.011513847623548845,height_median=0.011513847623548845,height_range={0.011513847623548845,0.011513847623548845}]:0.047693403686148494,C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.05920725130969734)S9[&height_95%HPD={0.05920725130969734,0.05920725130969734},height_mean=0.05920725130969734,height_median=0.05920725130969734,height_range={0.05920725130969734,0.05920725130969734}]:0.015963259976522592,B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.07517051128621993)S7[&height_95%HPD={0.07517051128621993,0.07517051128621993},height_mean=0.07517051128621993,height_median=0.07517051128621993,height_range={0.07517051128621993,0.07517051128621993}]:0.014152035623780423)S6[&height_95%HPD={0.08932254691000036,0.08932254691000036},height_mean=0.08932254691000036,height_median=0.08932254691000036,height_range={0.08932254691000036,0.08932254691000036}]:0.007258898402344061)S4[&height_95%HPD={0.09658144531234442,0.09658144531234442},height_mean=0.09658144531234442,height_median=0.09658144531234442,height_range={0.09658144531234442,0.09658144531234442}]:0.00739106174085874)S3[&height_95%HPD={0.10397250705320316,0.10397250705320316},height_mean=0.10397250705320316,height_median=0.10397250705320316,height_range={0.10397250705320316,0.10397250705320316}]:0.03827810967344902,(((J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.013721401562071023,K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.013721401562071023)S14[&height_95%HPD={0.013721401562071023,0.013721401562071023},height_mean=0.013721401562071023,height_median=0.013721401562071023,height_range={0.013721401562071023,0.013721401562071023}]:0.006198634807491221,I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.019920036369562244)S11[&height_95%HPD={0.019920036369562244,0.019920036369562244},height_mean=0.019920036369562244,height_median=0.019920036369562244,height_range={0.019920036369562244,0.019920036369562244}]:0.05311371821237321,(F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.05255315671467581,(H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03388403764749137,G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03388403764749137)S13[&height_95%HPD={0.03388403764749137,0.03388403764749137},height_mean=0.03388403764749137,height_median=0.03388403764749137,height_range={0.03388403764749137,0.03388403764749137}]:0.018669119067184436)S10[&height_95%HPD={0.05255315671467581,0.05255315671467581},height_mean=0.05255315671467581,height_median=0.05255315671467581,height_range={0.05255315671467581,0.05255315671467581}]:0.020480597867259642)S2[&height_95%HPD={0.07303375458193545,0.07303375458193545},height_mean=0.07303375458193545,height_median=0.07303375458193545,height_range={0.07303375458193545,0.07303375458193545}]:0.06921686214471673)S1[&height_95%HPD={0.14225061672665218,0.14225061672665218},height_mean=0.14225061672665218,height_median=0.14225061672665218,height_range={0.14225061672665218,0.14225061672665218}]:0.10774938327334782)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#1", "((((((P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0011488989617972623,O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0011488989617972623)S14[&height_95%HPD={0.0011488989617972623,0.0011488989617972623},height_mean=0.0011488989617972623,height_median=0.0011488989617972623,height_range={0.0011488989617972623,0.0011488989617972623}]:0.04526966444241476,C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.04641856340421202)S10[&height_95%HPD={0.04641856340421202,0.04641856340421202},height_mean=0.04641856340421202,height_median=0.04641856340421202,height_range={0.04641856340421202,0.04641856340421202}]:0.006274845054608957,(G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.012070070388067788,H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.012070070388067788)S11[&height_95%HPD={0.012070070388067788,0.012070070388067788},height_mean=0.012070070388067788,height_median=0.012070070388067788,height_range={0.012070070388067788,0.012070070388067788}]:0.04062333807075319)S6[&height_95%HPD={0.05269340845882098,0.05269340845882098},height_mean=0.05269340845882098,height_median=0.05269340845882098,height_range={0.05269340845882098,0.05269340845882098}]:0.04387321421515933,(D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.04389120388794043,E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.04389120388794043)S7[&height_95%HPD={0.04389120388794043,0.04389120388794043},height_mean=0.04389120388794043,height_median=0.04389120388794043,height_range={0.04389120388794043,0.04389120388794043}]:0.052675418786039874)S2[&height_95%HPD={0.0965666226739803,0.0965666226739803},height_mean=0.0965666226739803,height_median=0.0965666226739803,height_range={0.0965666226739803,0.0965666226739803}]:0.04981313095694023,(((K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.006965254230692247,L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.006965254230692247)S8[&height_95%HPD={0.006965254230692247,0.006965254230692247},height_mean=0.006965254230692247,height_median=0.006965254230692247,height_range={0.006965254230692247,0.006965254230692247}]:0.0879084200542716,A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.09487367428496385)S4[&height_95%HPD={0.09487367428496385,0.09487367428496385},height_mean=0.09487367428496385,height_median=0.09487367428496385,height_range={0.09487367428496385,0.09487367428496385}]:0.004749778461906495,((((M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0037538986501611527,N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0037538986501611527)S15[&height_95%HPD={0.0037538986501611527,0.0037538986501611527},height_mean=0.0037538986501611527,height_median=0.0037538986501611527,height_range={0.0037538986501611527,0.0037538986501611527}]:0.020967100594211557,F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02472099924437271)S13[&height_95%HPD={0.02472099924437271,0.02472099924437271},height_mean=0.02472099924437271,height_median=0.02472099924437271,height_range={0.02472099924437271,0.02472099924437271}]:0.02773798756706003,(J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007843084355129148,I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.007843084355129148)S12[&height_95%HPD={0.007843084355129148,0.007843084355129148},height_mean=0.007843084355129148,height_median=0.007843084355129148,height_range={0.007843084355129148,0.007843084355129148}]:0.04461590245630359)S9[&height_95%HPD={0.05245898681143274,0.05245898681143274},height_mean=0.05245898681143274,height_median=0.05245898681143274,height_range={0.05245898681143274,0.05245898681143274}]:0.012924417269580823,B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.06538340408101356)S5[&height_95%HPD={0.06538340408101356,0.06538340408101356},height_mean=0.06538340408101356,height_median=0.06538340408101356,height_range={0.06538340408101356,0.06538340408101356}]:0.03424004866585678)S3[&height_95%HPD={0.09962345274687034,0.09962345274687034},height_mean=0.09962345274687034,height_median=0.09962345274687034,height_range={0.09962345274687034,0.09962345274687034}]:0.046756300884050195)S1[&height_95%HPD={0.14637975363092054,0.14637975363092054},height_mean=0.14637975363092054,height_median=0.14637975363092054,height_range={0.14637975363092054,0.14637975363092054}]:0.10362024636907946)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//        beastStrings.put("#0", "(((((H[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.042096758475518414,G[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.042096758475518414)S14[&height_95%HPD={0.042096758475518414,0.042096758475518414},height_mean=0.042096758475518414,height_median=0.042096758475518414,height_range={0.042096758475518414,0.042096758475518414}]:0.015240097936022712,E[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.057336856411541126)S5[&height_95%HPD={0.057336856411541126,0.057336856411541126},height_mean=0.057336856411541126,height_median=0.057336856411541126,height_range={0.057336856411541126,0.057336856411541126}]:0.11985275847322663,((((O[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:5.032484717042429E-4,P[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:5.032484717042429E-4)S11[&height_95%HPD={5.032484717042429E-4,5.032484717042429E-4},height_mean=5.032484717042429E-4,height_median=5.032484717042429E-4,height_range={5.032484717042429E-4,5.032484717042429E-4}]:0.07798132398258217,(M[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0055096873167276295,N[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0055096873167276295)S12[&height_95%HPD={0.0055096873167276295,0.0055096873167276295},height_mean=0.0055096873167276295,height_median=0.0055096873167276295,height_range={0.0055096873167276295,0.0055096873167276295}]:0.07297488513755879)S10[&height_95%HPD={0.07848457245428642,0.07848457245428642},height_mean=0.07848457245428642,height_median=0.07848457245428642,height_range={0.07848457245428642,0.07848457245428642}]:0.016782609032906792,((D[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.058417226965202046,C[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.058417226965202046)S13[&height_95%HPD={0.058417226965202046,0.058417226965202046},height_mean=0.058417226965202046,height_median=0.058417226965202046,height_range={0.058417226965202046,0.058417226965202046}]:0.0041346663766346214,B[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.06255189334183667)S9[&height_95%HPD={0.06255189334183667,0.06255189334183667},height_mean=0.06255189334183667,height_median=0.06255189334183667,height_range={0.06255189334183667,0.06255189334183667}]:0.03271528814535654)S6[&height_95%HPD={0.09526718148719321,0.09526718148719321},height_mean=0.09526718148719321,height_median=0.09526718148719321,height_range={0.09526718148719321,0.09526718148719321}]:0.02610314948977943,A[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.12137033097697264)S4[&height_95%HPD={0.12137033097697264,0.12137033097697264},height_mean=0.12137033097697264,height_median=0.12137033097697264,height_range={0.12137033097697264,0.12137033097697264}]:0.055819283907795114)S2[&height_95%HPD={0.17718961488476775,0.17718961488476775},height_mean=0.17718961488476775,height_median=0.17718961488476775,height_range={0.17718961488476775,0.17718961488476775}]:0.024638419781038146,((L[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.015712162697942678,K[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.015712162697942678)S8[&height_95%HPD={0.015712162697942678,0.015712162697942678},height_mean=0.015712162697942678,height_median=0.015712162697942678,height_range={0.015712162697942678,0.015712162697942678}]:0.0983393651689416,(F[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0473085236285527,(I[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.023739673167588682,J[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.023739673167588682)S15[&height_95%HPD={0.023739673167588682,0.023739673167588682},height_mean=0.023739673167588682,height_median=0.023739673167588682,height_range={0.023739673167588682,0.023739673167588682}]:0.02356885046096402)S7[&height_95%HPD={0.0473085236285527,0.0473085236285527},height_mean=0.0473085236285527,height_median=0.0473085236285527,height_range={0.0473085236285527,0.0473085236285527}]:0.06674300423833157)S3[&height_95%HPD={0.11405152786688427,0.11405152786688427},height_mean=0.11405152786688427,height_median=0.11405152786688427,height_range={0.11405152786688427,0.11405152786688427}]:0.08777650679892163)S1[&height_95%HPD={0.2018280346658059,0.2018280346658059},height_mean=0.2018280346658059,height_median=0.2018280346658059,height_range={0.2018280346658059,0.2018280346658059}]:0.0481719653341941)[&height_95%HPD={0.25,0.25},height_mean=0.25,height_median=0.25,height_range={0.25,0.25},topologySupport=1.0E-4];");
//
//        Map<String, String> paths = new HashMap<>();
//        paths.put("#14", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run17");
//        paths.put("#13", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run21");
//        paths.put("#12", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run22");
//        paths.put("#11", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run24");
//        paths.put("#10", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run25");
//        paths.put("#9", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run26");
//        paths.put("#8", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run27");
//        paths.put("#7", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run28");
//        paths.put("#6", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run29");
//        paths.put("#5", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run30");
//        paths.put("#4", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run31");
//        paths.put("#3", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run32");
//        paths.put("#2", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run33");
//        paths.put("#1", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run34");
//        paths.put("#0", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run35");
//
//        String tag = "#10";
//        CheckSimulationResult(tag, beastStrings.get(tag), paths.get(tag));
//    }

    public static void main(String[] args) {
        System.out.println(args.length);
        for(int i = 0 ; i < args.length ; i++) {
            System.out.println(args[i]);
        }

        //WalkThroughNetworks();
        //prepare(args);


        //measureCorrectness();
        check24(args);
    }
}
