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
        String root = args.length > 1 ? args[1] : ".";
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
                    out.println("SN_SEQ -cl 2000000 -bl 1000000 -sf 5000 -ee -sd 123456 -pl 2 -gtburnin -pre 10 ");
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
        String reduceRoot = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/Reduce/";
        int requiredIndex = 12;
        if(args.length > 1) {
            networkfilename = args[1];
            resultRoot = args[2].equals("null") ? null : args[2];
            reduceRoot = args[3].equals("null") ? null : args[3];
            if(!args[4].equals("def"))
                SuperNetwork3.eps = Double.parseDouble(args[4]);
            requiredIndex = Integer.parseInt(args[5]);
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

    public static void main(String[] args) {
        System.out.println(args.length);
        for(int i = 0 ; i < args.length ; i++) {
            System.out.println(args[i]);
        }

        if(args.length > 0) {
            if(args[0].equals("PrepareNetworkJson")) {
                WalkThroughNetworks();
            }
            else if(args[0].equals("PrepareNex")) {
                prepare(args);
            } else if(args[0].equals("Check24")) {
                check24(args);
            }
        }
    }
}
