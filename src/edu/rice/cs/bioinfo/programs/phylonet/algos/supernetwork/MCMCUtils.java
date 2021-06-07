package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.summary.SummaryBL;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/29/18
 * Time: 9:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class MCMCUtils {
    // net needs to be ultrametric
    public static Map<String, Map<String, Double>> getPairwiseDistance(Network net) {
        SuperNetwork3.initNetHeights(net);

        List<String> taxa = new ArrayList<>();
        for(Object leafObj : net.getLeaves()) {
            NetNode leaf = (NetNode) leafObj;
            taxa.add(leaf.getName());
        }
        Collections.sort(taxa);

        Map<String, Map<String, Double>> dist = new HashMap<>();
        for(int i = 0 ; i < taxa.size() ; i++){
            dist.put(taxa.get(i), new HashMap<>());
            for (int j = 0; j < taxa.size(); j++) {
                double curdist = 0.0;
                if(i != j) {
                    curdist = Double.MAX_VALUE;
                    Set<NetNode> MCRAs = SuperNetwork3.getMCRAs(net, taxa.get(i), taxa.get(j));
                    for(NetNode node : MCRAs) {
                        double value = ((NetNodeInfo) node.getData()).getHeight();
                        curdist = Math.min(value, curdist);
                    }
                }
                dist.get(taxa.get(i)).put(taxa.get(j), curdist);
            }
        }
        return dist;
    }

    public static Map<String, Map<String, Double>> getPairwiseDistance(Tree tree, boolean useParentDistance) {
        List<String> taxa = new ArrayList<>();
        taxa.addAll(Arrays.asList(tree.getLeaves()) );
        Collections.sort(taxa);

        Map<String, Map<String, Double>> dist = new HashMap<>();
        for(int i = 0 ; i < taxa.size() ; i++){
            dist.put(taxa.get(i), new HashMap<>());
            for (int j = 0; j < taxa.size(); j++) {
                double curdist = 0.0;

                TNode n1 = tree.getNode(taxa.get(i));
                TNode n2 = tree.getNode(taxa.get(j));

                Map<TNode, Double> d1 = new HashMap<>();
                d1.put(n1, 0.0);
                TNode p1 = n1;
                double dd = 0;
                while(!p1.isRoot()) {
                    dd += useParentDistance ? p1.getParentDistance() : 1;
                    d1.put(p1.getParent(), dd);
                    p1 = p1.getParent();
                }

                TNode p2 = n2;
                double d2 = 0;

                while(p2 != null) {
                    if(d1.containsKey(p2)) {
                        curdist = d1.get(p2) + d2;
                        break;
                    }
                    d2 += useParentDistance ? p2.getParentDistance() : 1;
                    p2 = p2.getParent();
                }

                dist.get(taxa.get(i)).put(taxa.get(j), curdist);
            }
        }
        return dist;
    }

    static List<String> readLociList(String filename) {
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String s;
            int index = 0;
            List<String> result = new ArrayList<>();

            while((s = in.readLine()) != null) {
                s = s.trim();

                result.add(s);
                index++;
            }

            System.out.println("Read loci list: " + result.size());

            return result;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    public static void analyzeIQTREE(String resultFolder, String lociListFile) {
        List<String> lociList = readLociList(lociListFile);

        File path = new File(resultFolder);

        List<String> filenames = new ArrayList<>();

        File[] files = path.listFiles();

        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){
                if(files[i].getName().endsWith(".treefile")) {
                    String ss[] = files[i].getName().split("_");
                    if(lociList.contains(ss[1] + "_" + ss[2])) {
                        filenames.add(files[i].toString());
                    }
                }
            }
        }

        System.out.println(filenames.size());
        int count0 = 0;
        int count1 = 0;
        int count2 = 0;

        for(String filename : filenames) {
            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                s = in.readLine();
                Tree tree = Trees.readTree(s);

                STITree stitree = new STITree(tree);
                List<String> leaves = new ArrayList<>();
                leaves.add("SP03_indexing25_h0");
                leaves.add("SP03_indexing25_h1");
                leaves.add("SP09_indexing3_h0");
                leaves.add("SP09_indexing3_h1");
                //leaves.add("SP07_indexing9_h0");
                //leaves.add("SP07_indexing9_h1");
                leaves.add("SP10_indexing20_h0");
                leaves.add("SP10_indexing20_h1");

                stitree.rerootTreeAtEdge("SP07_indexing28_h0");
                stitree.constrainByLeaves(leaves);
                //System.out.println(stitree);

                tree.rerootTreeAtEdge("SP07_indexing28_h0");

                //System.out.println(tree.toNewick());
                Map<String, Map<String, Double>> pd = getPairwiseDistance(stitree, false);

                List<String> taxa = new ArrayList<>();
                taxa.addAll(Arrays.asList(tree.getLeaves()) );
                Collections.sort(taxa);

//                for(String t1 : taxa) {
//                    System.out.print("\t" + t1);
//                }
//                System.out.println();
//
//                for(String t1 : taxa) {
//                    System.out.print(t1);
//                    for(String t2 : taxa) {
//                        System.out.print("\t" + pd.get(t1).get(t2));
//                    }
//                    System.out.println();
//                }

//                double d0 = pd.get("SP03_indexing25_h0").get("SP09_indexing3_h0");
//                double d1 = pd.get("SP03_indexing25_h0").get("SP07_indexing9_h0");
//                double d2 = pd.get("SP09_indexing3_h0").get("SP07_indexing9_h0");



                double d0 = pd.get("SP03_indexing25_h0").get("SP09_indexing3_h0");
                double d1 = pd.get("SP03_indexing25_h0").get("SP10_indexing20_h0");
                double d2 = pd.get("SP09_indexing3_h0").get("SP10_indexing20_h0");
                System.out.println(d0 + " " + d1 + " " + d2);

                if(d0 < d1 && d0 < d2) {
                    count0++;
                } else if(d1 < d0 && d1 < d2){
                    count1++;
                } else if(d2 < d0 && d2 < d1) {
                    count2++;
                } else {
                    System.out.println(stitree.toNewick());
                }

                //System.out.println(tree.toNewick());


                //break;
            }catch (IOException e) {
                e.printStackTrace();
            }
        }

        System.out.println(count0 + " " + count1 + " " + count2);
    }

    public static void analyzeGeneTrees(String resultFolder) {
        File path = new File(resultFolder);

        List<String> filenames = new ArrayList<>();

        File[] files = path.listFiles();

        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){ //this line weeds out other directories/folders
                if(files[i].getName().startsWith("tree_") && files[i].getName().endsWith(".log")) {
                    //System.out.println(files[i]);
                    filenames.add(files[i].toString());
                }
            }
        }
        Collections.sort(filenames);
        Map<String, Tuple<String, Double>> file2bestgt = new HashMap<>();

        for(String filename : filenames) {
            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                int index = 0;
                String topSample = null;
                Map<Network, Integer> count = new HashMap<>();
                Map<Network, SummaryBL> summaryBL = new HashMap<>();
                boolean begin = true;
                int total = 0;
                double lastESS = 0;

                while((s = in.readLine()) != null) {

                    if (s.startsWith("(")) {
                        Network curSample = Networks.readNetwork(s);
                        total++;
                        boolean exist = false;
                        for (Network net : count.keySet()) {
                            if (Networks.hasTheSameTopology(net, curSample)) {
                                count.put(net, count.get(net) + 1);
                                summaryBL.get(net).addNetwork(curSample.clone());
                                exist = true;
                                break;
                            }
                        }
                        if (!exist) {
                            count.put(curSample, 1);
                            summaryBL.put(curSample, new SummaryBL(curSample.toString()));
                            summaryBL.get(curSample).addNetwork(curSample.clone());
                        }
                    } else {
                        break;
                    }
                    index++;
                }

                int totalValue = 0;
                int maxValue = 0;
                for(Network net : count.keySet()) {
                    if(maxValue < count.get(net)) {
                        topSample = Networks.getFullString(net);
                        maxValue = count.get(net);
                    }
                    totalValue += count.get(net);
                }


                in.close();


                if(topSample == null) {
                    continue;
                }

                file2bestgt.put(filename, new Tuple<>(topSample, 1.0 * maxValue / totalValue));


            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        int count = 0;
        Mean avg0 = new Mean();
        Mean avg1 = new Mean();
        for(String filename : file2bestgt.keySet()) {
            Tuple<String, Double> tuple = file2bestgt.get(filename);
            Map<String, Map<String, Double>> dist = getPairwiseDistance(Networks.readNetworkWithRootPop(tuple.Item1));
//            double d0 = dist.get("SP03_indexing25_h0").get("SP09_indexing3_h0");
//            double d1 = dist.get("SP03_indexing25_h0").get("SP07_indexing9_h0");
//            double d2 = dist.get("SP09_indexing3_h0").get("SP07_indexing9_h0");

            double d0 = dist.get("SP03_indexing25_h0").get("SP09_indexing3_h0");
            double d1 = dist.get("SP03_indexing25_h0").get("SP10_indexing20_h0");
            double d2 = dist.get("SP09_indexing3_h0").get("SP10_indexing20_h0");

            System.out.println(d0 + " " + d1 + " " + d2);

            if(d0 < d1 && d0 < d2) {
                count++;
                avg0.increment(tuple.Item2);
            } else {
                avg1.increment(tuple.Item2);
            }
            //System.out.println(tuple.Item1 + " " + tuple.Item2);
        }

        System.out.println(count);
        System.out.println(avg0.getResult() + " " + avg1.getResult());

    }

    public static void analyzeNetwork(String args[]) {
        //String resultFolder = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run47/Reti3_B/";
        String resultFolder = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run57/";

        int chainlen = 3000000;
        int burnin = 1000000;
        int sample_freq = 5000;

        int start = (int)(burnin / sample_freq) + 1;
        int end = (int)(chainlen / sample_freq);

        File path = new File(resultFolder);

        List<String> filenames = new ArrayList<>();

        File[] files = path.listFiles();
        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){ //this line weeds out other directories/folders
                if(files[i].toString().endsWith(".out")) {
                    //System.out.println(files[i]);
                    filenames.add(files[i].toString());
                }
            }
        }

        Collections.sort(filenames);

        Map<String, List<String>> file2samples = new HashMap<>();

        for(String filename : filenames) {
            file2samples.put(filename, new ArrayList<>());

            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                int index = 0;
                String topSample = null;
                Map<Network, Integer> count = new HashMap<>();
                boolean begin = false;
                int total = 0;
                double lastESS = 0.0;

                while((s = in.readLine()) != null) {

                    if(begin) {
                        if (s.startsWith("[")) {
                            Network curSample = Networks.readNetworkWithRootPop(s);
                            total++;
                            if(total >= start) {
                                file2samples.get(filename).add(s);
                                boolean exist = false;
                                for (Network net : count.keySet()) {
                                    if (Networks.hasTheSameTopology(net, curSample)) {
                                        count.put(net, count.get(net) + 1);
                                        exist = true;
                                        break;
                                    }
                                }
                                if (!exist) {
                                    count.put(curSample, 1);
                                }
                            }
                        } else {
                            String ss[] = s.split("\\s+");
                            if(ss.length == 7 && Character.isDigit(ss[2].charAt(0)) && ss[2].charAt(ss[2].length() - 1) == ';')
                                lastESS = Double.parseDouble(ss[2].substring(0, ss[2].length() - 1));

                        }
                    } else {
                        if(s.contains("Logger")) {
                            begin = true;
                        }
                    }
                    index++;
                }

                int totalValue = 0;
                int maxValue = 0;
                for(Network net : count.keySet()) {
                    if(maxValue < count.get(net)) {
                        topSample = Networks.getFullString(net);
                        maxValue = count.get(net);
                    }
                    totalValue += count.get(net);
                }


                in.close();



                if(topSample == null) {
                    continue;
                }

                Network<Object> topNet = Networks.readNetworkWithRootPop(topSample);
                //if(topNet.findNode("Carlia_vivax") != null && topNet.findNode("Carlia_amax") != null) {
                //if(topNet.findNode("Carlia_rhomboidalis") != null && topNet.findNode("Carlia_longipes") != null) {
                if(topNet.findNode("Lampropholis_coggeri") != null && topNet.findNode("Lampropholis_guichenoti") != null) {
                    for(NetNode leaf : topNet.getLeaves()) {
                        System.out.print(leaf.getName() + " ");
                    }
                    System.out.println();
                    System.out.println(filename);
                    System.out.println(Networks.getDendroscopeCompatibleString(Networks.readNetworkWithRootPop(topSample)) );
                    System.out.println(Networks.getFullString(Networks.readNetworkWithRootPop(topSample)) );

                }

//                if(lastESS < 20 || 1.0 * maxValue / totalValue < 0.5) {
//                    file2samples.remove(filename);
//                    continue;
//                }
//
//                List<String> toRemove = new ArrayList<>();
//                for(String sample : file2samples.get(filename)) {
//                    if(!Networks.hasTheSameTopology(Networks.readNetworkWithRootPop(sample), Networks.readNetworkWithRootPop(topSample))) {
//                        toRemove.add(sample);
//                    }
//                }
//                file2samples.get(filename).removeAll(toRemove);

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public static void main(String args[]) {
        //analyzeGeneTrees("/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run50/6443091_8");
        //analyzeGeneTrees("/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/debug/run_9_1");

        analyzeNetwork(args);
        //analyzeIQTREE("/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/lizard/data/iqtree", "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/lizard/data/lists/locilist_100.txt");
    }
}
