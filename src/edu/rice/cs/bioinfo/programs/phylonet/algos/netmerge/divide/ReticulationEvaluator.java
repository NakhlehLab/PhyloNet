package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.divide;
/*
 * @ClassName:   Reticulation_Evaluator
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        11/10/23 4:13 PM
 */

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.summarize.majorTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.*;

import java.util.*;
import java.util.stream.Collectors;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.NJMergeTopology.getMRCA;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.summarize.tripartition2.getLeavesUnderNode;

public class ReticulationEvaluator {
    public Network trueNetwork;
    /* Constructor */
    public ReticulationEvaluator(String trunetString) {
        trueNetwork = Networks.readNetwork(trunetString);
    }


    public static List<Tuple3<Set<String>, Set<String>, Set<String>>> getRetiSets(Network net){
        List<Tuple3<Set<String>, Set<String>, Set<String>>> retiLeafSets = new ArrayList<>();
        List<NetNode> node2process = new ArrayList<>();
        int reticount = 0;
        for (Object o: net.getNetworkNodes()){
            reticount += 1;
            NetNode reticulation = (NetNode) o;
            Iterator<NetNode> parents = reticulation.getParents().iterator();
            NetNode mom = parents.next();
            NetNode dad = parents.next();
            Set<String> leaves3 = getLeavesUnderNode(reticulation);
            Set<String> leaves1 = getLeavesUnderNode(mom);
            Set<String> leaves2 = getLeavesUnderNode(dad);
            leaves1.removeAll(leaves3);
            leaves2.removeAll(leaves3);

            System.out.println(leaves1.size()+","+leaves2.size()+","+leaves3.size());
            // check if leaves size are more than 10
            if(leaves1.size() > 10 || leaves2.size() > 10){

//                System.out.println("large set");

                node2process.add(reticulation);
            }
//            System.out.println(leaves1);
//            System.out.println(leaves2);
            retiLeafSets.add(new Tuple3<>(leaves1, leaves2, leaves3));

        }
        //todo
        if  (reticount == node2process.size()){

            System.out.println("to replace");
        }
        else if (0 < node2process.size()){
            for (NetNode reticulation: node2process) {
                Iterator<NetNode> parents = reticulation.getParents().iterator();
                NetNode mom = parents.next();
                NetNode dad = parents.next();
                // remove non major parent
                if (reticulation.getParentProbability(mom) > reticulation.getParentProbability(dad) && !dad.isRoot()){

                    Utils.removeReticulation(dad, reticulation);
                }
                else{
                    Utils.removeReticulation(mom, reticulation);
                }
                Networks.removeBinaryNodes(net);
            }
            System.out.println(net.toString());
            System.out.println("todo: remove reticulation");
        }
        return retiLeafSets;
    }

    public static List<List<Tuple3<Set<String>, Set<String>, Set<String>>>> readNewicks(String filepath) {
        List<List<Tuple3<Set<String>, Set<String>, Set<String>>>> retiLeafSetsList = new ArrayList<>();
        System.out.println(filepath);
        int i = 0;
        try (BufferedReader reader = new BufferedReader(new FileReader(filepath))) {
            String line;

            while ((line = reader.readLine()) != null) {
                i += 1;
                System.out.println("net:"+i);
                Network network = Networks.readNetwork(line.trim());
//                System.out.println(network.toString());
                List<Tuple3<Set<String>, Set<String>, Set<String>>> retiLeafSets = getRetiSets(network);
                retiLeafSetsList.add(retiLeafSets);

//                System.out.println(line);
            }
        } catch (IOException e) {
            System.err.println("An error occurred while reading the file: " + e.getMessage());
        }
        return retiLeafSetsList;
    }

    public static List<String> readNewicks2MajorTrees(String filepath) {
        List<String> majortreelist = new ArrayList<>();
        System.out.println(filepath);
        int i = 0;
        try (BufferedReader reader = new BufferedReader(new FileReader(filepath))) {
            String line;

            while ((line = reader.readLine()) != null) {
                i += 1;
                System.out.println("net:"+i);
                Network network = Networks.readNetwork(line.trim());
                Map<Network, Double> networkmap = new HashMap<>();
                networkmap.put(network, 1.0);
                Map<Tree, Double> majortreemap = majorTree.summarize(networkmap);
                majortreelist.add(majortreemap.keySet().iterator().next().toString());

            }
        } catch (IOException e) {
            System.err.println("An error occurred while reading the file: " + e.getMessage());
        }
        return majortreelist;
    }
    private static String setToString(Set<String> set) {
        return set.stream().collect(Collectors.joining(","));
    }

    public static List<String> writeRetiLeafSetsListToFile(List<List<Tuple3<Set<String>, Set<String>, Set<String>>>> retiLeafSetsList, String outputprefix) {
        int fileIndex = 1;
        List<String> largetsetfilename = new ArrayList<>();

        for (List<Tuple3<Set<String>, Set<String>, Set<String>>> list : retiLeafSetsList) {
            String fileName = outputprefix + fileIndex+"/reti_leaf_sets3.txt";
//            System.out.println(fileName);
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(fileName))) {
                for (Tuple3<Set<String>, Set<String>,  Set<String>> tuple : list) {
                    String line = setToString(tuple.Item1) + ";" + setToString(tuple.Item2)+";"+setToString(tuple.Item3) ;;
                    writer.write(line);
                    writer.newLine();
//                    if(tuple.Item1.size() > 10 || tuple.Item2.size() > 10){
//                        System.out.println("large set");
//                        if (!largetsetfilename.contains(fileName)){
//                            largetsetfilename.add(fileName);
//                        }
//
//                    }
//                    System.out.println(tuple.Item1.size()+","+tuple.Item2.size());
//                    System.out.println(line);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            fileIndex++;
        }
        return largetsetfilename;
    }

    public static void readAllNets(String directory){
        List<String> checkfilenames = new ArrayList<>();

        for (int numtips = 20; numtips <= 160; numtips*=2) {
            String filepath = directory + "/" + numtips + "/network_newicks_outgroup.txt";
            String outputprefix = directory + "/" + numtips+"/net_";
//            System.out.println(outputprefix);
            List<List<Tuple3<Set<String>, Set<String>, Set<String>>>> retiLeafSetsList = readNewicks(filepath);
            List<String> largesetFileName = writeRetiLeafSetsListToFile(retiLeafSetsList, outputprefix);
            checkfilenames.addAll(largesetFileName);


        }
//            readNewicks(filepath);
//        for (String filename: checkfilenames){
//            System.out.println(filename.substring(filename.indexOf("multinet") + 10, filename.indexOf("reti_leaf_sets")));
//        }
//        System.out.println(checkfilenames);
        
    }

    /** shits below... might be useful later
     * Checks if there are reticulations between any two of the given clades in a network.
     *
     * @param clades A list of sets of strings, each set representing a clade.
     * @param net The phylogenetic network to check for reticulations.
     * @return true if there are reticulations between any two clades, false otherwise.
     */
//    public static boolean isRetiBetweenClades(List<Set<String>> clades, Network net) {
//        // Implement the logic to check for reticulations between clades
//        // This will likely involve traversing the network and checking the relationships
//        // between nodes corresponding to the clades.
//        for(Set<String> clade: clades){
//            NetNode mrca = getMRCA(net, clade);
//            System.out.println(mrca.getName());
//            NetNode parent = (NetNode) mrca.getParents().iterator().next();
//            parent.removeChild(mrca);
//            System.out.println(net.toString());
//        }
//        System.out.println(net.getReticulationCount());
//        System.out.println(net.toString());
//
//
//
//        // The exact implementation depends on the details of how networks and nodes
//        // are represented in PhyloNet, and the methods available for examining them.
//
//        // Return true if a reticulation is found, false otherwise
//        return false;
//    }
//    public static void test(){
//        Network net = Networks.readNetwork("(((((t1:4.5,#H1:3.5::0.5):7.366,((t2:0.763,t3:0.763):9.145,((t4:1.0)#H1:2.737::0.5,(t5:1.458,t6:1.458):2.279):6.171):1.958):11.752,(((t7:3.607,t8:3.607):11.077,((t9:2.123,t10:2.123):10.987,(((((t11:3.442,t12:3.442):1.0)#H3:1.595::0.6,(t13:1.156,t14:1.156):4.881):1.251,(t15:6.988,(t16:5.0,#H3:0.558::0.4):1.988):0.3):0.521,t17:7.809):5.301):1.574):3.449,((((t18:3.106,t19:3.106):2.245,t20:5.351):5.253,((t21:1.661,t22:1.661):5.733,(t23:1.332,t24:1.332):6.062):3.21):1.984,(t25:5.782,t26:5.782):6.806):5.545):5.485):2.859,(t27:13.386,((t28:1.5,#H2:0.5::0.4):2.234,(t29:2.586,(t30:1.0)#H2:1.586::0.6):1.148):9.652):13.091):73.523,Z:100.0);");
//        Networks.autoLabelNodes(net);
//        List<Set<String>> clades = new ArrayList<>();
//        clades.add(new HashSet<>(Arrays.asList("t1", "t2", "t3", "t4", "t5", "t6")));
//        clades.add(new HashSet<>(Arrays.asList("t7", "t8", "t9", "t10", "t11", "t12")));
//        clades.add(new HashSet<>(Arrays.asList("t13", "t14", "t15", "t16", "t17", "t18")));
//        clades.add(new HashSet<>(Arrays.asList("t19", "t20", "t21", "t22", "t23", "t24")));
//        isRetiBetweenClades(clades, net);
//    }

    public static void getMajorTreeAllNets(String directory){
        List<String> checkfilenames = new ArrayList<>();

        for (int numtips = 20; numtips <= 160; numtips*=2) {
            String outputfile = directory + numtips + "/major_tree.txt";
            String filepath = directory + "/" + numtips + "/network_newicks_outgroup.txt";
//            System.out.println(outputprefix);
            List<String> majortreeList = readNewicks2MajorTrees(filepath);
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputfile))) {
                for (String treestr: majortreeList) {
                    writer.write(treestr);
                    writer.newLine();
//                    if(tuple.Item1.size() > 10 || tuple.Item2.size() > 10){
//                        System.out.println("large set");
//                        if (!largetsetfilename.contains(fileName)){
//                            largetsetfilename.add(fileName);
//                        }
//
//                    }
//                    System.out.println(tuple.Item1.size()+","+tuple.Item2.size());
//                    System.out.println(line);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            for (String tree : majortreeList){
                System.out.println(tree);

            }
        }

    }


    public static void main(String[] args) {
//        Network net = Networks.readNetwork("(((((t1:4.5,#H1:3.5::0.5):7.366,((t2:0.763,t3:0.763):9.145,((t4:1.0)#H1:2.737::0.5,(t5:1.458,t6:1.458):2.279):6.171):1.958):11.752,(((t7:3.607,t8:3.607):11.077,((t9:2.123,t10:2.123):10.987,(((((t11:3.442,t12:3.442):1.0)#H3:1.595::0.6,(t13:1.156,t14:1.156):4.881):1.251,(t15:6.988,(t16:5.0,#H3:0.558::0.4):1.988):0.3):0.521,t17:7.809):5.301):1.574):3.449,((((t18:3.106,t19:3.106):2.245,t20:5.351):5.253,((t21:1.661,t22:1.661):5.733,(t23:1.332,t24:1.332):6.062):3.21):1.984,(t25:5.782,t26:5.782):6.806):5.545):5.485):2.859,(t27:13.386,((t28:1.5,#H2:0.5::0.4):2.234,(t29:2.586,(t30:1.0)#H2:1.586::0.6):1.148):9.652):13.091):73.523,Z:100.0);");
//        getRetiSets(net);
//        String directory = "/home/zc36/merge/data/multinet/";
        String directory = "/Users/zhen/Desktop/Zhen/research/phylogenetics/merge/data/multinet/";
        readAllNets(directory);
//        getMajorTreeAllNets(directory);
    }
}
