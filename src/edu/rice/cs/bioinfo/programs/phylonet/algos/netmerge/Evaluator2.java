package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge;

/*
 * @ClassName:   Evaluator
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        1/2/24 10:28 PM
 */

import com.google.common.collect.Lists;
import com.opencsv.CSVReader;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.NetworkPrior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.MULTreeUtils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.Pipeline;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.NetworkPrior.getReticulationNodeDiameter;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.summarize.majorTree.removeReti;
import java.util.stream.Collectors;


public class Evaluator2 {
    /* Constructor */
    public Evaluator2() {

    }

    private static List<List<String>> parseClades(String cladesStr) {
        // Implement parsing logic here for cladesStr into List<List<String>>
        // This is a simplified version. You might need a more robust parser depending on the format.
        cladesStr = cladesStr.replaceAll("\\[\\[", "").replaceAll("\\]\\]", "");
        String[] groups = cladesStr.split("\\], \\[");
        List<List<String>> clades = new ArrayList<>();
        for (String group : groups) {
            clades.add(Arrays.asList(group.split(", ")));
        }
        return clades;
    }

    public static Map<String, List<List<List<String>>>> divided_reader(String path, int numcluster, String solver) throws Exception {
        try (CSVReader reader = new CSVReader(new FileReader(path))) {
            // Skip the header
            reader.readNext();
            String[] line;
            Map<String, List<List<List<String>>>> dataMap = new HashMap<>();
            while ((line = reader.readNext()) != null) {
                String scale = line[1];
                Integer index = Integer.parseInt(line[2]);
                if (!solver.equals(line[3]) || numcluster != Integer.parseInt(line[4])) {
                    continue;
                }
                String cladesStr = line[9];

                // Parse the clades as a list of lists.
                List<List<String>> clades = parseClades(cladesStr);

                // Build the nested map structure.
                dataMap.computeIfAbsent(scale, k -> new ArrayList<>())
                        .add(clades);
                return dataMap;

            }

        }catch (Exception e){
            System.out.println(e);
        }
        return null;
    }


    //read divided sets, numtip -> netid -> replica->clades
    public static Map<String, List<List<List<List<String>>>>> divided_reader_multinet(String path, String solver) throws Exception {
        try (CSVReader reader = new CSVReader(new FileReader(path))) {
            // Skip the header
            String[] line = reader.readNext();
            Map<String, List<List<List<List<String>>>>> dataMap = new HashMap<>();
            while ((line = reader.readNext()) != null) {
//                String[] arr = line.trim().split(",");
                String numtip = line[1];
                Integer netid = Integer.parseInt(line[2]);
                if (!solver.equals(line[6])) {
                    continue;
                }
                String cladesStr = line[12];

                // Parse the clades as a list of lists.
                List<List<String>> clades = parseClades(cladesStr);

                // Build the nested map structure.
                dataMap.putIfAbsent(numtip, new ArrayList<>());
                if (dataMap.get(numtip).size() < netid){
                    dataMap.get(numtip).add(new ArrayList<>());
                }
                dataMap.get(numtip).get(netid-1).add(clades);


            }
            return dataMap;

        }catch (Exception e){
            System.out.println(e);
        }
        return null;
    }



    public static void evaluateReplicaTruegtTruenet() throws Exception {
        Network net_original = Networks.readNetwork("(((((t1:1.8,#H1:1.4::0.5):2.9464,((t2:0.3052,t3:0.3052):3.658,((t4:0.4)#H1:1.0948::0.5,(t5:0.5832,t6:0.5832):0.9116):2.4684):0.7832):4.7008,(((t7:1.4428,t8:1.4428):4.4308,((t9:0.8492,t10:0.8492):4.3948,(((((t11:1.3768,t12:1.3768):0.4)#H3:0.638::0.6,(t13:0.4624,t14:0.4624):1.9524):0.5004,(t15:2.7952,(t16:2.0,#H3:0.2232::0.4):0.7952):0.12):0.2084,t17:3.1236):2.1204):0.6296):1.3796,((((t18:1.2424,t19:1.2424):0.898,t20:2.1404):2.1012,((t21:0.6644,t22:0.6644):2.2932,(t23:0.5328,t24:0.5328):2.4248):1.284):0.7936,(t25:2.3128,t26:2.3128):2.7224):2.218):2.194):1.1436,(t27:5.3544,((t28:0.6,#H2:0.2::0.4):0.8936,(t29:1.0344,(t30:0.4)#H2:0.6344::0.6):0.4592):3.8608):5.2364):89.4092,Z:100.0);");
        String dir_path = "/shared/zc36/merge/data/replica/";
        String[] scales = {"short", "medium", "long"};

        int numsites = 1000;
        String csvpath = "/shared/zc36/merge/data/replica/cluster_our_score_truegt.csv";
//        String csvpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/merge/data/replica/cluster_our_score_truegt.csv";
        Map<String, List<List<List<String>>>> divided_set = divided_reader(csvpath, 5, "nmf");

        for (String scale: scales){
            for (int i = 1; i <= 10; i ++){
                String distpath = dir_path+"/"+scale+"/"+i+"/"+numsites+ "/dist_matrix.txt";
                List<List<String>> subsets = divided_set.get(scale).get(i-1);
                List<Network> subnetworks = new ArrayList<>();
                for (List<String> subset: subsets){
                    for(int k = 0; k < subset.size(); k++) {
                        String updatedString = subset.get(k).replace("\'", "");
                        subset.set(k, updatedString);
                    }
                    Tuple<Network, Map<NetNode, NetNode>> subNet = SuperNetwork3.getSubNetwork(net_original, subset, true);
                    subnetworks.add(subNet.Item1);
                    NetNJMerge2 merger = new NetNJMerge2(distpath, subnetworks);
                    List<Network> netlist = merger.mergePairs();
                    for (Network net: netlist){
                        System.out.println(net.toString());
                        net.resetRoot("Z");
                        System.out.println(net);
                        System.out.println(scale+","+Networks.computeDistanceBetweenTwoNetworks(net, net_original));
                    }

                }

            }
        }

    }

    // read the matrix from a tree, where the distance between two species i and j are the entry of matrix[i][j]
    public static double[][] getMatrixFromTree(STITree tree, List<String> taxonlist) throws Exception {

        // get pair wise distance for every pair of species
        for (Object leaf: tree.getLeaves()){
            taxonlist.add(leaf.toString());
        }
        int leafcount = taxonlist.size();

        double [][] matrix = new double[leafcount][leafcount];
        for (int i = 0; i < leafcount; i++){
            for (int j = 0; j < leafcount; j++){
                matrix[i][j] = (double) (tree.getNodeDistanceFromDistance(taxonlist.get(i), taxonlist.get(j))).Item2;
//                matrix[i][j] = (double) (tree.getNodeDistance(taxonlist.get(i), taxonlist.get(j))).Item2;
            }
        }

        return matrix;
    }
    public static void compareOneReplica(Network truentwork, STITree tree, List<List<String>> subsets) throws Exception {
        List<Network> subnetworks = new ArrayList<>();
        for (List<String> subset: subsets) {
            for (int k = 0; k < subset.size(); k++) {
                String updatedString = subset.get(k).replace("\'", "");
                subset.set(k, updatedString);
            }
            Tuple<Network, Map<NetNode, NetNode>> subNet = SuperNetwork3.getSubNetwork(truentwork, subset, true);
            subnetworks.add(subNet.Item1);
        }

        List<String> taxonlist = new ArrayList<>();
        double [][] matrix = getMatrixFromTree(tree, taxonlist);
        NetNJMerge2 merger = new NetNJMerge2(matrix, subnetworks, taxonlist);

//        Network net = merger.mergePairs();
        List<Network> netlist = merger.mergePairs();
        for (Network net: netlist){
            System.out.println(net.toString());
            net.resetRoot("Z");
            System.out.println(net);
            System.out.println(Networks.computeDistanceBetweenTwoNetworks(net, truentwork));
        }



//        return new Tuple<Double, Network>(Networks.computeDistanceBetweenTwoNetworks(net, truentwork), net);
    }

    //read the path of truenetworks, each line is a network
    public static List<Network> readTrueNetworks(String path) throws Exception {
        List<Network> truenetworks = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(path));
        String line;
        while ((line = reader.readLine()) != null) {
            truenetworks.add(Networks.readNetwork(line));
        }
        return truenetworks;
    }

    public static List<STITree> readTrees(String path) throws Exception {
        List<STITree> treelist = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(path));
        String line;
        while ((line = reader.readLine()) != null) {
            treelist.add(new STITree(line));
        }
        return treelist;
    }


    //todo
    public static void evaluateMultiNetTruegtTruesubnet() throws Exception {
//        Network net_original = Networks.readNetwork("(((((t1:1.8,#H1:1.4::0.5):2.9464,((t2:0.3052,t3:0.3052):3.658,((t4:0.4)#H1:1.0948::0.5,(t5:0.5832,t6:0.5832):0.9116):2.4684):0.7832):4.7008,(((t7:1.4428,t8:1.4428):4.4308,((t9:0.8492,t10:0.8492):4.3948,(((((t11:1.3768,t12:1.3768):0.4)#H3:0.638::0.6,(t13:0.4624,t14:0.4624):1.9524):0.5004,(t15:2.7952,(t16:2.0,#H3:0.2232::0.4):0.7952):0.12):0.2084,t17:3.1236):2.1204):0.6296):1.3796,((((t18:1.2424,t19:1.2424):0.898,t20:2.1404):2.1012,((t21:0.6644,t22:0.6644):2.2932,(t23:0.5328,t24:0.5328):2.4248):1.284):0.7936,(t25:2.3128,t26:2.3128):2.7224):2.218):2.194):1.1436,(t27:5.3544,((t28:0.6,#H2:0.2::0.4):0.8936,(t29:1.0344,(t30:0.4)#H2:0.6344::0.6):0.4592):3.8608):5.2364):89.4092,Z:100.0);");
        String dir_path = "/shared/zc36/merge/data/multinet/";
        int numsites = 1000;

        String csvpath = "/shared/zc36/merge/data/multinet/cluster_score_filtered_relaxed.csv";

        for(int numtip = 20; numtip <= 20; numtip *= 2){
            String truenetworkpath = dir_path + numtip+"/network_newicks_outgroup.txt";
            String majortreepath = dir_path + numtip+"/major_tree.txt";
            List<Network> truenetworklist  = readTrueNetworks(truenetworkpath);
            List<STITree> treelist = readTrees(majortreepath);
            Map<String, List<List<List<List<String>>>>> subsets = divided_reader_multinet(csvpath, "nmf");
            System.out.println(subsets);
            for (int i = 1; i <= 10; i++){
                for (int j = 1; j <= 10; j++){
//                    String distpath = dir_path+"/net_"+i+"/"+j+"/"+numsites+ "/dist_matrix.txt";


//                    Tuple<Double, Network> result = compareOneReplica(truenetworklist.get(i), treelist.get(i), subsets);
//                    System.out.println(result.Item1);
//                    System.out.println(result.Item2);
                }

            }

        }

    }


    // break the network into subnetworks without breaking any cycles
    public static void getIdealSubSets(Network trunetwork){
        List<NetNode> nodes = Networks.postTraversal(trunetwork);
        for (NetNode node: nodes){
            if (node.isNetworkNode()){
                System.out.println(node.getName());
            }
        }

    }

    // get the leaves under reticulation node, and their parents. These are species that are ideally to be divided together
    public static List<Tuple3<Set<String>, Set<String>, Set<String>>> getNetNodeParentLeaves(Network network){
        List<Tuple3<Set<String>, Set<String>, Set<String>>> leavesTuples = new ArrayList<>();

        for (Object o: Networks.postTraversal(network)){
            NetNode node = (NetNode) o;
            if (node.isNetworkNode()){
                Set<String> leavesmom = new HashSet<>();
                Set<String> leavesdad = new HashSet<>();
                Set<String> leavesnode = new HashSet<>();

                Utils.getLeavesUnderNode(node, leavesnode);
                Iterator parentsIT = node.getParents().iterator();
                NetNode mom = (NetNode) parentsIT.next();
                NetNode dad = (NetNode) parentsIT.next();


                Utils.getLeavesUnderNode(mom, leavesmom);
                Utils.getLeavesUnderNode(dad, leavesdad);

                leavesmom.removeAll(leavesnode);
                leavesdad.removeAll(leavesnode);
                Tuple3<Set<String>, Set<String>, Set<String>> tuple3set = new Tuple3<>(leavesnode, leavesmom, leavesdad);
                leavesTuples.add(tuple3set);
            }
        }
        return leavesTuples;
    }


    // species related to a reticulation node are in the same cluster, other species are in a cluster
    public static List<Set<String>> divideSpecies(Network network){
        //species that are ideally to be divided together, transform to non-overlapping clusters
        List<Tuple3<Set<String>, Set<String>, Set<String>>> leavesTuples = getNetNodeParentLeaves(network);
        List<Set<String>> clusters = new ArrayList<>();
        for(Tuple3<Set<String>, Set<String>, Set<String>> tuple: leavesTuples){
            Set<String> tmp = new HashSet<>();
            tmp.addAll(tuple.Item1);
            tmp.addAll(tuple.Item2);
            tmp.addAll(tuple.Item3);

            boolean overlap = false;
            List<Integer> overlapset = new ArrayList<>();
            for (int i = 0; i < clusters.size(); i++){
                // check if tmp overlap with cluster
                Set<String> cluster = clusters.get(i);
                if (!Collections.disjoint(cluster, tmp)){
//                    cluster.addAll(tmp);
                    overlapset.add(i);
                    overlap = true;
//                    break;
                }
            }
//            for (Set<String> cluster: clusters){
//                // check if tmp overlap with cluster
//                if (!Collections.disjoint(cluster, tmp)){
//                    cluster.addAll(tmp);
//                    overlap = true;
//                    break;
//                }
//            }
            for (int i: overlapset){
                tmp.addAll(clusters.get(i));
            }
            Collections.sort(overlapset, Collections.reverseOrder());
            for (int i: overlapset){
                clusters.remove(i);
            }
//            if (!overlap){
            clusters.add(tmp);
//            }
//
        }

        //get the full set of species in the network
        Set<String> fullset = new HashSet<>();
        Utils.getLeavesUnderNode(network.getRoot(), fullset);


        // the rest of species are formed a cluster
        for (Set<String> cluster: clusters){
            fullset.removeAll(cluster);
        }
        clusters.add(fullset);

        return clusters;
    }


    // merge the subnetworks obtained from the truenetwork in the ideal subsets
    public static Tuple3<Integer, Network, Integer> mergeIdealSubsets(Network trunetwork, STITree tree, String outgroup) throws Exception{
        List<Set<String>> clusters = divideSpecies(trunetwork);
        List<Network> subnetworks = new ArrayList<>();
        int i = 0;
        for (Set<String> cluster: clusters){
            i += 1;
            System.out.println(cluster);
            if (cluster.size() < 3){
                System.out.println("cluster size < 3");
            return new Tuple3<>(100, null, 100);
            }
            //test
            Network subNet = SuperNetwork3.getSubNetwork(Networks.readNetwork(trunetwork.toString()), new ArrayList<>(cluster), true).Item1;

//            List<NetNode> netnodeList = Lists.newArrayList(subNet.getNetworkNodes());
//
//
//            removeReti(subNet, netnodeList, 0);
            subnetworks.add(subNet);

//            subnetworks.add(multreetuple.Item1);
            System.out.println(i+","+subNet.toString());
//            System.out.println(multreetuple.Item1.toString());
        }


        List<String> taxonlist = new ArrayList<>();
        double [][] matrix = getMatrixFromTree(tree, taxonlist);
//        for (int j = 0; j < taxonlist.size(); j++){
////            System.out.println(taxonlist.get(j));
//            System.out.print(j+" ");
//            for (int k = 0; k < taxonlist.size(); k++){
//                System.out.print(matrix[j][k]+" ");
//            }
//
//            System.out.println();
//        }
//        System.out.println(matrix);
        NetNJMerge2 merger = new NetNJMerge2(matrix, subnetworks, taxonlist);
        merger.setTrueNetwork(trunetwork);
        merger.setBackboneTree(Networks.readNetwork(tree.toString()));
        merger.setOutgroup(outgroup);
        int mindist = 100;
        int minsubDist = 0;
//        Network net = merger.mergePairs();
        List<Network> netlist = merger.mergePairs();
        System.out.println("true network:"+Networks.readNetwork(trunetwork.toString()));
        System.out.println("mergednetworks:"+netlist.size());
        Network bestnet = null;
        for(Network net: netlist){
            net.resetRoot(outgroup);
            System.out.println(net.toString());
            double distance = Networks.computeDistanceBetweenTwoNetworks(net, trunetwork);
            System.out.println("distance="+distance);
            if (distance < mindist){
                mindist = (int) distance;
                bestnet = net;
                minsubDist = 0;
                for (Network subnettrue: subnetworks){
                    List<String> leaves = new ArrayList<>();
                    for (Object o: subnettrue.getLeaves()){
                        NetNode node = (NetNode) o;
                        leaves.add(node.getName());
                    }
                    Network subnetinferred = (Network) SuperNetwork3.getSubNetwork(net, leaves, true).Item1;

                    minsubDist += Networks.computeDistanceBetweenTwoNetworks(subnetinferred, subnettrue);
                }
            }


        }
//        return false;
        return new Tuple3<>(mindist, bestnet, minsubDist);
//        System.out.println("merged network:");
    }

    public static Tuple3<Integer, Network, Integer> mergeIdealSubtrees(Network trunetwork, STITree tree, String outgroup) throws Exception{
        List<Set<String>> clusters = divideSpecies(trunetwork);
        List<Network> subnetworks = new ArrayList<>();
        int i = 0;
        for (Set<String> cluster: clusters){
            i += 1;
            System.out.println(cluster);
            //test
            Network subNet = SuperNetwork3.getSubNetwork(Networks.readNetwork(trunetwork.toString()), new ArrayList<>(cluster), true).Item1;

            List<NetNode> netnodeList = Lists.newArrayList(subNet.getNetworkNodes());


            removeReti(subNet, netnodeList, 0);
            subnetworks.add(subNet);

//            subnetworks.add(multreetuple.Item1);
            System.out.println(i+","+subNet.toString());
//            System.out.println(multreetuple.Item1.toString());
        }


        List<String> taxonlist = new ArrayList<>();
        double [][] matrix = getMatrixFromTree(tree, taxonlist);
//        for (int j = 0; j < taxonlist.size(); j++){
////            System.out.println(taxonlist.get(j));
//            System.out.print(j+" ");
//            for (int k = 0; k < taxonlist.size(); k++){
//                System.out.print(matrix[j][k]+" ");
//            }
//
//            System.out.println();
//        }
//        System.out.println(matrix);
        NetNJMerge2 merger = new NetNJMerge2(matrix, subnetworks, taxonlist);
        merger.setOutgroup(outgroup);
//        Network net = merger.mergePairs();
        List<Network> netlist = merger.mergePairs();
        System.out.println("true network:"+Networks.readNetwork(trunetwork.toString()));
        System.out.println("mergednetworks:"+netlist.size());
        int mindist = 100;
        Network bestnet = null;
        int minsubDist = 0;
        for(Network net: netlist) {
            net.resetRoot(outgroup);
            System.out.println(net.toString());
            double distance = Networks.computeDistanceBetweenTwoNetworks(net, Networks.readNetwork(tree.toString()));
            System.out.println("distance=" + distance);
            if (distance < mindist) {
                mindist = (int) distance;
                bestnet = net;
                minsubDist = 0;
                for (Network subnettrue : subnetworks) {
                    List<String> leaves = new ArrayList<>();
                    for (Object o : subnettrue.getLeaves()) {
                        NetNode node = (NetNode) o;
                        leaves.add(node.getName());
                    }
                    Network subnetinferred = (Network) SuperNetwork3.getSubNetwork(net, leaves, true).Item1;
                    minsubDist += Networks.computeDistanceBetweenTwoNetworks(subnetinferred, subnettrue);
                }
            }

            if (distance < mindist) {
                mindist = (int) distance;
                bestnet = net;
                minsubDist = 0;
                for (Network subnettrue : subnetworks) {
                    List<String> leaves = new ArrayList<>();
                    for (Object o : subnettrue.getLeaves()) {
                        NetNode node = (NetNode) o;
                        leaves.add(node.getName());
                    }
                    Network subnetinferred = (Network) SuperNetwork3.getSubNetwork(net, leaves, true).Item1;
                    minsubDist += Networks.computeDistanceBetweenTwoNetworks(subnetinferred, subnettrue);
                }
            }
        }
//        return false;
        return new Tuple3<>(mindist, bestnet, minsubDist);
//        System.out.println("merged network:");
    }


    public static boolean mergeIdealMulTree(Network trunetwork, STITree tree, String outgroup) throws Exception{
        List<Set<String>> clusters = divideSpecies(trunetwork);
        List<Network> subnetworks = new ArrayList<>();
        int i = 0;
        for (Set<String> cluster: clusters){
            i += 1;
            System.out.println(cluster);
            //test
            Tuple<Network, Map<NetNode, NetNode>> subNet = SuperNetwork3.getSubNetwork(Networks.readNetwork(trunetwork.toString()), new ArrayList<>(cluster), true);
            Tuple<Network, Map<String, Double>> multreetuple = MULTreeUtils.GetMULTree(subNet.Item1);
//            Tuple<Network, Map<NetNode, NetNode>> subNet = SuperNetwork3.getSubNetwork(trunetwork, new ArrayList<>(cluster), true);
            subnetworks.add(multreetuple.Item1);
            System.out.println(i+","+subNet.Item1.toString());
        }


        List<String> taxonlist = new ArrayList<>();
        double [][] matrix = getMatrixFromTree(tree, taxonlist);
//        for (int j = 0; j < taxonlist.size(); j++){
////            System.out.println(taxonlist.get(j));
//            System.out.print(j+" ");
//            for (int k = 0; k < taxonlist.size(); k++){
//                System.out.print(matrix[j][k]+" ");
//            }
//
//            System.out.println();
//        }
//        System.out.println(matrix);
        NetNJMerge2 merger = new NetNJMerge2(matrix, subnetworks, taxonlist);
        merger.setOutgroup(outgroup);
//        Network net = merger.mergePairs();
        List<Network> netlist = merger.mergePairs();
        System.out.println("true network:"+Networks.readNetwork(trunetwork.toString()));
        System.out.println("mergednetworks:"+netlist.size());
        for(Network net: netlist){
            net.resetRoot(outgroup);
            System.out.println(net.toString());
            double distance = Networks.computeDistanceBetweenTwoNetworks(net, Networks.readNetwork(trunetwork.toString()));
            System.out.println("distance="+distance);
            if (distance <0.001){
                return true;

            }
        }
        return false;
//        System.out.println("merged network:");
////        System.out.println(net.toString());
//        net.resetRoot(outgroup);
//        System.out.println(net);
//        System.out.println(Networks.readNetwork(trunetwork.toString()));
//        System.out.println("distance="+Networks.computeDistanceBetweenTwoNetworks(net, Networks.readNetwork(trunetwork.toString())));
    }


//    public static void matchIdealSubsets(Network trunetwork, STITree tree, String outgroup) throws Exception{
//        List<Set<String>> clusters = divideSpecies(trunetwork);
//        List<Network> subnetworks = new ArrayList<>();
//        int i = 0;
//        for (Set<String> cluster: clusters){
//            i += 1;
//            System.out.println(cluster);
//            //test
//            Tuple<Network, Map<NetNode, NetNode>> subNet = SuperNetwork3.getSubNetwork(Networks.readNetwork(trunetwork.toString()), new ArrayList<>(cluster), true);
//
//            if (subNet.Item1.getReticulationCount() >=1){
//                subnetworks.add(subNet.Item1);
//                System.out.println(i+","+subNet.Item1.toString());
//            }
//
//
//        }
//
//
//        List<String> taxonlist = new ArrayList<>();
////        double [][] matrix = getMatrixFromTree(tree, taxonlist);
//        Network backboneTree = Networks.readNetwork(tree.toString());
//
//        NetMatch matcher = new NetMatch(subnetworks, taxonlist, backboneTree);
////        merger.setOutgroup(outgroup);
////        Network net = merger.mergePairs();
//        List<Network> netlist = matcher.matchPairs();
//        System.out.println("true network:"+Networks.readNetwork(trunetwork.toString()));
//        System.out.println("mergednetworks:"+netlist.size());
//        for(Network net: netlist){
//            net.resetRoot(outgroup);
//            System.out.println(net.toString());
//            System.out.println("distance="+Networks.computeDistanceBetweenTwoNetworks(net, Networks.readNetwork(trunetwork.toString())));
//        }
//
//    }




    public static void evaluateMultiNetTruesubnetTrueDivide() throws Exception {
//        Network net_original = Networks.readNetwork("(((((t1:1.8,#H1:1.4::0.5):2.9464,((t2:0.3052,t3:0.3052):3.658,((t4:0.4)#H1:1.0948::0.5,(t5:0.5832,t6:0.5832):0.9116):2.4684):0.7832):4.7008,(((t7:1.4428,t8:1.4428):4.4308,((t9:0.8492,t10:0.8492):4.3948,(((((t11:1.3768,t12:1.3768):0.4)#H3:0.638::0.6,(t13:0.4624,t14:0.4624):1.9524):0.5004,(t15:2.7952,(t16:2.0,#H3:0.2232::0.4):0.7952):0.12):0.2084,t17:3.1236):2.1204):0.6296):1.3796,((((t18:1.2424,t19:1.2424):0.898,t20:2.1404):2.1012,((t21:0.6644,t22:0.6644):2.2932,(t23:0.5328,t24:0.5328):2.4248):1.284):0.7936,(t25:2.3128,t26:2.3128):2.7224):2.218):2.194):1.1436,(t27:5.3544,((t28:0.6,#H2:0.2::0.4):0.8936,(t29:1.0344,(t30:0.4)#H2:0.6344::0.6):0.4592):3.8608):5.2364):89.4092,Z:100.0);");
        String dir_path = "/shared/zc36/merge/data/multinet/";
        int numsites = 1000;
        int correctCnt = 0;
        int totalCnt = 0;
        boolean tree = false;
        int numnet = 1;
        List<Double> scorelist = new ArrayList<>();
        List<Double> nrdistList = new ArrayList<>();
        List<Double> closestDistList = new ArrayList<>();
        List<Double> diameterlist = new ArrayList<>();
        List<Double> numretiList = new ArrayList<>();
        int unabletoinfer = 0;
        List<String> filteredNetList = new ArrayList<>();
        List<String> filteredTreelist = new ArrayList<>();
//        List<Integer> parallelList = new ArrayList<>();
        for(int numtip = 20; numtip <= 20; numtip *= 2){

//            String majortreepath = dir_path + numtip+"/major_tree2000.txt";
//            String truenetworkpath = dir_path + numtip+ "/network_newicks_outgroup2000.txt";
//            String majortreepath = dir_path + numtip+"/major_tree1500.txt";
//            String truenetworkpath = dir_path + numtip+ "/network_newicks_outgroup1500.txt";

            String majortreepath = dir_path + numtip+"/major_tree1.txt";
            String truenetworkpath = dir_path + numtip+ "/network_newicks_outgroup1.txt";
//            String majortreepath = dir_path + numtip+"/major_tree.txt";
//            String truenetworkpath = dir_path + numtip+ "/network_newicks_outgroup.txt";


            List<STITree> treelist = readTrees(majortreepath);
            List<Network> truenetworklist  = readTrueNetworks(truenetworkpath);

            Utils._debug = false;

            for (int i = 0; i <= numnet-1; i++){
                System.out.println(i+"-----------------------");
                if (i == 11){
                    Utils._debug = true;
                }
//                else {
//                    continue;
//                }
                if(tree){
                    Tuple3<Integer, Network, Integer> tuple = mergeIdealSubtrees(truenetworklist.get(i), treelist.get(i), "Z");
                    int mindist = tuple.Item1;


                    Network inferrednet = tuple.Item2;
                    if (inferrednet!=null){
                        double nrdist =  mindist*1.0 /(numtip+truenetworklist.get(i).getReticulationCount()+inferrednet.getReticulationCount());
                        nrdistList.add(nrdist);
                        scorelist.add(tuple.Item3*1.0);
                        Utils.addInheritanceProb(inferrednet);
                        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(inferrednet, truenetworklist.get(i));
                        closestDistList.add(closest.Item3);
                        filteredNetList.add(truenetworklist.get(i).toString());
                        filteredTreelist.add(treelist.get(i).toString());

                    }
                    else {
                        nrdistList.add(100.0);
                        scorelist.add(-1*1.0);
                        closestDistList.add(-1.0);
                        unabletoinfer += 1;
                    }

                    if (mindist < 0.001){
                        correctCnt+=1;
                    }
                    else{
                        System.out.println("wrong net:"+numtip+" "+i);
                    }
                    totalCnt ++;
                }
                else{
                    Tuple3<Integer, Network, Integer> tuple = mergeIdealSubsets(truenetworklist.get(i), treelist.get(i), "Z");
                    int mindist = tuple.Item1;
                    Network net = truenetworklist.get(i);
                    int numReti = net.getReticulationCount();
                    double diameter = 0;
                    Map<NetNode, Double> distMap = NetworkPrior.getDiameterMap(net.toString());
                    for(NetNode key: distMap.keySet()) {
                        diameter += distMap.get(key);
//                        System.out.println(distMap.get(key));
                    }
                    System.out.println(numReti+","+diameter);
                    diameterlist.add(diameter);
                    numretiList.add(numReti*1.0);

                    Network inferrednet = tuple.Item2;
                    if (inferrednet!=null){
                        double nrdist =  mindist*1.0 /(numtip+truenetworklist.get(i).getReticulationCount()+inferrednet.getReticulationCount());
                        nrdistList.add(nrdist);
                        scorelist.add(tuple.Item3*1.0);
                        Utils.addInheritanceProb(inferrednet);
                        Tuple3<Network, Network, Double>  closest = Pipeline.CheckWithTrueNetwork(inferrednet, truenetworklist.get(i));
                        closestDistList.add(closest.Item3);
                        filteredNetList.add(truenetworklist.get(i).toString());
                        filteredTreelist.add(treelist.get(i).toString());

                    }
                    else {
                        nrdistList.add(100.0);
                        scorelist.add(-1.0);
                        closestDistList.add(-1.0);
                        unabletoinfer += 1;
                    }

                    if (mindist < 0.001){
                        correctCnt+=1;
                    }
                    else{
                        System.out.println("wrong net:"+numtip+" "+i);
                    }
                    totalCnt ++;

                }


            }
            String filePath = dir_path+numtip+"/network_newick_outgroup1.txt";
            Files.write(Paths.get(filePath), filteredNetList);
//            String treePath = dir_path+numtip+"/major_tree1.txt";
//            Files.write(Paths.get(treePath), filteredTreelist);

            System.out.println("correct:"+correctCnt);
            System.out.println("total:"+totalCnt);
            System.out.println(nrdistList);
            System.out.println(unabletoinfer);
            System.out.println(scorelist);
            System.out.println(closestDistList);
            System.out.println(diameterlist);
            System.out.println(numretiList);

            String suffix = "";
            if (tree){
                suffix = "_tree";
            }
            List<List<Double>> multipleDoubleLists = new ArrayList<>();
            multipleDoubleLists.add(nrdistList);
            multipleDoubleLists.add(scorelist);
            multipleDoubleLists.add(closestDistList);
            multipleDoubleLists.add(diameterlist);
            multipleDoubleLists.add(numretiList);

            String correctFilePath = dir_path+numtip+"/net_dist"+suffix+".txt";
            List<String> lines = multipleDoubleLists.stream()
                    .map(list -> list.stream()
                            .map(String::valueOf)
                            .collect(Collectors.joining(",")))
                    .collect(Collectors.toList());
            Files.write(Paths.get(correctFilePath), lines);


        }

    }


    public static void test(){
        Network inferred = Networks.readNetwork("((((t3,t20),((t2,t17),(((t11,t34),((t28,t16))#H100),(t38,(((t22,(t33)#H2),(t27,(t18,#H100))),#H2))))),((t26,(t39,t31)),(t19,(t29,(t10,t4))))),Z);");

        Network truenet = Networks.readNetwork("(Z:100.0,(((t19:0.6916443648,((t4:0.04115588418,t10:0.04115588418):0.1650568111,t29:0.2062126953):0.4854316696):2.461066606,(t26:0.6288571168,(t39:0.1486444659,t31:0.1486444659):0.4802126509):2.523853855):4.532398833,((t20:0.07990613658,t3:0.07990613658):6.728437711,((t17:0.5189419599,t2:0.5189419599):5.987387876,((t38:4.652012987,((t33:1.479894037)#H1:2.473201678::0.5535111291,((#H1:1.45797114::0.4464888709,t22:2.937865176):0.4574710767,(t27:3.188080979,(((t28:0.7278426847,t16:0.7278426847):0.2254566606)#H2:1.481069149::0.4658388651,t18:2.434368494):0.7537124848):0.2072552739):0.5577594616):0.698917272):0.7315794042,(#H2:0.1945802142::0.5341611349,(t34:0.42338387,t11:0.42338387):0.7244956895):4.235712831):1.122737445):0.302014012):0.8767659562):92.3148901963);");
        System.out.println(inferred.toString());
        System.out.println(Networks.computeDistanceBetweenTwoNetworks(inferred, truenet));
        System.out.println("distance="+Networks.computeDistanceBetweenTwoNetworks(inferred, truenet));
    }

    public static void evaluateMultiNets() throws Exception {
        String dir_path = "/shared/zc36/merge/data/multinet/";
        int numsites = 1000;

        int correctCnt = 0;
        int totalCnt = 0;
        for (int numtip = 20; numtip <= 160; numtip *= 2) {
            List<String> outputmultrees = new ArrayList<>();
            String truenetworkpath = dir_path + numtip + "/network_newicks_outgroup.txt";
//            String majortreepath = dir_path + numtip+"/major_tree.txt";
//            String outputfilepath = dir_path + numtip+"/multrees_true_subnet.txt";
            String inferrednetspath = dir_path + numtip + "/nets_true_subnet.txt";

            List<Network> truenetworklist = readTrueNetworks(truenetworkpath);
            List<Network> inferrednetworklist = readTrueNetworks(inferrednetspath);
            for (int i = 0; i <= 9; i++) {
//                System.out.println(i + "-----------------------");
                double distance = Networks.computeDistanceBetweenTwoNetworks(truenetworklist.get(i), inferrednetworklist.get(i));
                System.out.println(truenetworklist.get(i));
                System.out.println(inferrednetworklist.get(i));
                if (distance < 0.001) {

                    correctCnt += 1;
                } else {
                    System.out.println("wrong net:" + numtip + " " + i);
                }
                totalCnt++;

            }
        }
        System.out.println("correct:"+correctCnt);
        System.out.println("total:"+totalCnt);
    }
    public static void writeMajortrees(){
        String truetreepath = "/shared/zc36/merge/data/multinet/80/network_newicks_outgroup1500.txt";
        String majortreepath = "/shared/zc36/merge/data/multinet/80/major_tree1500.txt";
        try{
            List<Network> majortreelist = new ArrayList<>();
            List<Network> networklist = readTrueNetworks(truetreepath);
            for (Network net: networklist){
                List<NetNode> netnodeList = Lists.newArrayList(net.getNetworkNodes());
//                int retiIndex = 0;
                removeReti(net, netnodeList, 0);
                majortreelist.add(net);

            }
            BufferedWriter writer = new BufferedWriter(new FileWriter(majortreepath));
            for (Network net: majortreelist){
                writer.write(net.toString());
                writer.newLine();
            }
            writer.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }



    }

    public static Map<NetNode, Double> getDiameterMap(String netStr) {
        Network net = Networks.readNetwork(netStr);
        Map<NetNode, Double> distMap = new HashMap<NetNode, Double>();
        Map<NetNode, Tuple<NetNode, NetNode>> parentsMap = new HashMap<>();
        for(Object o : net.getNetworkNodes()) {
            NetNode node = (NetNode) o;
            List<NetNode> parents = IterableHelp.toList(node.getParents());
            if(parents.size() != 2) throw new IllegalArgumentException("Invalid network " + net.toString());
            distMap.put(node, 0.0);
            parentsMap.put(node, new Tuple<>(parents.get(0), parents.get(1)));
            for(NetNode par : parents) {
                distMap.put(node, distMap.get(node) + node.getParentDistance(par));
                if(node.getParentProbability(par) < 0.50) {
                    par.removeChild(node);
                }
            }
        }
        for(NetNode key: distMap.keySet()) {
            distMap.put(key, distMap.get(key) + getReticulationNodeDiameter(net, parentsMap.get(key)));
        }
        return distMap;
    }

    public static void judgeNets(){
        String dir_path = "/shared/zc36/merge/data/multinet/";
        int numsites = 1000;

        int correctCnt = 0;
        int totalCnt = 0;
        boolean tree = false;
        int numnet = 72;
        List<Integer> scorelist = new ArrayList<>();
        List<Double> nrdistList = new ArrayList<>();
        List<Double> closestDistList = new ArrayList<>();
        int unabletoinfer = 0;
        for(int numtip = 20; numtip <= 20; numtip *= 2) {

            String truenetworkpath = dir_path + numtip + "/network_newicks_outgroup1000.txt";
            try{
                List<Network> networklist = readTrueNetworks(truenetworkpath);
                for (Network net: networklist){
                    int numReti = net.getReticulationCount();
                    double diameter = 0;
                    Map<NetNode, Double> distMap = NetworkPrior.getDiameterMap(net.toString());
                    for(NetNode key: distMap.keySet()) {
                        diameter += distMap.get(key);
//                        System.out.println(distMap.get(key));
                    }
                    System.out.println(numReti+","+diameter);

                }

            }
            catch (Exception e){
                e.printStackTrace();
            }

        }



    }

    public static void testExample(){
        Network net1 = Networks.readNetwork("((((((t3)#H1:1.0::0.5,t24),t4),(((t18,t22),#H1:1.0::0.5),(t6:1.427090615,t9:1.427090615):5.113928554999999):0.3578702352):2.701196324,((t7:1.215284447,t14:1.215284447):8.190006208,((t8:0.8348553228,t20:0.8348553228):3.595846815,(((((t23:0.04817188842,t19:0.04817188842):0.1441017709,t27:0.1922736593):0.04564114961,t15:0.2379148089):1.733149218,t17:1.971064027):0.4534118845,t26:2.424475911):2.006226226):4.974588517):0.1947950738),Z:100.0);");
        Network net2 = Networks.readNetwork("((((((t18,t22),(t6:1.427090615,t9:1.427090615):5.113928554999999),(t3)#H1:1.0::0.5),((#H1:1.0::0.5,t24),t4)):2.701196324,((t7:1.215284447,t14:1.215284447):8.190006208,((t8:0.8348553228,t20:0.8348553228):3.595846815,(((((t23:0.04817188842,t19:0.04817188842):0.1441017709,t27:0.1922736593):0.04564114961,t15:0.2379148089):1.733149218,t17:1.971064027):0.4534118845,t26:2.424475911):2.006226226):4.974588517):0.1947950738),Z:100.0);");
        Network net3 = Networks.readNetwork("((((((t18,(t3)#H1:1.0::0.5),t22),(t6:1.427090615,t9:1.427090615):5.113928554999999),((#H1:1.0::0.5,t24),t4)):2.701196324,((t7:1.215284447,t14:1.215284447):8.190006208,((t8:0.8348553228,t20:0.8348553228):3.595846815,(((((t23:0.04817188842,t19:0.04817188842):0.1441017709,t27:0.1922736593):0.04564114961,t15:0.2379148089):1.733149218,t17:1.971064027):0.4534118845,t26:2.424475911):2.006226226):4.974588517):0.1947950738),Z:100.0);");
        Network netTrue = Networks.readNetwork("(((((t20:0.8348553228,t8:0.8348553228):3.595846815,((t17:1.971064027,(t15:0.2379148089,((t23:0.04817188842,t19:0.04817188842):0.1441017709,t27:0.1922736593):0.04564114961):1.733149218):0.4534118845,t26:2.424475911):2.006226226):4.974588517,(t14:1.215284447,t7:1.215284447):8.190006208):0.1947950738,((((t3:0.3130595985)#H1:2.277022167::0.5447802767,t24:2.590081765):2.451115297,t4:5.041197062):1.857692343,(((t6:1.427090615,t9:1.427090615):3.498980165,((t12:0.6709872192)#H2:2.739182801::0.6720805128,((t11:0.002951596819,t21:0.002951596819):0.6680356224,#H2:1.0::0.3279194872):2.739182801):1.51590076):1.61494839,((t18:0.3130595985,#H1:1.0::0.4552197233):3.036403515,t22:3.349463113):3.191556057):0.3578702352):2.701196324):90.39991427140001,Z:100.0);");
        List<NetNode> leafList = Lists.newArrayList(net1.getLeaves());
        Network subTrue = SuperNetwork3.getSubNetwork(netTrue, leafList.stream().map(NetNode::getName).collect(Collectors.toList()), true).Item1;
        System.out.println(Networks.computeDistanceBetweenTwoNetworks(net1, netTrue));
        System.out.println(Networks.computeDistanceBetweenTwoNetworks(net1, subTrue));
        System.out.println(Networks.computeDistanceBetweenTwoNetworks(net2, subTrue));
        System.out.println(Networks.computeDistanceBetweenTwoNetworks(net3, subTrue));

        Network net31 = Networks.readNetwork("(((((((t21,t11)I1,t12)I5,(t6:1.427090615,t9:1.427090615)I2)I6,((t18,(t3)#H1:1.0::0.5),t22)I4),((#H1:1.0::0.5,t24)I3,t4)):2.701196324,((t7:1.215284447,t14:1.215284447)I0:8.190006208,((t8:0.8348553228,t20:0.8348553228):3.595846815,(((((t23:0.04817188842,t19:0.04817188842):0.1441017709,t27:0.1922736593):0.04564114961,t15:0.2379148089):1.733149218,t17:1.971064027):0.4534118845,t26:2.424475911):2.006226226):4.974588517):0.1947950738),Z:100.0);");
        List<NetNode> leafList2 = Lists.newArrayList(net31.getLeaves());
        Network subTrue2 = SuperNetwork3.getSubNetwork(netTrue, leafList2.stream().map(NetNode::getName).collect(Collectors.toList()), true).Item1;
        Tuple3<Network, Network, Double> cloest = Pipeline.CheckWithTrueNetwork(net31, subTrue2);
        System.out.println(cloest.Item3);


    }


    public static void main(String[] args) throws Exception{
//        evaluateMultiNetTruesubnetMultreeTrueDivide();
//        writeMajortrees();
        evaluateMultiNetTruesubnetTrueDivide();
//        judgeNets();

//        evaluateMultiNets();
//        test();

//        testExample();
    }

}
