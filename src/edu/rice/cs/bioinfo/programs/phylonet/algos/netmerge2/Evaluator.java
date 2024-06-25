package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge2;

/*
 * @ClassName:   Evaluator
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        1/2/24 10:28 PM
 */

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.MULTreeUtils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge2.NetMatch;
import edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge2.NetNJMerge;
import edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge2.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

import com.opencsv.CSVReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;

public class Evaluator {
    /* Constructor */
    public Evaluator() {

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
                    NetNJMerge merger = new NetNJMerge(distpath, subnetworks);
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
        NetNJMerge merger = new NetNJMerge(matrix, subnetworks, taxonlist);

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

        for(int numtip = 20; numtip <= 160; numtip *= 2){
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
            for (Set<String> cluster: clusters){
                // check if tmp overlap with cluster
                if (!Collections.disjoint(cluster, tmp)){
                    cluster.addAll(tmp);
                    overlap = true;
                    break;
                }
            }
            if (!overlap){
                clusters.add(tmp);
            }
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
    public static Network mergeIdealSubsets(Network trunetwork, STITree tree, String outgroup) throws Exception{
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
//            subnetworks.add(subNet.Item1);
            subnetworks.add(multreetuple.Item1);
//            System.out.println(i+","+subNet.Item1.toString());
            System.out.println(multreetuple.Item1.toString());
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
        NetNJMerge merger = new NetNJMerge(matrix, subnetworks, taxonlist);
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
//            if (distance <0.001){
//                return true;
//
//            }
        }
//        return false;
        return netlist.get(0);
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
        NetNJMerge merger = new NetNJMerge(matrix, subnetworks, taxonlist);
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


    public static void matchIdealSubsets(Network trunetwork, STITree tree, String outgroup) throws Exception{
        List<Set<String>> clusters = divideSpecies(trunetwork);
        List<Network> subnetworks = new ArrayList<>();
        int i = 0;
        for (Set<String> cluster: clusters){
            i += 1;
            System.out.println(cluster);
            //test
            Tuple<Network, Map<NetNode, NetNode>> subNet = SuperNetwork3.getSubNetwork(Networks.readNetwork(trunetwork.toString()), new ArrayList<>(cluster), true);

            if (subNet.Item1.getReticulationCount() >=1){
                subnetworks.add(subNet.Item1);
                System.out.println(i+","+subNet.Item1.toString());
            }


        }


        List<String> taxonlist = new ArrayList<>();
//        double [][] matrix = getMatrixFromTree(tree, taxonlist);
        Network backboneTree = Networks.readNetwork(tree.toString());

        NetMatch matcher = new NetMatch(subnetworks, taxonlist, backboneTree);
//        merger.setOutgroup(outgroup);
//        Network net = merger.mergePairs();
        List<Network> netlist = matcher.matchPairs();
        System.out.println("true network:"+Networks.readNetwork(trunetwork.toString()));
        System.out.println("mergednetworks:"+netlist.size());
        for(Network net: netlist){
            net.resetRoot(outgroup);
            System.out.println(net.toString());
            System.out.println("distance="+Networks.computeDistanceBetweenTwoNetworks(net, Networks.readNetwork(trunetwork.toString())));
        }

    }




    public static void evaluateMultiNetTruesubnetTrueDivide() throws Exception {
//        Network net_original = Networks.readNetwork("(((((t1:1.8,#H1:1.4::0.5):2.9464,((t2:0.3052,t3:0.3052):3.658,((t4:0.4)#H1:1.0948::0.5,(t5:0.5832,t6:0.5832):0.9116):2.4684):0.7832):4.7008,(((t7:1.4428,t8:1.4428):4.4308,((t9:0.8492,t10:0.8492):4.3948,(((((t11:1.3768,t12:1.3768):0.4)#H3:0.638::0.6,(t13:0.4624,t14:0.4624):1.9524):0.5004,(t15:2.7952,(t16:2.0,#H3:0.2232::0.4):0.7952):0.12):0.2084,t17:3.1236):2.1204):0.6296):1.3796,((((t18:1.2424,t19:1.2424):0.898,t20:2.1404):2.1012,((t21:0.6644,t22:0.6644):2.2932,(t23:0.5328,t24:0.5328):2.4248):1.284):0.7936,(t25:2.3128,t26:2.3128):2.7224):2.218):2.194):1.1436,(t27:5.3544,((t28:0.6,#H2:0.2::0.4):0.8936,(t29:1.0344,(t30:0.4)#H2:0.6344::0.6):0.4592):3.8608):5.2364):89.4092,Z:100.0);");
        String dir_path = "/shared/zc36/merge/data/multinet/";
        int numsites = 1000;

        int correctCnt = 0;
        int totalCnt = 0;
        for(int numtip = 20; numtip <= 160; numtip *= 2){
            String truenetworkpath = dir_path + numtip+"/network_newicks_outgroup.txt";
            String majortreepath = dir_path + numtip+"/major_tree.txt";
            List<Network> truenetworklist  = readTrueNetworks(truenetworkpath);
            List<STITree> treelist = readTrees(majortreepath);

            Utils._debug = true;
            for (int i = 0; i <= 9; i++){
                System.out.println(i+"-----------------------");
//                if (i == 7){
//                    System.out.println("checkpoint");
//                    Utils._debug = true;
//                }
                if (truenetworklist.get(i).getReticulationCount() == 1){
                    continue;
                }
                Network mergedNet = mergeIdealSubsets(truenetworklist.get(i), treelist.get(i), "Z");
                if (Networks.computeDistanceBetweenTwoNetworks(mergedNet, truenetworklist.get(i)) < 0.001){
                    correctCnt+=1;
                }
                else{
                    System.out.println("wrong net:"+numtip+" "+i);
                }
                totalCnt ++;



            }

        }

        System.out.println("correct:"+correctCnt);
        System.out.println("total:"+totalCnt);

    }

    public static void evaluateMultiNetTruesubnetMultreeTrueDivide() throws Exception {
//        Network net_original = Networks.readNetwork("(((((t1:1.8,#H1:1.4::0.5):2.9464,((t2:0.3052,t3:0.3052):3.658,((t4:0.4)#H1:1.0948::0.5,(t5:0.5832,t6:0.5832):0.9116):2.4684):0.7832):4.7008,(((t7:1.4428,t8:1.4428):4.4308,((t9:0.8492,t10:0.8492):4.3948,(((((t11:1.3768,t12:1.3768):0.4)#H3:0.638::0.6,(t13:0.4624,t14:0.4624):1.9524):0.5004,(t15:2.7952,(t16:2.0,#H3:0.2232::0.4):0.7952):0.12):0.2084,t17:3.1236):2.1204):0.6296):1.3796,((((t18:1.2424,t19:1.2424):0.898,t20:2.1404):2.1012,((t21:0.6644,t22:0.6644):2.2932,(t23:0.5328,t24:0.5328):2.4248):1.284):0.7936,(t25:2.3128,t26:2.3128):2.7224):2.218):2.194):1.1436,(t27:5.3544,((t28:0.6,#H2:0.2::0.4):0.8936,(t29:1.0344,(t30:0.4)#H2:0.6344::0.6):0.4592):3.8608):5.2364):89.4092,Z:100.0);");
        String dir_path = "/shared/zc36/merge/data/multinet/";
        int numsites = 1000;

        int correctCnt = 0;
        int totalCnt = 0;
        for(int numtip = 20; numtip <= 160; numtip *= 2){
            List<String> outputmultrees = new ArrayList<>();
            String truenetworkpath = dir_path + numtip+"/network_newicks_outgroup.txt";
            String majortreepath = dir_path + numtip+"/major_tree.txt";
            String outputfilepath = dir_path + numtip+"/multrees_true_subnet.txt";
            List<Network> truenetworklist  = readTrueNetworks(truenetworkpath);
            List<STITree> treelist = readTrees(majortreepath);
            if (numtip == 40)
                Utils._debug = true;
            for (int i = 1; i <= 9; i++){
                System.out.println(i+"-----------------------");


                Network inferredNet = mergeIdealSubsets(truenetworklist.get(i), treelist.get(i), "Z");


                Network multree = MULTreeUtils.GetMULTree(truenetworklist.get(i)).Item1;
                outputmultrees.add(multree.toString());
//                if (Networks.computeDistanceBetweenTwoNetworks(multree, inferredNet) < 0.001){
//                    correctCnt+=1;
//                }
//                else{
//                    System.out.println("wrong net:"+numtip+" "+i);
//                }
//                totalCnt ++;

            }
            try {
                // Write the list of strings to the file
                Files.write(Paths.get(outputfilepath), outputmultrees);
                System.out.println("File written successfully.");
            } catch (IOException e) {
                // Handle any IO exceptions here
                e.printStackTrace();
            }


        }

        System.out.println("correct:"+correctCnt);
        System.out.println("total:"+totalCnt);

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

    public static void main(String[] args) throws Exception{
        evaluateMultiNetTruesubnetMultreeTrueDivide();
//        evaluateMultiNets();
//        test();

    }

}
