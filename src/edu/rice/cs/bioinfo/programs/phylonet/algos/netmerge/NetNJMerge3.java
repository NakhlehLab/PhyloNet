package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge;
/*
 * @ClassName:   NetMerge
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        9/12/23 11:34 PM
 */

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils._debug;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils.addInheritanceProb;


public class NetNJMerge3 {
    //    private List<Network> _subnetworks = new ArrayList<>();
    private List<List<Network>> _subnetworks = new ArrayList<>();
    private int _num_taxa = 0;
    private double[][] _matrix = null;
    private List<String> _taxonList = new ArrayList<>();
    private List<Set<String>> _leaves = null;
    private String _outgroup = "";
    private double _scale = 1;

    /* Constructor */
    public NetNJMerge3(String distPath, String netPath) {
        readMatrix(distPath);
        readSubNetworks(netPath);

    }
    public void setOutgroup(String outgroup){
        this._outgroup = outgroup;
    }
    public NetNJMerge3(String distPath, List<Network> subnetlist){
        readMatrix(distPath);
        for (Network subnet: subnetlist){
            List<Network> tmp = new ArrayList<>();
            tmp.add(subnet);
            _subnetworks.add(tmp);
        }
    }
    public NetNJMerge3(double[][] matrix, List<Network> subnetlist, List<String> taxonList){
        _matrix = matrix;
        for (Network subnet: subnetlist){
            List<Network> tmp = new ArrayList<>();
            tmp.add(subnet);
            _subnetworks.add(tmp);
        }
        _taxonList = taxonList;
    }



    public void readMatrix(String fileName){
        try {
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            String line;

            while ((line = br.readLine()) != null) {
                if (_num_taxa == 0){
                    _num_taxa = Integer.parseInt(line.trim());
                    _matrix = new double[_num_taxa][_num_taxa];
                }
                else{
                    String[] arr = line.split(",");
                    _taxonList.add(arr[0]);
                    int row = _taxonList.size() - 1;
                    for (int i = 1; i <= _num_taxa; i++){
                        _matrix[row][i - 1] = Double.parseDouble(arr[i].trim());
                    }
                }

            }
            if (Utils._debug){
                System.out.println(_matrix);
            }

        } catch (Exception e) {
            System.err.println(_matrix);
            e.printStackTrace();
        }
    }

    public void readSubNetworks(String fileName){
        try {
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.equals("")){
                    break;
                }
                List<Network> tmp = new ArrayList<>();
                tmp.add(Networks.readNetwork(line));
                _subnetworks.add(tmp);
//                _subnetworks.add(Networks.readNetwork(line));
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private double[][] restrainedDistantMatrix(List<String> taxonList){
        int n = taxonList.size();
        double[][] matrix = new double[n][n];
        List<Integer> taxonIndexList = new ArrayList<>();
        for (String taxon: taxonList){
            taxonIndexList.add(_taxonList.indexOf(taxon));
        }
        for (int i = 0; i < n; i ++ ){
//            System.out.print(taxonList.get(i)+",");
            for (int j = 0; j < n; j++){
                matrix[i][j] = _matrix[taxonIndexList.get(i)][taxonIndexList.get(j)];
//                System.out.print(matrix[i][j]+",");
            }
//            System.out.println();
        }

        return matrix;

    }

    private List<String> getLeafList(Network net){
        List<String> leafList = new ArrayList<>();
        for (Object o: net.getLeaves()){
            leafList.add(((NetNode) o).getName());
        }
        return leafList;
    }

//    private double calcScore(Network net){
//        double topologyDistanceSum = 0;
//        for (Network subnet: _subnetworks) {
//            List<String> selectedLeaves = getLeafList(subnet);
//            Tuple<Network, Map<NetNode, NetNode>> restrainedFullNet = SuperNetwork3.getSubNetwork(net, selectedLeaves, true);
//            topologyDistanceSum += Networks.computeDistanceBetweenTwoNetworks(restrainedFullNet.Item1, subnet);
//
//        }
//        return topologyDistanceSum;
//    }

    public List<Network> mergePairs(){
        int networkCnt = _subnetworks.size();
        while (networkCnt > 1){
            int i = Randomizer.getRandomInt(networkCnt);
//            Network net1 = _subnetworks.get(i);
            List<Network> net1list = _subnetworks.get(i);
            int j = Randomizer.getRandomInt(networkCnt);
            while(j == i){
                j = Randomizer.getRandomInt(networkCnt);
            }
//            Network net2 = _subnetworks.get(j);
            List<Network> net2list = _subnetworks.get(j);
            List<Network> newcandidatelist = new ArrayList<>();
            System.out.println("list size:"+net1list.size()+" "+net2list.size());
            for (Network net1: net1list){

                for (Network net2: net2list){
                    addInheritanceProb(net1);
                    Network net1_copy = net1.clone();
                    List<Network> subnetlist = new ArrayList<>();
                    addInheritanceProb(net2);
                    Network net2_copy = net2.clone();
                    subnetlist.add(net1_copy);
                    subnetlist.add(net2_copy);
                    List<String> taxonList = new ArrayList<>();
                    for (Object o:net1_copy.getLeaves()){
                        taxonList.add(((NetNode) o).getName());
                    }
                    for (Object o:net2_copy.getLeaves()){
                        taxonList.add(((NetNode) o).getName());
                    }
//                    if (!taxonList.contains(_outgroup)){
//                        taxonList.add(_outgroup);
//                    }
                    // remove outgroup
//                    if (taxonList.contains(_outgroup)){
//                        taxonList.remove(_outgroup);
//                        NetNode outgroup = net1_copy.findNode(_outgroup);
//                        if(outgroup != null){
//                            ((NetNode)outgroup.getParents().iterator().next()).removeChild(outgroup);
//                            Networks.removeBinaryNodes(net1_copy);
//                        }
//                        else {
//                            outgroup = net2_copy.findNode(_outgroup);
//                            ((NetNode)outgroup.getParents().iterator().next()).removeChild(outgroup);
//                            Networks.removeBinaryNodes(net2_copy);
//                        }
//                    }
                    double[][] matrix = restrainedDistantMatrix(taxonList);

                    NJMergeTopology3 mergepair = new NJMergeTopology3(subnetlist, matrix, taxonList, _outgroup);
                    mergepair.setOutgroup(_outgroup);
                    System.out.println("merge pair:");
                    System.out.println(subnetlist.get(0).toString());
                    System.out.println(subnetlist.get(1).toString());
                    List<Network> candidatelist = mergepair.mergeNetsViaNJ();
                    if (candidatelist != null){
                        for (Network net: candidatelist){
                            Utils.emptyNodeLabel(net);
                        }
                        newcandidatelist.addAll(candidatelist);
                    }


                }
            }
            if (_debug){
                if (newcandidatelist.isEmpty()){
                    System.out.println("empty merge pair");

                }
                else {
                    System.out.println("new candidate list:");
                    for (Network net: newcandidatelist){
                        System.out.println(net.toString());
                    }
                }
            }

            if(i < j){
                _subnetworks.remove(j);
                _subnetworks.remove(i);
            }
            else{
                _subnetworks.remove(i);
                _subnetworks.remove(j);
            }
            _subnetworks.add(newcandidatelist);

            if (Utils._debug){
                System.out.println("********************");
            }

            networkCnt--;

        }
        if (Utils._debug){
            System.out.println(_subnetworks.get(0).toString());
        }

        return _subnetworks.get(0);

    }

    public static List<List<String>> readClusterInfo(String fileName){
        List<List<String>> clusters = new ArrayList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            String line;
            int cnt = 0;

            while ((line = br.readLine()) != null) {
                clusters.add(new ArrayList<>());
                String[] arr = line.split(",");
                for (String taxon: arr){
                    clusters.get(cnt).add(taxon);
                }
                cnt += 1;
            }

        } catch (Exception e) {

        }
        return clusters;
    }






    public static void preprocess(String treepath, String subsetspath){
        List<List<String>> subsetList = new ArrayList<>();
        List<List<String>> speciesList = new ArrayList<>();
        // read subsets
        try (BufferedReader reader = new BufferedReader(new FileReader(subsetspath))) {
            String line = reader.readLine();
            String[] parts = line.split(";");
            for (String part : parts) {
                subsetList.add(Arrays.asList(part.split(",")));
            }

            line = reader.readLine();
            String[] parts2 = line.split(";");
            for (String part : parts2) {
                speciesList.add(Arrays.asList(part.split(",")));
            }


        } catch (IOException e) {
            e.printStackTrace();
        }



        //read tree
        STITree astraltree = null;
        try (BufferedReader reader = new BufferedReader(new FileReader(treepath))) {
            String line = reader.readLine();
            astraltree = new STITree(line);
        }catch (Exception e) {
            e.printStackTrace();
        }

        // handle tree leaves



    }



//    public static void run(String distpath, String subnetworkpath){
//        NetNJMerge merge = new NetNJMerge(distpath, subnetworkpath);
//        List<Network> netlist = merge.mergePairs();
//        for (Network net: netlist){
//            net.resetRoot("Z");
//            System.out.println(net);
//        }
//
//
//    }



    public static void main(String[] args) {
//        testCase1();
//        testCase2_idealSubnet();
//        idealTests();
//        testCase3();
//        testCase4();
//        System.out.println(Randomizer.getRandomInt(100));
//        System.out.println(Randomizer.getRandomInt(100));
//        System.out.println(Randomizer.getRandomInt(100));


    }

}
