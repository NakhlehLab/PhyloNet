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

import static edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.NJMergeTopology2.match2Networks;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils.addInheritanceProb;


public class NetNJMerge2 {

    private List<List<Network>> _subnetworks = new ArrayList<>();

    private final List<Network> _originalSubNetworks = new ArrayList<>();
    private int _num_taxa = 0;
    private double[][] _matrix = null;
    private List<String> _taxonList = new ArrayList<>();
    private List<Set<String>> _leaves = null;
    private String _outgroup = "";
    private double _scale = 1;

    private Network _truenetwork = null;

    private Network _backboneTree = null;

    /* Constructor */
    public NetNJMerge2(String distPath, String netPath) {
        readMatrix(distPath);
        readSubNetworks(netPath);

    }
    public void setOutgroup(String outgroup){
        this._outgroup = outgroup;
    }
    public NetNJMerge2(String distPath, List<Network> subnetlist){
        readMatrix(distPath);
        for (Network subnet: subnetlist){
            List<Network> tmp = new ArrayList<>();
            tmp.add(subnet);
            _subnetworks.add(tmp);
            _originalSubNetworks.add(subnet.clone());
        }
    }
    public NetNJMerge2(double[][] matrix, List<Network> subnetlist, List<String> taxonList){
        _matrix = matrix;
        for (Network subnet: subnetlist){
            List<Network> tmp = new ArrayList<>();
            tmp.add(subnet);
            _subnetworks.add(tmp);
            _originalSubNetworks.add(subnet.clone());
        }
        _taxonList = taxonList;
    }

    public void setTrueNetwork(Network network){
        this._truenetwork = network;
    }

    public void setBackboneTree(Network network){
        this._backboneTree = network;
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
                Network subnet = Networks.readNetwork(line);
                tmp.add(subnet);
                _subnetworks.add(tmp);
                _originalSubNetworks.add(subnet.clone());
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
            for (int j = 0; j < n; j++){
                matrix[i][j] = _matrix[taxonIndexList.get(i)][taxonIndexList.get(j)];

            }
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
                    if (!taxonList.contains(_outgroup)){
                        taxonList.add(_outgroup);
                    }
                    double[][] matrix = restrainedDistantMatrix(taxonList);

                    NJMergeTopology2 mergepair = new NJMergeTopology2(subnetlist, matrix, taxonList, _outgroup);
                    mergepair.setOutgroup(_outgroup);
                    System.out.println("merge pair:");
                    System.out.println(subnetlist.get(0).toString());
                    System.out.println(subnetlist.get(1).toString());
                    List<Network> candidatelist = mergepair.mergeNetsViaNJ();
                    for (Network net: candidatelist){
                        Utils.emptyNodeLabel(net);
                    }
                    newcandidatelist.addAll(candidatelist);


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


        List<Network> matchedNets = matchNetworks(_subnetworks.get(0), _originalSubNetworks);


        return matchedNets;

    }
    public static List<Network> matchNetworks(List<Network> networklist, List<Network> subnetlist){
        List<Network> matchedNets = new ArrayList<>();
        Queue<Network> q = new LinkedList();
        Queue<Network> nextlevel = new LinkedList<>();


        for (Network net: networklist){
            matchedNets.add(net);
            q.add(net);
        }

        for (Network originalNet: subnetlist){
            nextlevel = new LinkedList<>();
            while (!q.isEmpty()){
                Network net = q.poll();
                List<Network> matchedlist = match2Networks(originalNet, net);
                for (Network toaddnet: matchedlist){
                    boolean inside = false;
                    for (Network inlistnet: nextlevel){
                        if (Networks.hasTheSameTopology(inlistnet, toaddnet)){
                            inside = true;
                            break;
                        }
                    }
                    if (!inside){
                        nextlevel.add(toaddnet);
                    }
                }
            }

            q = nextlevel;

        }
        matchedNets.addAll(q);
        return matchedNets;
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



    public static void writeSubNets(){
        Network net_original = Networks.readNetwork("(((((t1:1.8,#H1:1.4::0.5):2.9464,((t2:0.3052,t3:0.3052):3.658,((t4:0.4)#H1:1.0948::0.5,(t5:0.5832,t6:0.5832):0.9116):2.4684):0.7832):4.7008,(((t7:1.4428,t8:1.4428):4.4308,((t9:0.8492,t10:0.8492):4.3948,(((((t11:1.3768,t12:1.3768):0.4)#H3:0.638::0.6,(t13:0.4624,t14:0.4624):1.9524):0.5004,(t15:2.7952,(t16:2.0,#H3:0.2232::0.4):0.7952):0.12):0.2084,t17:3.1236):2.1204):0.6296):1.3796,((((t18:1.2424,t19:1.2424):0.898,t20:2.1404):2.1012,((t21:0.6644,t22:0.6644):2.2932,(t23:0.5328,t24:0.5328):2.4248):1.284):0.7936,(t25:2.3128,t26:2.3128):2.7224):2.218):2.194):1.1436,(t27:5.3544,((t28:0.6,#H2:0.2::0.4):0.8936,(t29:1.0344,(t30:0.4)#H2:0.6344::0.6):0.4592):3.8608):5.2364):89.4092,Z:100.0);");
        String dir_path = "";
        String[] scales = {"short", "medium", "long"};

        int numsites = 1000;
        for (String scale: scales){
            for (int i = 1; i <= 10; i ++){
                String distpath = dir_path+"/"+scale+"/"+i+"/"+numsites+ "/dist_matrix.txt";
                String subnetworkpath = dir_path+"/"+scale+"/"+i+"/"+numsites+ "/subnets.txt";
            }
        }
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
