//package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge;
///*
// * @ClassName:   NetMerge
// * @Description:
// * @Author:      Zhen Cao
// * @Date:        9/12/23 11:34 PM
// */
//
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.util.Randomizer;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils;
//import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
//import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
//
//import java.io.BufferedReader;
//import java.io.FileReader;
//import java.io.IOException;
//import java.util.*;
//
//import edu.rice.cs.bioinfo.library.programming.Tuple;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork3;
//import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
//import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
//
//
//public class NetNJMerge {
////    private List<Network> _subnetworks = new ArrayList<>();
//    private List<List<Network>> _subnetworks = new ArrayList<>();
//    private int _num_taxa = 0;
//    private double[][] _matrix = null;
//    private List<String> _taxonList = new ArrayList<>();
//    private List<Set<String>> _leaves = null;
//    private String _outgroup = "";
//    private double _scale = 1;
//
//    /* Constructor */
//    public NetNJMerge(String distPath, String netPath) {
//        readMatrix(distPath);
//        readSubNetworks(netPath);
//
//    }
//    public void setOutgroup(String outgroup){
//        this._outgroup = outgroup;
//    }
//    public NetNJMerge(String distPath, List<Network> subnetlist){
//        readMatrix(distPath);
//        for (Network subnet: subnetlist){
//            List<Network> tmp = new ArrayList<>();
//            tmp.add(subnet);
//            _subnetworks.add(tmp);
//        }
//    }
//    public NetNJMerge(double[][] matrix, List<Network> subnetlist, List<String> taxonList){
//        _matrix = matrix;
//        for (Network subnet: subnetlist){
//            List<Network> tmp = new ArrayList<>();
//            tmp.add(subnet);
//            _subnetworks.add(tmp);
//        }
//        _taxonList = taxonList;
//    }
//
//
//
//    public void readMatrix(String fileName){
//        try {
//            BufferedReader br = new BufferedReader(new FileReader(fileName));
//            String line;
//
//            while ((line = br.readLine()) != null) {
//                if (_num_taxa == 0){
//                    _num_taxa = Integer.parseInt(line.trim());
//                    _matrix = new double[_num_taxa][_num_taxa];
//                }
//                else{
//                    String[] arr = line.split(",");
//                    _taxonList.add(arr[0]);
//                    int row = _taxonList.size() - 1;
//                    for (int i = 1; i <= _num_taxa; i++){
//                        _matrix[row][i - 1] = Double.parseDouble(arr[i].trim());
//                    }
//                }
//
//            }
//            if (Utils._debug){
//                System.out.println(_matrix);
//            }
//
//        } catch (Exception e) {
//            System.err.println(_matrix);
//            e.printStackTrace();
//        }
//    }
//
//    public void readSubNetworks(String fileName){
//        try {
//            BufferedReader br = new BufferedReader(new FileReader(fileName));
//            String line;
//            while ((line = br.readLine()) != null) {
//                line = line.trim();
//                if (line.equals("")){
//                    break;
//                }
//                List<Network> tmp = new ArrayList<>();
//                tmp.add(Networks.readNetwork(line));
//                _subnetworks.add(tmp);
////                _subnetworks.add(Networks.readNetwork(line));
//            }
//
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//    }
//
//    private double[][] restrainedDistantMatrix(List<String> taxonList){
//        int n = taxonList.size();
//        double[][] matrix = new double[n][n];
//        List<Integer> taxonIndexList = new ArrayList<>();
//        for (String taxon: taxonList){
//            taxonIndexList.add(_taxonList.indexOf(taxon));
//        }
//        for (int i = 0; i < n; i ++ ){
//            for (int j = 0; j < n; j++){
//                matrix[i][j] = _matrix[taxonIndexList.get(i)][taxonIndexList.get(j)];
//            }
//        }
//
//        return matrix;
//
//    }
//
//    private List<String> getLeafList(Network net){
//        List<String> leafList = new ArrayList<>();
//        for (Object o: net.getLeaves()){
//            leafList.add(((NetNode) o).getName());
//        }
//        return leafList;
//    }
//
////    private double calcScore(Network net){
////        double topologyDistanceSum = 0;
////        for (Network subnet: _subnetworks) {
////            List<String> selectedLeaves = getLeafList(subnet);
////            Tuple<Network, Map<NetNode, NetNode>> restrainedFullNet = SuperNetwork3.getSubNetwork(net, selectedLeaves, true);
////            topologyDistanceSum += Networks.computeDistanceBetweenTwoNetworks(restrainedFullNet.Item1, subnet);
////
////        }
////        return topologyDistanceSum;
////    }
//
//    public List<Network> mergePairs(){
//        int networkCnt = _subnetworks.size();
//        while (networkCnt > 1){
//            int i = Randomizer.getRandomInt(networkCnt);
////            Network net1 = _subnetworks.get(i);
//            List<Network> net1list = _subnetworks.get(i);
//            int j = Randomizer.getRandomInt(networkCnt);
//            while(j == i){
//                j = Randomizer.getRandomInt(networkCnt);
//            }
////            Network net2 = _subnetworks.get(j);
//            List<Network> net2list = _subnetworks.get(j);
//            List<Network> newcandidatelist = new ArrayList<>();
//            System.out.println("list size:"+net1list.size()+" "+net2list.size());
//            for (Network net1: net1list){
//
//                for (Network net2: net2list){
//                    Network net1_copy = net1.clone();
//                    List<Network> subnetlist = new ArrayList<>();
//                    Network net2_copy = net2.clone();
//                    subnetlist.add(net1_copy);
//                    subnetlist.add(net2_copy);
//                    List<String> taxonList = new ArrayList<>();
//                    for (Object o:net1_copy.getLeaves()){
//                        taxonList.add(((NetNode) o).getName());
//                    }
//                    for (Object o:net2_copy.getLeaves()){
//                        taxonList.add(((NetNode) o).getName());
//                    }
//                    double[][] matrix = restrainedDistantMatrix(taxonList);
//                    NJMergeTopology mergepair = new NJMergeTopology(subnetlist, matrix, taxonList);
//                    mergepair.setOutgroup(_outgroup);
//                    List<Network> candidatelist = mergepair.mergeNetsViaNJ();
//                    for (Network net: candidatelist){
//                        Utils.emptyNodeLabel(net);
//                    }
//                    newcandidatelist.addAll(candidatelist);
//
//
//
//                }
//            }
//            if(i < j){
//                _subnetworks.remove(j);
//                _subnetworks.remove(i);
//            }
//            else{
//                _subnetworks.remove(i);
//                _subnetworks.remove(j);
//            }
//            _subnetworks.add(newcandidatelist);
//
//            if (Utils._debug){
//                System.out.println("********************");
//            }
//
//            networkCnt--;
//
//        }
//        if (Utils._debug){
//            System.out.println(_subnetworks.get(0).toString());
//        }
//
//        return _subnetworks.get(0);
//
//    }
//
//    public static List<List<String>> readClusterInfo(String fileName){
//        List<List<String>> clusters = new ArrayList<>();
//        try {
//            BufferedReader br = new BufferedReader(new FileReader(fileName));
//            String line;
//            int cnt = 0;
//
//            while ((line = br.readLine()) != null) {
//                clusters.add(new ArrayList<>());
//                String[] arr = line.split(",");
//                for (String taxon: arr){
//                    clusters.get(cnt).add(taxon);
//                }
//                cnt += 1;
//            }
//
//        } catch (Exception e) {
//
//        }
//        return clusters;
//    }
//
//
////    public void readDistance()
//    public static void testCase1(){
//        String distPath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/merge/code/njmerge-master/example/100-taxon-dataset/distance-renamed-rows.mat";
//        String subnetPath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/merge/code/njmerge-master/example/100-taxon-dataset/subset-5.txt";
//        NetNJMerge netmerge = new NetNJMerge(distPath, subnetPath);
//        netmerge.mergePairs();
//    }
//
//    public static void testCase2_idealSubnet(){
//        String distPath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/merge/code/distance_matrix.csv";
//        String cluster_info_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/merge/code/cluster_info.csv";
//        Network truenetwork = Networks.readNetwork("((((M:1.2,N:1.2)S14:0.528,K:1.728)S5:8.971320537907197,(((H:2.48832,G:2.48832)S11:1.8114969599999995,(F:2.9859839999999997,E:2.9859839999999997)S10:1.3138329599999996)S6:4.616283488255998,((((J:2.0736,I:2.0736)S12:1.5095807999999997,(L:1.44,(P:1.0,O:1.0)S15:0.43999999999999995)S13:2.1431807999999997)S9:1.5765995519999998,D:5.159780351999999)S8:1.0319560703999997,C:6.191736422399999)S7:2.7243640258559987)S4:1.7832200896511985)S3:2.139864107581438,(A:7.430083706879999,B:7.430083706879999)S2:5.409100938608636)S1;");
//        List<Network> subnetworklist = new ArrayList<>();
//        List<List<String>> clusters = readClusterInfo(cluster_info_path);
//        for (List<String> cluster: clusters){
//            List<String> selectedLeaves = new ArrayList<>();
//            for (String taxon: cluster){
//                selectedLeaves.add(taxon);
//            }
//            Tuple<Network, Map<NetNode, NetNode>> tuple = SuperNetwork3.getSubNetwork(truenetwork, selectedLeaves, true);
//            subnetworklist.add(tuple.Item1);
//        }
//
//        NetNJMerge netmerge = new NetNJMerge(distPath, subnetworklist);
////        Network mergednet = netmerge.mergePairs();
//        List<Network> mergednet = netmerge.mergePairs();
////        System.out.println("Distance from the true network: "+Networks.computeDistanceBetweenTwoNetworks(truenetwork, mergednet));
//    }
//
//    public static void idealTests(){
//        String dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/merge/data/";
//        String[] network_names = {"Reti0_A", "Reti0_B", "Reti0_C", "Reti0_D", "Reti1_A", "Reti1_B", "Reti1_C", "Reti1_D"};
//        String[] true_networks = {"((((M:1.2,N:1.2)S14:0.528,K:1.728)S5:8.971320537907197,(((H:2.48832,G:2.48832)S11:1.8114969599999995,(F:2.9859839999999997,E:2.9859839999999997)S10:1.3138329599999996)S6:4.616283488255998,((((J:2.0736,I:2.0736)S12:1.5095807999999997,(L:1.44,(P:1.0,O:1.0)S15:0.43999999999999995)S13:2.1431807999999997)S9:1.5765995519999998,D:5.159780351999999)S8:1.0319560703999997,C:6.191736422399999)S7:2.7243640258559987)S4:1.7832200896511985)S3:2.139864107581438,(A:7.430083706879999,B:7.430083706879999)S2:5.409100938608636)S1;",
//                "(Z:100.0,((((E:2.9859839999999997,F:2.9859839999999997)S8:4.4440997068799994,((N:1.2,M:1.2)S15:0.528,J:1.728)S7:5.702083706879999)S6:1.486016741375999,((H:2.0736,I:2.0736)S11:3.0861803519999995,((K:1.44,L:1.44)S14:1.04832,G:2.48832)S12:2.6714603519999995)S5:3.7563200962559984)S3:3.9230841972326367,((((P:1.0,O:1.0)S13:3.2998169599999994,B:4.299816959999999)S10:1.8919194623999998,(D:3.5831807999999996,C:3.5831807999999996)S9:2.6085556223999995)S4:4.507584115507197,A:10.699320537907196)S2:2.139864107581438)S1:87.16081535451137);",
//                "(Z:100.0,((((N:6.191736422399999,B:6.191736422399999)S15:1.2383472844799996,((O:1.0,P:1.0)S14:2.5831807999999996,D:3.5831807999999996)S13:3.846902906879999)S12:1.486016741375999,(C:5.159780351999999,(((G:2.0736,(I:1.728,H:1.728)S11:0.3455999999999999)S10:0.9123839999999999,(J:1.44,K:1.44)S9:1.5459839999999998)S8:1.3138329599999996,(M:1.2,L:1.2)S7:3.099816959999999)S6:0.859963392)S5:3.7563200962559984)S4:3.9230841972326367,(A:10.699320537907196,(F:2.48832,E:2.48832)S3:8.211000537907196)S2:2.139864107581438)S1:87.16081535451137);",
//                "(Z:100.0,(((E:2.9859839999999997,((L:1.44,K:1.44)S15:0.6335999999999999,H:2.0736)S14:0.9123839999999999)S9:3.2057524223999994,(G:2.48832,F:2.48832)S8:3.7034164223999992)S2:6.647448223088635,((A:8.916100448255998,((C:4.299816959999999,(D:3.5831807999999996,(J:1.728,I:1.728)S13:1.8551807999999996)S12:0.7166361599999997)S11:0.859963392,(O:1.0,P:1.0)S10:4.159780351999999)S6:3.7563200962559984)S5:1.7832200896511985,(B:7.430083706879999,(N:1.2,M:1.2)S7:6.230083706879999)S4:3.2692368310271975)S3:2.139864107581438)S1:87.16081535451137);",
//                "(Z:100.0,(((((E:4.299816959999999)#H1:3.1302667468799994::0.4,(D:5.159780351999999,#H1:0.859963392::0.6)S16:2.2703033548799993)S15:3.2692368310271975,(F:3.5831807999999996,(H:2.9859839999999997,G:2.9859839999999997)S14:0.5971967999999999)S13:7.116139737907197)S12:2.139864107581438,A:12.839184645488634)S11:5.649241244014997,(((M:1.2,N:1.2)S4:0.24,L:1.44)S3:13.967021574586362,(((((O:1.0,P:1.0)S10:0.728,K:1.728)S9:0.3455999999999999,J:2.0736)S8:0.41472,I:2.48832)S7:6.427780448255998,(C:6.191736422399999,B:6.191736422399999)S6:2.7243640258559987)S5:6.490921126330363)S2:3.08140431491727)S1:81.51157411049637);",
//                "(Z:100.0,(A:18.48842588950363,(((((K:1.728,L:1.728)S13:4.463736422399999,B:6.191736422399999)S12:1.2383472844799996,(M:1.44,(N:1.2)#H1:0.24::0.4)S11:5.990083706879998)S10:1.486016741375999,(#H1:3.9597803519999992::0.6,((J:2.0736,I:2.0736)S15:1.5095807999999997,E:3.5831807999999996)S14:1.5765995519999998)S9:3.7563200962559984)S3:6.490921126330363,((((H:2.48832,G:2.48832)S16:0.4976639999999999,F:2.9859839999999997)S8:7.713336537907196,(O:1.0,P:1.0)S7:9.699320537907196)S6:2.139864107581438,(C:4.299816959999999,D:4.299816959999999)S5:8.539367685488635)S4:2.5678369290977265)S2:3.08140431491727)S1:81.51157411049637);",
//                "(Z:100.0,(((C:7.430083706879999,((N:1.2,M:1.2)S13:2.3831808)#H1:3.846902906879999::0.4)S7:3.2692368310271975,A:10.699320537907196)S2:7.789105351596435,(((((P:1.0,O:1.0)S16:0.728,J:1.728)S14:1.2579839999999998,(K:1.44,L:1.44)S15:1.5459839999999998)S11:2.1737963519999997,(G:2.48832,F:2.48832)S10:2.6714603519999995)S4:10.24724122258636,(((((H:2.0736,I:2.0736)S12:2.2262169599999995,E:4.299816959999999)S9:1.8919194623999998,D:6.191736422399999)S8:2.7243640258559987,B:8.916100448255998)S6:3.9230841972326367,#H1:9.256003845488635::0.6)S5:2.5678369290977265)S3:3.08140431491727)S1:81.51157411049637);",
//                "(Z:100.0,((((((C:7.430083706879999,(D:6.191736422399999,(J:1.728,(K:1.44,L:1.44)S16:0.28800000000000003)S8:4.463736422399999)S7:1.2383472844799996)S6:1.486016741375999,(E:5.159780351999999,(((((P:1.0,O:1.0)S15:1.0735999999999999,(N:1.2,M:1.2)S14:0.8735999999999999)S13:0.41472,I:2.48832)S12:0.4976639999999999,H:2.9859839999999997)S11:1.3138329599999996,(F:3.5831807999999996,G:3.5831807999999996)S10:0.7166361599999997)S9:0.859963392)S5:3.7563200962559984)S4:1.7832200896511985,B:10.699320537907196)S3:2.139864107581438)#H1:2.5678369290977265::0.4,A:15.407021574586361)S2:3.08140431491727,#H1:5.649241244014997::0.6)S1:81.51157411049637);"};
//        for (int i = 1; i < network_names.length; i++){
//
//            String netName = network_names[i];
//            String trueNet = true_networks[i];
//            String dist_path = dir_path+netName+"/distance_matrix.csv";
//            String cluster_info_path = dir_path+netName+"/cluster_info.txt";
//            List<Network> subnetworklist = new ArrayList<>();
//            List<List<String>> clusters = readClusterInfo(cluster_info_path);
//            Network truenetwork = Networks.readNetwork(trueNet);
//
//            System.out.println(netName);
//            for (List<String> cluster: clusters){
//                List<String> selectedLeaves = new ArrayList<>();
//                for (String taxon: cluster){
//                    selectedLeaves.add(taxon);
//                }
//                Tuple<Network, Map<NetNode, NetNode>> tuple = SuperNetwork3.getSubNetwork(truenetwork, selectedLeaves, true);
//                subnetworklist.add(tuple.Item1);
//            }
//
//            NetNJMerge netmerge = new NetNJMerge(dist_path, subnetworklist);
////            Network mergednet = netmerge.mergePairs();
//            List<Network> mergednetlist = netmerge.mergePairs();
//            for (Network mergednet: mergednetlist){
//                System.out.println("Distance from the true network: "+Networks.computeDistanceBetweenTwoNetworks(truenetwork, mergednet));
//            }
//
//
//        }
//    }
//
//    public static void testCase3(){
//        Network net_original = Networks.readNetwork("(((((t1:1.8,#H1:1.4::0.5):2.9464,((t2:0.3052,t3:0.3052):3.658,((t4:0.4)#H1:1.0948::0.5,(t5:0.5832,t6:0.5832):0.9116):2.4684):0.7832):4.7008,(((t7:1.4428,t8:1.4428):4.4308,((t9:0.8492,t10:0.8492):4.3948,(((((t11:1.3768,t12:1.3768):0.4)#H3:0.638::0.6,(t13:0.4624,t14:0.4624):1.9524):0.5004,(t15:2.7952,(t16:2.0,#H3:0.2232::0.4):0.7952):0.12):0.2084,t17:3.1236):2.1204):0.6296):1.3796,((((t18:1.2424,t19:1.2424):0.898,t20:2.1404):2.1012,((t21:0.6644,t22:0.6644):2.2932,(t23:0.5328,t24:0.5328):2.4248):1.284):0.7936,(t25:2.3128,t26:2.3128):2.7224):2.218):2.194):1.1436,(t27:5.3544,((t28:0.6,#H2:0.2::0.4):0.8936,(t29:1.0344,(t30:0.4)#H2:0.6344::0.6):0.4592):3.8608):5.2364):89.4092,Z:100.0);");
//        List<List<String>> cladeslist = new ArrayList<>();
//        List<String> taxonlist = new ArrayList();
//        for(Object o: net_original.getLeaves()){
//            NetNode n = (NetNode) o;
//            taxonlist.add(n.getName());
//        }
//        List<Network> netlist = new ArrayList<>();
//
////        "[['t1', 't2', 't3', 't4', 't5', 't6'], ['Z', 't27', 't28', 't29', 't30'], ['t10', 't11', 't12', 't13', 't14', 't15', 't16', 't17', 't7', 't8', 't9'], ['t18', 't19', 't20', 't21', 't22', 't23', 't24', 't25', 't26']]"
//        cladeslist.add(Arrays.asList("t1", "t2", "t3", "t4", "t5", "t6"));
//        cladeslist.add(Arrays.asList("Z", "t27", "t28", "t29", "t30"));
//        cladeslist.add(Arrays.asList("t10", "t11", "t12", "t13", "t14", "t15", "t16", "t17", "t7", "t8", "t9"));
//        cladeslist.add(Arrays.asList("t18", "t19", "t20", "t21", "t22", "t23", "t24", "t25", "t26"));
//        for (List<String> clade: cladeslist){
//            Tuple<Network, Map<NetNode, NetNode>> subnet1 = SuperNetwork3.getSubNetwork(net_original, clade, true);
//            System.out.println(subnet1.Item1.toString());
//            netlist.add(subnet1.Item1);
//        }
//        System.out.println(netlist);
//
//
//        String[] taxonarr = {"t7", "t8", "t9", "t10", "t17", "t15", "t13", "t14", "t16", "t11", "t12", "t26", "t25", "t20", "t19", "t18", "t24", "t23", "t22", "t21", "t1", "t3", "t2", "t4", "t5", "t6", "Z", "t27", "t29", "t28", "t30"};
//        taxonlist = Arrays.asList(taxonarr);
//        double [][] matrix = {{0, 3, 6, 6, 6, 7, 9, 9, 9, 10, 10, 7, 7, 8, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9, 10, 10, 7, 8, 9, 10, 10},
//                {3, 0, 5, 5, 5, 6, 8, 8, 8, 9, 9, 6, 6, 7, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 9, 9, 6, 7, 8, 9, 9},
//                {6, 5, 0, 2, 4, 5, 7, 7, 7, 8, 8, 7, 7, 8, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9, 10, 10, 7, 8, 9, 10, 10},
//                {6, 5, 2, 0, 4, 5, 7, 7, 7, 8, 8, 7, 7, 8, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9, 10, 10, 7, 8, 9, 10, 10},
//                {6, 5, 4, 4, 0, 3, 5, 5, 5, 6, 6, 7, 7, 8, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9, 10, 10, 7, 8, 9, 10, 10},
//                {7, 6, 5, 5, 3, 0, 4, 4, 4, 5, 5, 8, 8, 9, 10, 10, 10, 10, 10, 10, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 7, 7, 5, 4, 0, 2, 4, 5, 5, 10, 10, 11, 12, 12, 12, 12, 12, 12, 10, 12, 12, 12, 13, 13, 10, 11, 12, 13, 13}, {9, 8, 7, 7, 5, 4, 2, 0, 4, 5, 5, 10, 10, 11, 12, 12, 12, 12, 12, 12, 10, 12, 12, 12, 13, 13, 10, 11, 12, 13, 13}, {9, 8, 7, 7, 5, 4, 4, 4, 0, 3, 3, 10, 10, 11, 12, 12, 12, 12, 12, 12, 10, 12, 12, 12, 13, 13, 10, 11, 12, 13, 13}, {10, 9, 8, 8, 6, 5, 5, 5, 3, 0, 2, 11, 11, 12, 13, 13, 13, 13, 13, 13, 11, 13, 13, 13, 14, 14, 11, 12, 13, 14, 14}, {10, 9, 8, 8, 6, 5, 5, 5, 3, 2, 0, 11, 11, 12, 13, 13, 13, 13, 13, 13, 11, 13, 13, 13, 14, 14, 11, 12, 13, 14, 14}, {7, 6, 7, 7, 7, 8, 10, 10, 10, 11, 11, 0, 2, 5, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 9, 9, 6, 7, 8, 9, 9}, {7, 6, 7, 7, 7, 8, 10, 10, 10, 11, 11, 2, 0, 5, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 9, 9, 6, 7, 8, 9, 9}, {8, 7, 8, 8, 8, 9, 11, 11, 11, 12, 12, 5, 5, 0, 3, 3, 5, 5, 5, 5, 7, 9, 9, 9, 10, 10, 7, 8, 9, 10, 10}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 3, 0, 2, 6, 6, 6, 6, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 3, 2, 0, 6, 6, 6, 6, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 5, 6, 6, 0, 2, 4, 4, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 5, 6, 6, 2, 0, 4, 4, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 5, 6, 6, 4, 4, 0, 2, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 5, 6, 6, 4, 4, 2, 0, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {7, 6, 7, 7, 7, 8, 10, 10, 10, 11, 11, 6, 6, 7, 8, 8, 8, 8, 8, 8, 0, 4, 4, 4, 5, 5, 4, 5, 6, 7, 7}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 8, 8, 9, 10, 10, 10, 10, 10, 10, 4, 0, 2, 4, 5, 5, 6, 7, 8, 9, 9}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 8, 8, 9, 10, 10, 10, 10, 10, 10, 4, 2, 0, 4, 5, 5, 6, 7, 8, 9, 9}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 8, 8, 9, 10, 10, 10, 10, 10, 10, 4, 4, 4, 0, 3, 3, 6, 7, 8, 9, 9}, {10, 9, 10, 10, 10, 11, 13, 13, 13, 14, 14, 9, 9, 10, 11, 11, 11, 11, 11, 11, 5, 5, 5, 3, 0, 2, 7, 8, 9, 10, 10}, {10, 9, 10, 10, 10, 11, 13, 13, 13, 14, 14, 9, 9, 10, 11, 11, 11, 11, 11, 11, 5, 5, 5, 3, 2, 0, 7, 8, 9, 10, 10}, {7, 6, 7, 7, 7, 8, 10, 10, 10, 11, 11, 6, 6, 7, 8, 8, 8, 8, 8, 8, 4, 6, 6, 6, 7, 7, 0, 3, 4, 5, 5}, {8, 7, 8, 8, 8, 9, 11, 11, 11, 12, 12, 7, 7, 8, 9, 9, 9, 9, 9, 9, 5, 7, 7, 7, 8, 8, 3, 0, 3, 4, 4}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 8, 8, 9, 10, 10, 10, 10, 10, 10, 6, 8, 8, 8, 9, 9, 4, 3, 0, 3, 3}, {10, 9, 10, 10, 10, 11, 13, 13, 13, 14, 14, 9, 9, 10, 11, 11, 11, 11, 11, 11, 7, 9, 9, 9, 10, 10, 5, 4, 3, 0, 2}, {10, 9, 10, 10, 10, 11, 13, 13, 13, 14, 14, 9, 9, 10, 11, 11, 11, 11, 11, 11, 7, 9, 9, 9, 10, 10, 5, 4, 3, 2, 0}};
//        for (int i = 0; i < matrix.length; i++){
//            for (int j =0; j < matrix.length; j++){
//                matrix[i][j] *= 2;
//            }
//        }
//
//        NetNJMerge merge = new NetNJMerge(matrix, netlist, taxonlist);
//        List<Network> networklist = merge.mergePairs();
//        for (Network net: networklist){
//            System.out.println(net.toString());
//            net.resetRoot("Z");
//            System.out.println("rerooted network");
//            System.out.println(net.toString());
//            System.out.println(Networks.computeDistanceBetweenTwoNetworks(net, net_original));
//        }
//
//    }
//
//    public static void testCase4(){
//        String distpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/merge/data/replica/long/1/1000/dist_matrix.txt";
//        String subnetworkpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/merge/data/replica/long/1/1000/subnets.txt";
//        NetNJMerge merge = new NetNJMerge(distpath, subnetworkpath);
////        Network net = merge.mergePairs();
//        List<Network> networklist = merge.mergePairs();
//        for (Network net: networklist){
//            net.resetRoot("Z");
//            System.out.println(net);
//            Network net_original = Networks.readNetwork("(((((t1:1.8,#H1:1.4::0.5):2.9464,((t2:0.3052,t3:0.3052):3.658,((t4:0.4)#H1:1.0948::0.5,(t5:0.5832,t6:0.5832):0.9116):2.4684):0.7832):4.7008,(((t7:1.4428,t8:1.4428):4.4308,((t9:0.8492,t10:0.8492):4.3948,(((((t11:1.3768,t12:1.3768):0.4)#H3:0.638::0.6,(t13:0.4624,t14:0.4624):1.9524):0.5004,(t15:2.7952,(t16:2.0,#H3:0.2232::0.4):0.7952):0.12):0.2084,t17:3.1236):2.1204):0.6296):1.3796,((((t18:1.2424,t19:1.2424):0.898,t20:2.1404):2.1012,((t21:0.6644,t22:0.6644):2.2932,(t23:0.5328,t24:0.5328):2.4248):1.284):0.7936,(t25:2.3128,t26:2.3128):2.7224):2.218):2.194):1.1436,(t27:5.3544,((t28:0.6,#H2:0.2::0.4):0.8936,(t29:1.0344,(t30:0.4)#H2:0.6344::0.6):0.4592):3.8608):5.2364):89.4092,Z:100.0);");
//            System.out.println(Networks.computeDistanceBetweenTwoNetworks(net, net_original));
//        }
//
//    }
//
//    public static void writeSubNets(){
//        Network net_original = Networks.readNetwork("(((((t1:1.8,#H1:1.4::0.5):2.9464,((t2:0.3052,t3:0.3052):3.658,((t4:0.4)#H1:1.0948::0.5,(t5:0.5832,t6:0.5832):0.9116):2.4684):0.7832):4.7008,(((t7:1.4428,t8:1.4428):4.4308,((t9:0.8492,t10:0.8492):4.3948,(((((t11:1.3768,t12:1.3768):0.4)#H3:0.638::0.6,(t13:0.4624,t14:0.4624):1.9524):0.5004,(t15:2.7952,(t16:2.0,#H3:0.2232::0.4):0.7952):0.12):0.2084,t17:3.1236):2.1204):0.6296):1.3796,((((t18:1.2424,t19:1.2424):0.898,t20:2.1404):2.1012,((t21:0.6644,t22:0.6644):2.2932,(t23:0.5328,t24:0.5328):2.4248):1.284):0.7936,(t25:2.3128,t26:2.3128):2.7224):2.218):2.194):1.1436,(t27:5.3544,((t28:0.6,#H2:0.2::0.4):0.8936,(t29:1.0344,(t30:0.4)#H2:0.6344::0.6):0.4592):3.8608):5.2364):89.4092,Z:100.0);");
//        String dir_path = "";
//        String[] scales = {"short", "medium", "long"};
//
//        int numsites = 1000;
//        for (String scale: scales){
//            for (int i = 1; i <= 10; i ++){
//                String distpath = dir_path+"/"+scale+"/"+i+"/"+numsites+ "/dist_matrix.txt";
//                String subnetworkpath = dir_path+"/"+scale+"/"+i+"/"+numsites+ "/subnets.txt";
//            }
//        }
//    }
//
//    public static void preprocess(String treepath, String subsetspath){
//        List<List<String>> subsetList = new ArrayList<>();
//        List<List<String>> speciesList = new ArrayList<>();
//        // read subsets
//        try (BufferedReader reader = new BufferedReader(new FileReader(subsetspath))) {
//            String line = reader.readLine();
//            String[] parts = line.split(";");
//            for (String part : parts) {
//                subsetList.add(Arrays.asList(part.split(",")));
//            }
//
//            line = reader.readLine();
//            String[] parts2 = line.split(";");
//            for (String part : parts2) {
//                speciesList.add(Arrays.asList(part.split(",")));
//            }
//            
//
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//
//
//
//        //read tree
//        STITree astraltree = null;
//        try (BufferedReader reader = new BufferedReader(new FileReader(treepath))) {
//            String line = reader.readLine();
//            astraltree = new STITree(line);
//        }catch (Exception e) {
//            e.printStackTrace();
//        }
//
//        // handle tree leaves
//
//
//
//    }
//
//
//
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
//
//
//
//    public static void main(String[] args) {
////        testCase1();
////        testCase2_idealSubnet();
////        idealTests();
////        testCase3();
////        testCase4();
//
//        
//    }
//
//}
