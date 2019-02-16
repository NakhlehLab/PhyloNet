package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.summary.SummaryBL;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 9/19/18
 * Time: 11:26 AM
 * To change this template use File | Settings | File Templates.
 */
public class PairwiseDistribution {
    // taxa1, taxa2, node, height
    private Map<Tuple3<String, String, String>, List<Double>> heights = new HashMap<>();
    private Map<Tuple<String, String>, List<List<Double>>> allheights = new HashMap<>();

    private void initNetHeights(Network<NetNodeInfo> network) {
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(network)) {
            if(node.getData() == null) {
                node.setData(new NetNodeInfo(0.0));
            }
            for(NetNode<NetNodeInfo> par : node.getParents()) {
                double dist = node.getParentDistance(par);
                if(par.getData() == null) {
                    par.setData(new NetNodeInfo(node.getData().getHeight() + dist));
                }
            }
        }
    }

    private static Set<NetNode> getAllAncestors(NetNode node) {
        Queue<NetNode> q = new ArrayDeque<>();
        Set<NetNode> result = new HashSet<>();
        q.add(node);
        while(!q.isEmpty()) {
            NetNode curnode = q.poll();
            for(Object parentObj : curnode.getParents()) {
                NetNode parent = (NetNode) parentObj;
                result.add(parent);
                q.add(parent);
            }
        }
        return result;
    }

    private static Set<NetNode> getMCRAs(List<NetNode> leaves) {
        Set<NetNode> curSet = null;
        for(NetNode leaf : leaves) {
            Set<NetNode> s = getAllAncestors(leaf);
            if(curSet == null) {
                curSet = s;
            } else {
                curSet.retainAll(s);
            }
        }
        return curSet;
    }

    private static Set<NetNode> getMCRAs(Network net, String s1, String s2) {
        NetNode l1 = net.findNode(s1);
        NetNode l2 = net.findNode(s2);
        List<NetNode> leaves = new ArrayList<>();
        leaves.add(l1);
        leaves.add(l2);
        return getMCRAs(leaves);
    }

    public void LoadFromFiles(List<String> filenames, String[] taxa, int chainlen, int burnin, int sample_freq) {
        List<Network> networks = new ArrayList<>();
        List<Network> samples = new ArrayList<>();

        for(String filename : filenames) {

            Network bestnet = null;

            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                int index = 0;
                String topSample = null;
                while((s = in.readLine()) != null) {
                    if(s.startsWith("Rank = 0")) {
                        topSample = s;
                    }
                    index++;
                }

                if(topSample == null) {
                    System.out.println("No samples in the file! " + filename);
                    continue;
                }

                String str = topSample.substring(topSample.indexOf(":") + 1);
                String bestNetworkString = str.substring(0, str.indexOf(";") + 1);

                bestnet = Networks.readNetwork(bestNetworkString);
                Networks.autoLabelNodes(bestnet);
                initNetHeights(bestnet);
                // networks.add(Networks.readNetwork(Networks.getTopologyString(net)));

                in.close();

                SummaryBL sbl = new SummaryBL(bestNetworkString);
                int start = (int)(burnin / sample_freq) + 1;
                int end = (int)(chainlen / sample_freq);
                sbl.addFile(filename, true, start, end);
                sbl.report(1.0, 1.0);

                samples = sbl.getSamples();
            } catch (IOException e) {
                e.printStackTrace();
            }

            for(Network net : samples) {
                if(!Networks.hasTheSameTopology(net, bestnet)) continue;
                initNetHeights(net);
                Map<NetNode, NetNode> nodemap = Networks.mapTwoNetworks(net, bestnet);
                for(int i = 0 ; i < taxa.length ; i++) {
                    if(net.findNode(taxa[i]) == null) continue;
                    for(int j = i + 1 ; j < taxa.length ; j++) {
                        if(net.findNode(taxa[j]) == null) continue;
                        Set<NetNode> MCRAs = getMCRAs(net, taxa[i], taxa[j]);
                        for(NetNode node : MCRAs) {
                            if(node.isNetworkNode()) continue;
                            Tuple3<String, String, String> key = new Tuple3<>(taxa[i], taxa[j], nodemap.get(node).getName());
                            double value = ((NetNodeInfo) node.getData()).getHeight();
                            if(!heights.containsKey(key))
                                heights.put(key, new ArrayList<>());
                            heights.get(key).add(value);
                        }
                    }
                }
            }

            for(Tuple3<String, String, String> tuple : heights.keySet()) {
                System.out.println(tuple.Item1 + " " + tuple.Item2 + " " + tuple.Item3 + " " + heights.get(tuple).size());
                Tuple<String, String> key = new Tuple<>(tuple.Item1, tuple.Item2);
                if(!allheights.containsKey(key))
                    allheights.put(key, new ArrayList<>());
                allheights.get(key).add(heights.get(tuple));
            }
        }


    }

    void PrintDataPointsToFile(String outputfilename) {
        // Output to file
        // Mixed
        try {
            PrintWriter out = new PrintWriter(outputfilename + "2");
            int index = 0;
            for(Tuple<String, String> tuple : allheights.keySet()) {
                out.println("[KEY]" + tuple.Item1 + ":" + tuple.Item2 + ":" + index);
                for(List<Double> list : allheights.get(tuple)) {
                    for(Double d : list) {
                        out.println(d);
                    }
                }
                out.println("[ENDKEY]");
                index++;
            }
            out.println("[END]");
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        // Seperated
        try {
            PrintWriter out = new PrintWriter(outputfilename);
            int index = 0;
            for(Tuple<String, String> tuple : allheights.keySet()) {
                for(List<Double> list : allheights.get(tuple)) {
                    out.println("[KEY]" + tuple.Item1 + ":" + tuple.Item2 + ":" + index);
                    for(Double d : list) {
                        out.println(d);
                    }
                    out.println("[ENDKEY]");
                    index++;
                }
            }
            out.println("[END]");
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    Map<Tuple<String, String>, List<Double>> ComputeExtendedDistanceMatrix() {
        final double eps = 1e-6;
        Map<Tuple<String, String>, List<Double>> edm = new HashMap<>();
        for(Tuple<String, String> tuple : allheights.keySet()) {
            edm.put(tuple, new ArrayList<>());
            for(List<Double> list : allheights.get(tuple)) {
                double mean = 0.0;
                System.out.println(tuple.Item1 + ":" + tuple.Item2);
                for(Double d : list) {
                    mean += d;
                }
                mean /= list.size();

                boolean found = false;
                for(double d : edm.get(tuple)) {
                    if(Math.abs(d - mean) < eps) {
                        found = true;
                        break;
                    }
                }

                if(!found)
                    edm.get(tuple).add(mean);
            }
            Collections.sort(edm.get(tuple));
        }
        return edm;
    }

    public static void main(String args[]) {
        List<String> filenames = new ArrayList<>();
        String outputfilename = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/PhyloData0/pairwise.DP";

        Network trueNetwork = Networks.readNetwork("(((((((P:4.723,(O:2.0)#H1:2.723::0.6):7.87,N:12.593):2.829,(M:6.976,(L:3.0)#H2:3.976::0.8):8.446):3.272,((K:11.378,(#H1:4.0::0.4,J:6.0):5.378):4.977,((H:10.99,((G:6.106,((I:1.0)#H3:1.106::0.7,F:2.106):4.0):3.462,E:9.568):1.422):2.74,(#H2:1.0::0.2,D:4.0):9.73):2.625):2.339):2.535,(#H3:6.229::0.3,C:7.229):14.0):6.136,B:27.365):12.913,A:40.278);");
        List<String> leaves = new ArrayList<>();
        for(Object leafObj : trueNetwork.getLeaves()) {
            NetNode leaf = (NetNode) leafObj;
            leaves.add(leaf.getName());
        }
        Collections.sort(leaves);

        String[] taxa = new String[leaves.size()];
        leaves.toArray(taxa);
        for(int i = 0 ; i < taxa.length ; i++) {
            System.out.print("\"" + taxa[i] + "\",");
        }

        //File path = new File("/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/PhyloData0");
        File path = new File("/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run7");


        File [] files = path.listFiles();
        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){ //this line weeds out other directories/folders
                if(files[i].toString().endsWith(".out")) {
                    System.out.println(files[i]);
                    filenames.add(files[i].toString());
                }
            }
        }

        Collections.sort(filenames);
        PairwiseDistribution dist = new PairwiseDistribution();
        dist.LoadFromFiles(filenames, taxa, 10000000, 1000000, 5000);
        dist.PrintDataPointsToFile(outputfilename);
    }
}
