package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.summary;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

/**
 * Created by wendingqiao on 10/17/14.
 */
public class Summary<T> {

    private List<List<Tuple<String,Double>>> _samples;

    private Map<Network, List<Map<NetNode<T>, NetNode<T>>>> netTopologies;

    private Map<Network, OperatorLog> information;

    private PriorityQueue<Network> consensusTopologies;

    private boolean _ultrametric;

    private Map<Network, Map<NetNode, Double>> netHeightMap;

    //private double ratio;

    public Summary(List<List<Tuple<String,Double>>> samples, boolean ultrametric) {
        this._samples = samples;
        this._ultrametric = ultrametric;
        if(_ultrametric) {
            netHeightMap = new HashMap<Network, Map<NetNode, Double>>();
        }
        //this.ratio = ratio;
        netTopologies = new HashMap<Network, List<Map<NetNode<T>, NetNode<T>>>>();
        information = new HashMap<>();
        rankSpeciesTopologies();
    }

    /**
     * Given samples, rank species network topologies
     * - Build network topology map
     * - Construct priority queue of network topologies
     */
    private void rankSpeciesTopologies() {
        // build network topology map
        long timeHasCycle = 0;
        long timeMapping = 0;

        double maxPost = Double.NEGATIVE_INFINITY;
        String maxTopo = null;

        //for(int i = (int)(_samples.size() * ratio); i < _samples.size(); i++) {
        for(int i = 0; i < _samples.size(); i++) {
            Tuple<String,Double> tuple = _samples.get(i).get(0);
            if(tuple.Item2 > maxPost) {
                maxPost = tuple.Item2;
                maxTopo = tuple.Item1;
            }
            Network n = Networks.readNetwork(tuple.Item1);
            boolean find = false;
            // go over map to see whether this topology has exist
            for(Network network : netTopologies.keySet()) {
                // if find, add map to the list
                // compute time
                timeMapping -= System.currentTimeMillis();
                Map<NetNode<T>, NetNode<T>> tmpMap = Networks.mapTwoNetworks(network, n);
                boolean t1 = tmpMap != null;
                timeMapping += System.currentTimeMillis();

                timeHasCycle -= System.currentTimeMillis();
                boolean t2 = Networks.hasTheSameTopology(network, n);
                timeHasCycle += System.currentTimeMillis();
                // exception
                if(t1 != t2) {
                    System.out.println(network.toString());
                    System.out.println(n.toString());
                    System.out.println("!!!!!!!!Exception: hasSameTopology() != mapResult() !!!!!!!!");
                }
                // is same topology
                if(t1 && t2) {
                    netTopologies.get(network).add(tmpMap);
                    if(_ultrametric){
                        netHeightMap.put(network, addHeightMap(tmpMap, netHeightMap.get(network), getNodeHeightMap(n)));
                    }
                    find = true;
                    information.get(network).addSample(tuple);
                    break;
                }
            }
            // else create a new topology list.
            if(!find) {
                List<Map<NetNode<T>, NetNode<T>>> list = new ArrayList<Map<NetNode<T>, NetNode<T>>>();
                netTopologies.put(n, list);
                if(_ultrametric) {
                    netHeightMap.put(n, getNodeHeightMap(n));
                }
                information.put(n, new OperatorLog(tuple));
            }
        }
        System.out.println("time has cycle: " + timeHasCycle + "  time mapping: " + timeMapping);
        // build priority queue
        consensusTopologies = new PriorityQueue<Network>(netTopologies.size(),
                new Comparator<Network>() {
                    // the rank of networks is descending by list size
                    @Override
                    public int compare(Network o1, Network o2) {
                        return netTopologies.get(o2).size() - netTopologies.get(o1).size();
                    }
                });
        for(Network n : netTopologies.keySet()) {
            consensusTopologies.add(n);
        }
        System.out.println("Overall MAP = " + maxPost + "\n" + maxTopo);
    }

    /**
     * average branch length and probability of same topology given network
     * @param net   given network
     * @return    averaged network
     */
    private void averageNetworks(Network net) {
        List<Map<NetNode<T>, NetNode<T>>> list = netTopologies.get(net);
        for(Map<NetNode<T>, NetNode<T>> map : list) {
            bfs(net, "add", map); // traverse net - add branch length and probability
        }
        double num = list.size() + 1;
        bfs(net, "ave", num); // traverse net - average branch length and probability
        if(_ultrametric) {
            averageHeight(netHeightMap.get(net), num);
            setBranchLength(netHeightMap.get(net));
        }
    }

    private Map<NetNode<T>, Double> probabilityMap = new HashMap<>();

    public void bfs(Network net, String type, Object param) {
        Queue<NetNode<T>> queue = new LinkedList<NetNode<T>>();
        queue.add(net.getRoot());
        Set<NetNode<T>> visited = new HashSet<NetNode<T>>();
        while(!queue.isEmpty()) {
            NetNode<T> cur = queue.poll();
            if(visited.contains(cur)) continue;
            visited.add(cur);
            for(NetNode<T> c : cur.getChildren()) {
                if (type.equals("add")) {
                    Map<NetNode<T>, NetNode<T>> map = (Map<NetNode<T>, NetNode<T>>) param;
                    double pr = map.get(c).getParentProbability(map.get(cur));
                    double bl = map.get(c).getParentDistance(map.get(cur));
                    c.setParentDistance(cur, c.getParentDistance(cur) + bl);
                    // todo: This causes problems when parent probability is exactly 0.5
                    if(c.isNetworkNode() && c.getParentProbability(cur) < 0.50) {
                        if(!probabilityMap.containsKey(c)) {
                            probabilityMap.put(c, c.getParentProbability(cur));
                        } else {
                            probabilityMap.put(c, probabilityMap.get(c) + pr);
                        }
                    }
                } else if (type.equals("ave")) {
                    double size = (double)param;
                    c.setParentDistance(cur, c.getParentDistance(cur) / size);
                    if(c.isNetworkNode()) {
                        double prob = probabilityMap.get(c) / size;
                        if(c.getParentProbability(cur) > 0.50) prob = 1.0 - prob;
                        c.setParentProbability(cur, prob);
                    }
                }
                queue.add(c);
            }
        }
    }

    /**
     * Gets the top k number of average networks
     * @param k   top k
     * @return    the string presentation of the networks
     */
    public String getTopK(int k) {
        StringBuilder sb = new StringBuilder("");
        int i = 0;
        while(!consensusTopologies.isEmpty()) {
            Network net = consensusTopologies.poll();
            sb.append("Rank = " + i++ + "; ");
            int size = getSize(net);
            sb.append("Size = " + getSize(net) + "; ");
            String percent = Double.toString((double) size / (double) _samples.size());
            sb.append("Percent = " + percent + "; ");
            sb.append("MAP = " + information.get(net).toString(size) + "; ");
            sb.append("Topology = " + Networks.getTopologyString(net) + " ");
			if (getSize(net) > 1) {
				averageNetworks(net);
			}
            sb.append(net.toString());
            sb.append(" Dendroscope Visualization = " + Networks.getDendroscopeCompatibleString(net) + " ");
            sb.append("\n");
            if(consensusTopologies.size() == 0) break;
        }
        return sb.toString();
    }

    private int getSize(Network net) {
        return netTopologies.get(net).size()+1;
    }

    // get network height map
    private Map<NetNode, Double> getNodeHeightMap(Network network) {
        Map<NetNode, Double> networkHeights = new HashMap<NetNode, Double>();
        for(Object nodeO: Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeO;
            if(node.isLeaf()){
                networkHeights.put(node, 0.0);
                continue;
            }
            double height = 0.0;
            for(Object childO: node.getChildren()){
                NetNode childNode = (NetNode)childO;
                double ht = childNode.getParentDistance(node) + networkHeights.get(childNode);
                height = Math.max(height, ht);
            }
            networkHeights.put(node, height);
        }
        return networkHeights;
    }

    private void setBranchLength(Map<NetNode, Double> heightMap) {
        for(NetNode node : heightMap.keySet()) {
            for(Object c : node.getChildren()) {
                NetNode child = (NetNode)c;
                if(!heightMap.containsKey(child)) {
                    throw new IllegalArgumentException("Cannot find the child " + child.getName()
                            + " of " + node.getName() + " in heightMap.");
                }
                child.setParentDistance(node, heightMap.get(node) - heightMap.get(child));
            }
        }
    }

    private void averageHeight(Map<NetNode, Double> heightMap, double num) {
        for(NetNode node : heightMap.keySet()) {
            heightMap.put(node, heightMap.get(node) / num);
        }
    }


    private Map<NetNode, Double> addHeightMap(Map<NetNode<T>, NetNode<T>> n12n2, Map<NetNode, Double> heightMap1, Map<NetNode, Double> heightMap2) {
        Map<NetNode, Double> heightMap = new HashMap<NetNode, Double>();
        for(NetNode n1 : heightMap1.keySet()) {
            NetNode n2 = n12n2.get(n1);
            double height = heightMap1.get(n1) + heightMap2.get(n2);
            heightMap.put(n1, height);
        }
        return heightMap;
    }

    public static void main(String[] args) {
        //do the summary for unfinished job
        String filename = "/scratch/jz55/butterfly/run_2ndwave/slurm-3425287_3.out";

        List<List<Tuple<String,Double>>> netSamples = new ArrayList<>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line;
            String inputfile = "";
            int seqLength;
            String allowConstantLabel;
            String currentChain = "";
            List<Network> trueBackbones = new ArrayList<>();
            boolean ready = false;
            boolean finished = false;
            double bestLogPosterior = -1e18;

            while ((line = br.readLine()) != null) {
                if(line.contains("/") && inputfile.equals("")) {
                    inputfile = line.substring(line.lastIndexOf('/') + 1);
                    System.out.println("Last Run Input: " + inputfile);
                }
                else if(line.startsWith("Temp")) {
                    currentChain = line;
                    ready = true;
                }
                else if(ready) {
                    Scanner lineScanner = new Scanner(line);
                    lineScanner.useDelimiter(";    ");
                    int numSamples = lineScanner.nextInt();
                    double logPosterior = lineScanner.nextDouble();
                    if(currentChain.contains("main")) {
                        String netstring = br.readLine();
                        List<Tuple<String,Double>> netSample = new ArrayList<>();
                        netSample.add(new Tuple<>(netstring, new Double(logPosterior)));
                        netSamples.add(netSample);
                    }
                    ready = false;
                }
                else if(line.startsWith("Rank = 0")){
                    System.out.println("Last Run Finished ");
                }
            }

            System.out.println("Got: " + netSamples.size());

        } catch (Exception e) {
            e.printStackTrace();
        }

        Summary summary = new Summary(netSamples, true);
        System.out.println("         -------------- Top Topologies: --------------");
        System.out.println(summary.getTopK(Utils._TOPK_NETS));
    }
}