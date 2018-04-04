package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.summary;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/13/17
 * Time: 10:31 AM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkPlotter {

    private List<String> _leafOrder;
    private Map<NetNode, List<NetNode>> _edgeOrder;
    private Network _network;
    private Map<NetNode, Double> _XCoordinates;
    private Map<NetNode, Double> _YCoordinates;
    private double _maxX;
    private double _scaleY;


    NetworkPlotter(Network network, List<String> order, double maxX, double scaleY) {
        _network = network;
        _leafOrder = order;
        _scaleY = scaleY;
        _maxX = maxX;
        computeCoordinates();
    }

    List<String> getAllLeavesBelow(NetNode<Object> node) {
        Queue<NetNode<Object>> queue = new LinkedList<>();
        queue.add(node);
        List<String> result = new ArrayList<>();
        while(!queue.isEmpty()) {
            NetNode<Object> node1 = queue.poll();
            for(NetNode<Object> child : node1.getChildren() ) {
                if(child.isLeaf()) {
                    if(!result.contains(child.getName())) {
                        result.add(child.getName());
                    }
                } else {
                    if(!queue.contains(child)) {
                        queue.add(child);
                    }
                }
            }
        }
        Collections.sort(result, new Comparator<String>() {
            public int compare(String left, String right) {
                return Integer.compare(_leafOrder.indexOf(left), _leafOrder.indexOf(right));
            }
        });
        return result;
    }

    void determineReticulationEdgeOrders() {
        _edgeOrder = new HashMap<>();

        for(Object nodeObj : Networks.postTraversal(_network)) {
            NetNode node = (NetNode) nodeObj;
            if(node.isNetworkNode()) {
                int index = 0;
                NetNode parents[] = new NetNode[2];

                for(Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    parents[index] = parent;
                    index++;
                }

                List<String> l0 = getAllLeavesBelow(parents[0]);
                List<String> l1 = getAllLeavesBelow(parents[1]);
                List<String> l2 = getAllLeavesBelow(node);

                l0.removeAll(l2);
                l1.removeAll(l2);

                int s0 = 0, s1 = 0;
                for(String s : l0) {
                    s0 += _leafOrder.indexOf(s);
                }
                for(String s : l1) {
                    s1 += _leafOrder.indexOf(s);
                }

                _edgeOrder.put(node, new ArrayList<>());
                if(s0 > s1) {
                    _edgeOrder.get(node).add(parents[1]);
                    _edgeOrder.get(node).add(parents[0]);
                } else {
                    _edgeOrder.get(node).add(parents[0]);
                    _edgeOrder.get(node).add(parents[1]);
                }
            }
        }
    }

    void computeXCoordinates() {
        _XCoordinates = new HashMap<>();
        double leafInterval = _maxX / (_leafOrder.size() - 1);
        double currentX = 0.0;
        for(String leafLabel : _leafOrder) {
            NetNode node = _network.findNode(leafLabel);
            _XCoordinates.put(node, currentX);
            currentX += leafInterval;
        }

        for(Object nodeObj : Networks.postTraversal(_network)) {
            NetNode node = (NetNode) nodeObj;
            if(node.isLeaf()) {
                continue;
            }
            currentX = 0;
            for(Object childObj : node.getChildren()) {
                NetNode child = (NetNode) childObj;
                currentX += _XCoordinates.get(child);
            }
            currentX /= node.getChildCount();
            _XCoordinates.put(node, currentX);
        }
    }

    void computeYCoordinates() {
        _YCoordinates = new HashMap<>();
        Map<NetNode, Double> heightMap = getNodeHeights(_network);
        for(NetNode node : heightMap.keySet()) {
            _YCoordinates.put(node, heightMap.get(node) * _scaleY);
        }
    }

    static public Map<NetNode, Double> getNodeHeights(Network network) {
        Map<NetNode, Double> heightMap = new HashMap<>();
        for(Object nodeObj : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeObj;
            if(node.isLeaf()) {
                heightMap.put(node, 0.0);
            }
            // store height of parents
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                double parentHeight = node.getParentDistance(parent) + heightMap.get(node);
                // check if network is ultrametric
                if(heightMap.containsKey(parent)) {
                    if(Math.abs(heightMap.get(parent) - parentHeight) > 1e-6) {
                        System.out.println("Input network is not ultrametric!");
                        System.out.println(network.toString());
                        System.exit(1);
                    }
                } else {
                    heightMap.put(parent, parentHeight);
                }
            }
        }
        return heightMap;
    }

    public void computeCoordinates() {
        computeXCoordinates();
        computeYCoordinates();
    }

    public String producePythonLists() {
        StringBuilder sb = new StringBuilder();

        sb.append("    x=[");
        for(Object nodeObj : Networks.postTraversal(_network)) {
            NetNode node = (NetNode) nodeObj;
            sb.append(_XCoordinates.get(node));
            sb.append(",");
        }
        sb.append("]\n");

        sb.append("    y=[");
        for(Object nodeObj : Networks.postTraversal(_network)) {
            NetNode node = (NetNode) nodeObj;
            sb.append(_YCoordinates.get(node));
            sb.append(",");
        }
        sb.append("]\n");

        sb.append("    lines=[");
        for(Object nodeObj : Networks.postTraversal(_network)) {
            NetNode node = (NetNode) nodeObj;
            for(Object childObj : node.getChildren()) {
                NetNode child = (NetNode) childObj;
                double x1 = _XCoordinates.get(node);
                double y1 = _YCoordinates.get(node);
                double x2 = _XCoordinates.get(child);
                double y2 = _YCoordinates.get(child);

                sb.append("[(");
                sb.append(x1);
                sb.append(",");
                sb.append(y1);
                sb.append("),(");
                sb.append(x2);
                sb.append(",");
                sb.append(y2);
                sb.append(")],");
            }
        }
        sb.append("]\n");

        sb.append("    c=[");
        for(Object nodeObj : Networks.postTraversal(_network)) {
            NetNode node = (NetNode) nodeObj;
            for(Object childObj : node.getChildren()) {
                NetNode child = (NetNode) childObj;
                if (child.isNetworkNode()) {
                    sb.append("\"blue\"");
                } else {
                    sb.append("\"black\"");
                }
                sb.append(",");
            }
        }
        sb.append("]\n");

        return sb.toString();
    }

    public String producePythonAnnotations() {
        StringBuilder sb = new StringBuilder();
        for(String leafLabel : _leafOrder) {
            NetNode leaf = _network.findNode(leafLabel);
            double x = _XCoordinates.get(leaf);
            sb.append("    ax.annotate('");
            sb.append(leafLabel);
            sb.append("', xy=(");
            sb.append(x-0.05);
            sb.append(", -0.3), xytext=(");
            sb.append(x-0.05);
            sb.append(", -0.3))\n");
        }
        return sb.toString();
    }

    public static void main(String[] args) {
        List<String> order = new ArrayList<>();
        order.add("cae");
        order.add("HYB");
        order.add("mcccal");
        order.add("mcplac");
        order.add("mccmcc");

        String netO4 = "[0.031269291383741966]((cae:0.023498520577909784:0.11453515690711917,(HYB:0.009167547496130975:0.03805835880867663)#H1:0.01433097308177881:0.03423044524104861:0.5101383020949305):0.17382587887270695:0.024235535337852533,(mccmcc:0.18439345369566149:0.06925787527744415,(mcplac:0.1465258072175068:0.08235404912736036,(#H1:0.014564359287868745:0.03939682755966605:0.4898616979050695,mcccal:0.02373190678399972:0.04882348145145685):0.12279390043350707:0.04534154430390398):0.03786764647815469:0.05099477836461758):0.012930945754955264:0.06687609363382643);";

        Network network = Networks.readNetwork(netO4);
        NetworkPlotter plotter = new NetworkPlotter(network, order, 4.0, 5.0);
        plotter.computeCoordinates();
        System.out.println(plotter.producePythonLists());
        System.out.println(plotter.producePythonAnnotations());
    }
}
