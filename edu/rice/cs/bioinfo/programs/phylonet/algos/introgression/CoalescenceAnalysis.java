package edu.rice.cs.bioinfo.programs.phylonet.algos.introgression;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.*;

/**
 * Created by wendingqiao on 3/11/15.
 */
public class CoalescenceAnalysis {

    private Network _network;
    private List<List<MutableTuple<Tree,Double>>> _geneTrees;

    public CoalescenceAnalysis(Network net, List<List<MutableTuple<Tree,Double>>> gts) {
        this._network = net;
        this._geneTrees = gts;
    }


    public String computeIntrogressionRate() {
        Network net = Networks.readNetwork(_network.toString());
        StringBuilder sb = new StringBuilder("\n");
        sb.append(net.toString());
        sb.append("\n");
        // Process Network
        Map<Double, Double> rateMapping = new HashMap<>();
        Map<Double, String> edgeMapping = new HashMap<>();

        TreeMap<Double, NetNode<Double>> map = new TreeMap<>();
        for(Object o : net.getNetworkNodes())  {
            NetNode<Double> node = (NetNode<Double>)o;
            for(Object p : node.getParents()) {
                NetNode<Double> par = (NetNode<Double>)p;
                if(node.getParentProbability(par) > 0.5) continue;
                map.put(node.getParentProbability(par), node);
            }
        }
        int i = map.size();
        for(Double d : map.keySet()) {
            NetNode<Double> node = map.get(d);
            for(Object p : node.getParents()) {
                NetNode<Double> par = (NetNode<Double>)p;
                if(node.getParentProbability(par) > 0.5) {
                    node.setParentProbability(par, 1.0 + 0.000001 * i);
                } else {
                    node.setParentProbability(par, 1.0 - 0.000001 * i);
                    rateMapping.put(1.0 - 0.000001 * i, d);
                    edgeMapping.put(d, par.getName() + "->" + node.getName());
                }
            }
            i--;
        }
        // computation
        for(List<MutableTuple<Tree,Double>> gtPerLocus : _geneTrees) {
            List<Tree> input = new ArrayList<>();
            for(MutableTuple<Tree,Double> mt : gtPerLocus) {
                input.add(mt.Item1);
            }
            Map<Double, Double> result = new TreeMap<>();

            Coalescence coal = new Coalescence(net, input, null);
            List<Map<Double,Double>> introgressionList = coal.getInheritanceProbList();
            List<Double> gtProbList = coal.getGtProb();

            double totalProb = 0.0;
            for(i = 0; i < gtPerLocus.size(); i++) {
                double weight = gtPerLocus.get(i).Item2;
                totalProb += gtProbList.get(i) * weight;
                Map<Double,Double> tmp = introgressionList.get(i);
//                System.out.println(weight + " " + gtProbList.get(i) + " " + report(tmp, 1.0));
                for(Double key : tmp.keySet()) {
                    Double mappedKey = rateMapping.get(key);
                    if(!result.containsKey(mappedKey)) {
                        result.put(mappedKey, tmp.get(key) * weight);
                    } else {
                        result.put(mappedKey, result.get(mappedKey) + tmp.get(key) * weight);
                    }
                }
            }
            sb.append(report(result, totalProb, edgeMapping));
        }
        return sb.toString();
    }

    private String report(Map<Double,Double> map, double div, Map<Double, String> edges) {
        StringBuilder sb = new StringBuilder(" ");
        for(Double d : map.keySet()) {
            sb.append(String.format(" %s:%.4f; ", edges.get(d), map.get(d) / div));
        }
        sb.append("\n");
        return sb.toString();
    }


    public String computeXLG() {
        // Process Network
        Network net = Networks.readNetwork(_network.toString());
        for(Object o : net.getNetworkNodes()) {
            NetNode<Double> node = (NetNode<Double>) o;
            for(Object p : node.getParents()) {
                NetNode<Double> par = (NetNode<Double>)p;
                if(node.getParentProbability(par) > 0.5) {
                    node.setParentProbability(par, 1.0);
                } else {
                    node.setParentProbability(par, 1.0);
                }
            }
        }

        // computation
        StringBuilder sb = new StringBuilder("\n");
        for(List<MutableTuple<Tree,Double>> gtPerLocus : _geneTrees) {
            GeneTreeProbabilityDQ coal = new GeneTreeProbabilityDQ();
            List<Tree> input = new ArrayList<>();
            for(MutableTuple<Tree,Double> mt : gtPerLocus) {
                input.add(mt.Item1);
            }
            double xlg = 0.0;
            double prob = 0.0;
            MutableTuple<List<Double>, List<Double>> result = coal.calculateExpectXL(net, input, null);
            for(int i = 0; i < gtPerLocus.size(); i++) {
                xlg += gtPerLocus.get(i).Item2 * result.Item1.get(i);
                prob += gtPerLocus.get(i).Item2 * result.Item2.get(i);
            }
            sb.append(String.format("%.4f\n", xlg/prob));
        }
        return sb.toString();
    }

}