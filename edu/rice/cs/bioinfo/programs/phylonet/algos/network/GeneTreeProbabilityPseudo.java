package edu.rice.cs.bioinfo.programs.phylonet.algos.network;


import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 3/13/12
 * Time: 11:31 AM
 * To change this template use File | Settings | File Templates.
 */

public class GeneTreeProbabilityPseudo {
    private boolean _printDetails = false;
    private int _netNodeNum;
    private Set<NetNode> _articulationNodes;
    private Set<NetNode> _lowestArticulationNodes;
    private Map<String,Integer> _mergingMap;
    private Map<Integer,List<Tuple<Integer,Integer>>> _splittingMap;
    private Map<Integer, List<Tuple<Integer,Tuple<Integer,Integer>>>> _coalescingMap;
    private int _currentTripleID;
    private boolean _parallel;
    private int _batchSize = 10;
    //private List<MutableTuple<String, double[]>> _tripleFrequencies;

    /**
     * Constructor that initialize the variables.
     */
    public GeneTreeProbabilityPseudo(){
        _printDetails = false;
    }

    public void setBatchSize(int size){
        _batchSize = size;
    }

    public synchronized int getNextTripleID(){
        _currentTripleID = _currentTripleID + _batchSize;
        return _currentTripleID-_batchSize;

    }

    public void setParallel(boolean parallel){
        _parallel = parallel;
    }


    private void initializeMergingMap() {
        _mergingMap = new HashMap<>();
        /*
        _mergingMap.put("0|0",0);
        _mergingMap.put("0|1",1);
        _mergingMap.put("0|2",2);
        _mergingMap.put("0|3",3);
        _mergingMap.put("0|4",4);
        _mergingMap.put("0|5",5);
        _mergingMap.put("0|6",6);
        _mergingMap.put("0|7",7);
        _mergingMap.put("0|8",8);
        _mergingMap.put("0|9",9);
        _mergingMap.put("0|10",10);
        */
        _mergingMap.put("1|2",4);
        _mergingMap.put("1|3",5);
        _mergingMap.put("1|6",7);
        _mergingMap.put("2|3",6);
        _mergingMap.put("2|5",7);
        _mergingMap.put("3|4",7);
        _mergingMap.put("3|8",9);
    }

    private void initializeSplittingMap() {
        _splittingMap = new HashMap<>();
        _splittingMap.put(0,Arrays.asList(new Tuple<>(0,0)));
        _splittingMap.put(1,Arrays.asList(new Tuple<>(0,1)));
        _splittingMap.put(2,Arrays.asList(new Tuple<>(0,2)));
        _splittingMap.put(3,Arrays.asList(new Tuple<>(0,3)));
        _splittingMap.put(4,Arrays.asList(new Tuple<>(0,4),new Tuple<>(1,2)));
        _splittingMap.put(5,Arrays.asList(new Tuple<>(0,5),new Tuple<>(1,3)));
        _splittingMap.put(6,Arrays.asList(new Tuple<>(0,6),new Tuple<>(2,3)));
        _splittingMap.put(7,Arrays.asList(new Tuple<>(0,7),new Tuple<>(1,6),new Tuple<>(2,5),new Tuple<>(3,4)));
        _splittingMap.put(8,Arrays.asList(new Tuple<>(0,8)));
        _splittingMap.put(9,Arrays.asList(new Tuple<>(0,9),new Tuple<>(3,8)));
        _splittingMap.put(10,Arrays.asList(new Tuple<>(0,10)));
    }



    private void initializeCoalescingMap() {
        _coalescingMap = new HashMap<>();
        _coalescingMap.put(0,Arrays.asList(new Tuple<Integer, Tuple<Integer, Integer>>(0,new Tuple<>(0,0))));
        _coalescingMap.put(1,Arrays.asList(new Tuple<Integer, Tuple<Integer, Integer>>(1,new Tuple<>(1,1))));
        _coalescingMap.put(2,Arrays.asList(new Tuple<Integer, Tuple<Integer, Integer>>(2,new Tuple<>(1,1))));
        _coalescingMap.put(3,Arrays.asList(new Tuple<Integer, Tuple<Integer, Integer>>(3,new Tuple<>(1,1))));
        _coalescingMap.put(4,Arrays.asList(new Tuple<Integer, Tuple<Integer, Integer>>(4,new Tuple<>(2,2)),new Tuple<Integer, Tuple<Integer, Integer>>(8,new Tuple<>(2,1))));
        _coalescingMap.put(5,Arrays.asList(new Tuple<Integer, Tuple<Integer, Integer>>(5,new Tuple<>(2,2))));
        _coalescingMap.put(6,Arrays.asList(new Tuple<Integer, Tuple<Integer, Integer>>(6,new Tuple<>(2,2))));
        _coalescingMap.put(7,Arrays.asList(new Tuple<Integer, Tuple<Integer, Integer>>(7,new Tuple<>(3,3)),new Tuple<Integer, Tuple<Integer, Integer>>(9,new Tuple<>(3,2)),new Tuple<Integer, Tuple<Integer, Integer>>(10,new Tuple<>(3,1))));
        _coalescingMap.put(8,Arrays.asList(new Tuple<Integer, Tuple<Integer, Integer>>(8,new Tuple<>(1,1))));
        _coalescingMap.put(9,Arrays.asList(new Tuple<Integer, Tuple<Integer, Integer>>(9,new Tuple<>(2,2)),new Tuple<Integer, Tuple<Integer, Integer>>(10,new Tuple<>(2,1))));
        _coalescingMap.put(10,Arrays.asList(new Tuple<Integer, Tuple<Integer, Integer>>(10,new Tuple<>(1,1))));
    }


    public void initializeMaps(){
        initializeMergingMap();
        initializeSplittingMap();
        initializeCoalescingMap();
    }


    public void initialize(Network network){
        initializeMaps();
        computeArticulationNode(network);
        _netNodeNum = network.getReticulationCount();
    }



    public void computePseudoLikelihood(Network network, List<String> allTriplets, double[][] probs){
        if(_articulationNodes==null){
            initialize(network);
        }

        int tripleID = 0;
        if(_parallel){
            tripleID = getNextTripleID();
        }

        int totalTripleNum = allTriplets.size();
        int count = 0;
        while(tripleID < totalTripleNum){
            String[] triple = allTriplets.get(tripleID).split("&");
            //if(tripleID % 100 == 0)
            probs[tripleID][0] = calculateTripleProbability(network,triple);

            String temp = triple[1];
            triple[1] = triple[2];
            triple[2] = temp;
            probs[tripleID][1] = calculateTripleProbability(network, triple);

            probs[tripleID][2] = Math.max(1 - probs[tripleID][0] - probs[tripleID][1], 0);
/*
            temp = triple[0];
            triple[0] = triple[2];
            triple[2] = temp;
            probs[tripleID][2] = calculateTripleProbability(network,triple);


            double checkingOne = 0;


            for(int i=0; i<3; i++){
                checkingOne += probs[tripleID][i];

            }

            if(Math.abs(checkingOne-1)>0.0001){
                System.out.println(checkingOne);

            }

            for(double value: probs[tripleID]){
                System.out.print(value + ",");
            }
            */
/*
            temp = triple[0];
            triple[0] = triple[2];
            triple[2] = temp;
            probs[tripleID][2] = calculateTripleProbability(network,triple);
            double checkingOne = 0;


            for(int i=0; i<3; i++){
                checkingOne += probs[tripleID][i];

            }

            if(Math.abs(checkingOne-1)>0.0001){
                throw new RuntimeException(Arrays.toString(probs[tripleID]));
            }
*/

            if(_parallel){
                count++;
                if(count==_batchSize) {
                    tripleID = getNextTripleID();
                    count = 0;
                }
                else{
                    tripleID++;
                }
            }
            else{
                tripleID++;
            }

        }
    }




    private double calculateTripleProbability(Network network, String[] triple) {
        if (_printDetails) {
            System.out.println();
            System.out.println("Computing triple " + Arrays.toString(triple));
        }
        Map<Tuple<NetNode, NetNode>, Map<Integer, List<Configuration>>> edge2ACMinus = new HashMap<>();
        int netNodeID = 0;
        double totalProb = 0;
        for (Object o : Networks.postTraversal(network)) {
            NetNode node = (NetNode)o;
            if (_printDetails) {
                System.out.println();
                System.out.println("On node #" + node.getName());
            }
            Map<Integer, List<Configuration>> CACs = null;
            //set AC for a node
            if (node.isLeaf()) {
                Configuration config = new Configuration();
                String name = node.getName();
                int index = 0;
                for (; index < triple.length; index++) {
                    if (name.equals(triple[index])) {
                        break;
                    }
                }
                if(index<3){
                    CACs = new HashMap<>();
                    CACs.put(index + 1, Arrays.asList(config));
                }
                //t5 += (System.currentTimeMillis()-start)/1000.0;
            } else if (node.isNetworkNode()) {
                CACs = edge2ACMinus.get(new Tuple<>(node, (NetNode) (node.getChildren().iterator().next())));
            } else {
                Iterator children = node.getChildren().iterator();
                Map<Integer, List<Configuration>> AC1 = edge2ACMinus.get(new Tuple<>(node, (NetNode) (children.next())));
                Map<Integer, List<Configuration>> AC2 = edge2ACMinus.get(new Tuple<>(node, (NetNode) (children.next())));
                if (AC1 == null && AC2 != null || AC2 == null && AC1 != null) {
                    if(AC1 == null && AC2 != null) {
                        CACs = AC2;
                    }
                    else{
                        CACs = AC1;
                    }
                    if(_lowestArticulationNodes.contains(node)){
                        for(List<Configuration> configList: CACs.values()){
                            Configuration mergedConfig = null;
                            for(Configuration config: configList){
                                if(mergedConfig == null){
                                    mergedConfig = config;
                                }else {
                                    mergedConfig.addTotalProbability(config._totalProb);
                                }
                            }
                            mergedConfig.clearNetNodeChoice();
                            configList.clear();
                            configList.add(mergedConfig);
                        }
                    }

                } else if (AC1 != null && AC2 != null) {
                    CACs = new HashMap<>();
                    boolean isArticulation = _lowestArticulationNodes.contains(node);
                    //boolean isArticulation = false;
                    for (Map.Entry<Integer, List<Configuration>> entry1 : AC1.entrySet()) {
                        int configLineages1 = entry1.getKey();
                        for (Map.Entry<Integer, List<Configuration>> entry2 : AC2.entrySet()) {
                            int configLineages2 = entry2.getKey();
                            boolean canMerge = configLineages1 == 0 || configLineages2 == 0;
                            int targetConfigLineages = configLineages1;
                            if (targetConfigLineages == 0) {
                                targetConfigLineages = configLineages2;
                            }
                            if (!canMerge) {
                                String mergingType = configLineages1 < configLineages2 ? configLineages1 + "|" + configLineages2 : configLineages2 + "|" + configLineages1;
                                if (_mergingMap.containsKey(mergingType)) {
                                    targetConfigLineages = _mergingMap.get(mergingType);
                                    canMerge = true;

                                }
                            }
                            if (canMerge) {
                                List<Configuration> mergedConfigList = CACs.get(targetConfigLineages);
                                if (mergedConfigList == null) {
                                    mergedConfigList = new ArrayList<>();
                                    CACs.put(targetConfigLineages, mergedConfigList);

                                }
                                Configuration mergedConfig = null;

                                for (Configuration config1 : entry1.getValue()) {
                                    for (Configuration config2 : entry2.getValue()) {
                                        if (config1.isCompatible(config2)) {
                                            if (isArticulation) {
                                                if(mergedConfigList.size()==0){
                                                    mergedConfig = new Configuration(config1, config2);
                                                    mergedConfig.clearNetNodeChoice();
                                                    mergedConfigList.add(mergedConfig);
                                                }
                                                else {
                                                    if(mergedConfig==null){
                                                        mergedConfig = mergedConfigList.get(0);
                                                    }
                                                    double newProb = Math.max(0, config1._totalProb * config2._totalProb);
                                                    mergedConfig.addTotalProbability(newProb);
                                                }
                                            } else {
                                                mergedConfigList.add(new Configuration(config1, config2));
                                            }
                                        }
                                    }
                                }
                                if(mergedConfigList.size()==0){
                                    CACs.remove(targetConfigLineages);
                                }
                            }
                        }
                    }
                }
            }

            if (_printDetails) {
                System.out.print("AC: {");
                if(CACs!=null) {
                    for (Map.Entry<Integer, List<Configuration>> entry : CACs.entrySet()) {
                        //System.out.print(entry.getKey() + ": ");
                        if(entry.getValue().size()==0){
                            throw new RuntimeException("Wrong: " + entry.getKey());
                        }
                        for (Configuration config : entry.getValue()) {
                            System.out.print(config.toString(triple, entry.getKey())+",");
                        }

                    }
                }
                System.out.println("}");
            }



            //set AC- for a node
            if(CACs!=null) {
                if (CACs.containsKey(7) && _articulationNodes.contains(node)) {
                    totalProb = CACs.get(7).get(0)._totalProb/3;
                    if(CACs.containsKey(9)){
                        totalProb += CACs.get(9).get(0)._totalProb;
                    }
                    if(CACs.containsKey(10)){
                        totalProb += CACs.get(10).get(0)._totalProb;
                    }

                    break;
                }


                if (node.isNetworkNode()) {
                    Map<Integer, List<Configuration>> ACplus1 = new HashMap<>();
                    Map<Integer, List<Configuration>> ACplus2 = new HashMap<>();
                    splittingAtNetworkNode(CACs, netNodeID, ACplus1, ACplus2);
                    netNodeID++;
                    int index = 0;
                    for (Object parentO : node.getParents()) {
                        NetNode parent = (NetNode) parentO;
                        Map<Integer, List<Configuration>> ACMinus = new HashMap<>();
                        if (index == 0) {
                            computeACMinus(ACplus1, ACMinus, node.getParentDistance(parent), node.getParentProbability(parent));

                        } else {
                            computeACMinus(ACplus2, ACMinus, node.getParentDistance(parent), node.getParentProbability(parent));

                        }
                        edge2ACMinus.put(new Tuple<>(parent, node), ACMinus);
                        index++;
                        if (_printDetails) {
                            System.out.print("ACMinus: {");
                            for (Map.Entry<Integer, List<Configuration>> entry : ACMinus.entrySet()) {
                                //System.out.print(entry.getKey() + ": ");
                                for (Configuration config : entry.getValue()) {
                                    System.out.print(config.toString(triple, entry.getKey())+", ");
                                }

                            }

                            System.out.println("}");
                        }
                    }
                } else {
                    NetNode parent = (NetNode) node.getParents().iterator().next();
                    Map<Integer, List<Configuration>> ACMinus = new HashMap<>();
                    computeACMinus(CACs, ACMinus, node.getParentDistance(parent), 1);
                    edge2ACMinus.put(new Tuple<>(parent, node), ACMinus);

                    if (_printDetails) {
                        System.out.print("ACMinus: {");
                        for (Map.Entry<Integer, List<Configuration>> entry : ACMinus.entrySet()) {
                            //System.out.print(entry.getKey() + ": ");
                            for (Configuration config : entry.getValue()) {
                                System.out.print(config.toString(triple, entry.getKey())+",");
                            }

                        }

                        System.out.println("}");
                    }
                }


            }

        }
        //System.out.println(t1 + " " + t2 + " " + t3+ " " + t4+ " " + t5+ " " + t6);
        return totalProb;

    }

    private void computeACMinus(Map<Integer, List<Configuration>> CACs, Map<Integer, List<Configuration>> ACMinus, double bl, double inheritanceProb){

        Map<Integer, Map<Configuration,Configuration>> ACMinusMap = new HashMap<>();
        for(Map.Entry<Integer, List<Configuration>> entry: CACs.entrySet()){
            for(Tuple<Integer,Tuple<Integer,Integer>> coalesceInto: _coalescingMap.get(entry.getKey())){
                Tuple<Integer,Integer> ij = coalesceInto.Item2;

                double prob = 1;
                if(inheritanceProb!=1){
                    prob = Math.pow(inheritanceProb, ij.Item1);
                }
                if(ij.Item1>1){
                    prob *= gij(bl, ij.Item1, ij.Item2);
                    if(ij.Item1==3 && ij.Item2!=3){
                        prob = prob/3;
                    }
                }
                if(prob!=0) {
                    Map<Configuration,Configuration> configurationMap = ACMinusMap.get(coalesceInto.Item1);
                    if(configurationMap==null){
                        configurationMap = new HashMap<>();
                        ACMinusMap.put(coalesceInto.Item1, configurationMap);
                    }
                    for (Configuration config : entry.getValue()) {
                        Configuration copy = new Configuration(config);
                        double newProb = Math.max(0, copy._totalProb * prob);
                        Configuration exisitingConfig = configurationMap.get(copy);
                        if (exisitingConfig == null) {
                            copy.setTotalProbability(Math.max(0, newProb));
                            configurationMap.put(copy, copy);
                        } else {
                            exisitingConfig.addTotalProbability(newProb);
                        }

                    }
                }
            }
        }

        for(Map.Entry<Integer,Map<Configuration,Configuration>> entry: ACMinusMap.entrySet()){
            List<Configuration> configList = new ArrayList<>();
            configList.addAll(entry.getValue().keySet());
            ACMinus.put(entry.getKey(), configList);
        }
    }


    private void splittingAtNetworkNode(Map<Integer, List<Configuration>> CACs, int netNodeID, Map<Integer,List<Configuration>> ACplus1, Map<Integer,List<Configuration>> ACplus2){
        int netIndex = 1;
        for(Map.Entry<Integer,List<Configuration>> entry: CACs.entrySet()){
            for(Tuple<Integer,Integer> splitting: _splittingMap.get(entry.getKey())){
                List<Configuration> list1 = ACplus1.get(splitting.Item1);
                if(list1==null){
                    list1 = new ArrayList<>();
                    ACplus1.put(splitting.Item1, list1);
                }
                List<Configuration> list2 = ACplus2.get(splitting.Item2);
                if(list2==null){
                    list2 = new ArrayList<>();
                    ACplus2.put(splitting.Item2, list2);
                }
                List<Configuration> opList1 = null;
                List<Configuration> opList2 = null;
                if(entry.getKey()!=0) {
                    opList1 = ACplus1.get(splitting.Item2);
                    if (opList1 == null) {
                        opList1 = new ArrayList<>();
                        ACplus1.put(splitting.Item2, opList1);
                    }
                    opList2 = ACplus2.get(splitting.Item1);
                    if (opList2 == null) {
                        opList2 = new ArrayList<>();
                        ACplus2.put(splitting.Item1, opList2);
                    }
                }

                for(Configuration configToSplit: entry.getValue()){
                    double prob = Math.sqrt(configToSplit._totalProb);
                    Configuration newConfig1 = new Configuration(configToSplit);
                    newConfig1.addNetNodeChoice(netNodeID, netIndex);
                    newConfig1.setTotalProbability(prob);
                    Configuration newConfig2 = new Configuration(newConfig1);
                    newConfig2.setTotalProbability(prob);
                    list1.add(newConfig1);
                    list2.add(newConfig2);
                    netIndex++;
                    if(entry.getKey()!=0) {
                        newConfig1 = new Configuration(configToSplit);
                        newConfig1.addNetNodeChoice(netNodeID, netIndex);
                        newConfig1.setTotalProbability(prob);
                        newConfig2 = new Configuration(newConfig1);
                        newConfig2.setTotalProbability(prob);
                        opList1.add(newConfig1);
                        opList2.add(newConfig2);
                        netIndex++;
                    }
                }
            }
        }
    }




    /**
     * The function is to calculate the g_{ij} function.
     * @param	length	the branch length
     * @param 	i	the number of lineages in
     * @param	j	the number of lineages out
     * @return	the resulting probability
     */
    private double gij(double length, int i, int j){
        if(length == TNode.NO_DISTANCE || length == -1){
            if(j == 1){
                return 1;
            }else{
                return 0;
            }
        }
        if(length==0){
            if(i==j)
                return 1;
            else
                return 0;
        }
        if(i==0){
            return 1;
        }

        double result = 0;
        for(int k=j; k<=i; k++){

            double temp = Math.exp(0.5 * k * (1.0 - k) * length) * (2.0 * k - 1.0) * Math.pow(-1, k - j)*(fact(j, j + k - 2))*(fact(i - k + 1, i));
            double denom = fact(1, j)*(fact(1, k - j))*(fact(i, i + k - 1));
            result += temp/denom;
        }
        return result;
    }


    private double chooseD(int N, int K) {
        double ret = 1.0;
        for (int k = 0; k < K; k++) {
            //ret = ret.multiply(BigInteger.valueOf(N-k))
            //.divide(BigInteger.valueOf(k+1));
            ret = ret*((N-k+0.0)/(k+1));
        }

        return ret;
    }


    private double fact(int start, int end){
        double result = 1;
        for(int i=start; i<=end; i++){
            result = result*i;
        }

        return result;
    }


    private void computeArticulationNode(Network net){
        _articulationNodes = new HashSet<>();
        _lowestArticulationNodes = new HashSet<NetNode>();
        for(Object o: Networks.postTraversal(net)){
            NetNode node = (NetNode)o;
            if(node.isLeaf()){
                _articulationNodes.add(node);
            }

            else if(node.isRoot()){
                boolean fArticulate = true;
                for(Object childO: node.getChildren()){
                    if(!_articulationNodes.contains(childO)){
                        fArticulate = false;
                        break;
                    }
                }
                if(!fArticulate){
                    _lowestArticulationNodes.add(node);
                }
                _articulationNodes.add(node);

            }
            else if(node.isTreeNode()){
                boolean ftotal = true;
                for(Object childO: node.getChildren()){
                    if(!_articulationNodes.contains(childO)){
                        ftotal = false;
                        break;
                    }
                }
                if(ftotal){
                    _articulationNodes.add(node);
                }else{
                    NetNode parent = (NetNode)node.getParents().iterator().next();
                    double distance = node.getParentDistance(parent);
                    parent.removeChild(node);
                    boolean disconnect = isValidNetwork(net, parent);
                    parent.adoptChild(node, distance);
                    if (disconnect) {
                        _lowestArticulationNodes.add(node);
                        _articulationNodes.add(node);
                    }
                }

            }
        }
    }


    private static  boolean isValidNetwork(Network net, NetNode ignoreNode){
        Set<NetNode> visited = new HashSet<NetNode>();
        Set<NetNode> seen = new HashSet<NetNode>();
        for(Object o: net.bfs()){
            NetNode node = (NetNode)o;
            if(node.getIndeg()==1 && node.getOutdeg()==1 && node!=ignoreNode){
                return false;
            }
            visited.add(node);
            for(Object parent: node.getParents()){
                seen.add((NetNode)parent);
            }
            for(Object child: node.getChildren()){
                seen.add((NetNode)child);
            }
        }
        return visited.size()==seen.size();
    }


    private class Configuration{
        private double _totalProb;
        int[] _netNodeIndex;


        public Configuration(){
            _netNodeIndex = new int[_netNodeNum];
            Arrays.fill(_netNodeIndex, 0);
            _totalProb = 1;
        }

        public Configuration(Configuration config){
            _totalProb = config._totalProb;
            _netNodeIndex = config._netNodeIndex.clone();
        }

        public Configuration(Configuration config1, Configuration config2){
            _totalProb = Math.max(0,config1._totalProb*config2._totalProb);
            _netNodeIndex = new int[_netNodeNum];
            for(int i=0; i< _netNodeNum; i++){
                if(config1._netNodeIndex[i] == config2._netNodeIndex[i]){
                    _netNodeIndex[i] = config1._netNodeIndex[i];
                }
                else{
                    _netNodeIndex[i] = Math.max(config1._netNodeIndex[i], config2._netNodeIndex[i]);
                }
            }
        }

        public boolean isCompatible(Configuration config){
            boolean compatible = true;
            for(int i=0; i< _netNodeNum; i++){
                if(_netNodeIndex[i] != config._netNodeIndex[i] && _netNodeIndex[i]!=0 && config._netNodeIndex[i]!=0){
                    compatible = false;
                    break;
                }
            }
            return compatible;
        }


        public void setTotalProbability(double prob){
            _totalProb = prob;
        }



        public void addTotalProbability(double adds){
            _totalProb += adds;
        }



        public String toString(String[] triple, int id){
            String exp = "";
            switch(id){
                case 0: exp = "{}";
                    break;
                case 1: exp = "{"+triple[0]+"}";
                    break;
                case 2: exp = "{"+triple[1]+"}";
                    break;
                case 3: exp = "{"+triple[2]+"}";
                    break;
                case 4: exp = "{"+triple[0]+"}{"+triple[1]+"}";
                    break;
                case 5: exp = "{"+triple[0]+"}{"+triple[2]+"}";
                    break;
                case 6: exp = "{"+triple[1]+"}{"+triple[2]+"}";
                    break;
                case 7: exp = "{"+triple[0]+"}{"+triple[1]+"}{"+triple[2]+"}";
                    break;
                case 8: exp = "{"+triple[0]+","+triple[1]+"}";
                    break;
                case 9: exp = "{"+triple[0]+","+triple[1]+"}{"+triple[2]+"}";
                    break;
                case 10: exp = "{"+triple[0]+","+triple[1]+","+triple[2]+"}";
            }

            exp += "[";
            for(int i=0; i<_netNodeIndex.length; i++){
                exp = exp + _netNodeIndex[i];
                if(i!=_netNodeIndex.length-1){
                    exp = exp + ",";
                }
            }
            exp = exp + "]:" + _totalProb;
            return exp;
        }




        public void addNetNodeChoice(int net, int index){
            _netNodeIndex[net] = index;
        }

        public void setNetNodeChoice(int[] choice){
            _netNodeIndex = choice.clone();
        }

        public void clearNetNodeChoice(){
            Arrays.fill(_netNodeIndex, 0);
        }


        public boolean equals(Object o) {
            if(!(o instanceof Configuration)){
                return false;
            }

            Configuration config = (Configuration) o;
            return Arrays.equals(config._netNodeIndex, _netNodeIndex);
        }

        public int hashCode(){
            return Arrays.hashCode(_netNodeIndex);
        }

    }




}
