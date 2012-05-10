package edu.rice.cs.bioinfo.programs.phylonet.algos.network;


import java.util.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 5/9/12
 * Time: 1:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class MDCOnNetworkYF {
    String[] _netTaxa;
    Map<NetNode, Integer> _netNode2id;
    BitSet _totalCoverNode;
    boolean _printDetail = false;
    int _netNodeNum;
    int _totalNodeNum;
    int[][] _totalNetNodeLinNum;



    public List<Integer> countExtraCoal(Network<List<Configuration>> network, List<Tree> gts, Map<String, List<String>> species2alleles){
        List<Integer> xlList = new ArrayList<Integer>();
        processNetwork(network);
        //System.out.println(gts.size());
        //System.exit(0);
        //int pos = 0;

        for(Tree gt: gts){
            Trees.removeBinaryNodes((MutableTree)gt);
            List<STITreeCluster> gtClusters = new ArrayList<STITreeCluster>();
            String[] gtTaxa = gt.getLeaves();
            int[][] gtclConstitution = new int[gt.getNodeCount()-gtTaxa.length][2];
            processGT(gt,gtTaxa, gtclConstitution, gtClusters);

            if(!checkLeafAgreement(species2alleles, gtTaxa)){
                throw new RuntimeException("Gene tree " + gt + " has leaf that the network doesn't have.");
            }


            HashSet<String> gtTaxaSet = new HashSet<String>();
            for(String taxon: gtTaxa){
                gtTaxaSet.add(taxon);
            }

            cleanNetwork(network);
            HashMap<BitSet, List<Configuration>> edge2ACminus = new HashMap<BitSet, List<Configuration>>();

            int netNodeIndex = 0;
            int numLineages = gtClusters.size();
            int xl = Integer.MAX_VALUE;
            for(NetNode<List<Configuration>> node: walkNetwork(network)){
                //System.out.println(edge2ACminus);
                if(_printDetail){
                    System.out.println();
                    System.out.println("On node #" + _netNode2id.get(node) + " " + node.getName());
                }
                List<Configuration> CACs = node.getData();
                if(CACs == null){
                    CACs = new ArrayList<Configuration>();
                    node.setData(CACs);
                }

                //set AC for a node
                if(node.isLeaf()){

                    Configuration config = new Configuration(gtTaxa, numLineages);
                    if(species2alleles == null){
                        if(gtTaxaSet.contains(node.getName())){
                            STITreeCluster cl = new STITreeCluster(gtTaxa);
                            cl.addLeaf(node.getName());
                            config.addLineage(cl, gtClusters.indexOf(cl));
                        }
                    }
                    else{
                        for(String allele: species2alleles.get(node.getName())){
                            if(gtTaxaSet.contains(allele)){
                                STITreeCluster cl = new STITreeCluster(gtTaxa);
                                cl.addLeaf(allele);
                                config.addLineage(cl, gtClusters.indexOf(cl));
                            }
                        }
                    }
                    config.setExtraLineage(0);
                    CACs.add(config);
                }
                else{
                    if(node.getOutdeg() == 1){
                        Iterator<NetNode<List<Configuration>>> childNode = node.getChildren().iterator();
                        BitSet edge = new BitSet();
                        edge.set(_netNode2id.get(node));
                        edge.set(_netNode2id.get(childNode.next()));
                        CACs.addAll(edge2ACminus.get(edge));
                    }
                    else{
                        //TODO only on binary nodes

                        Iterator<NetNode<List<Configuration>>> childNode = node.getChildren().iterator();
                        BitSet edge1 = new BitSet();
                        edge1.set(_netNode2id.get(node));
                        edge1.set(_netNode2id.get(childNode.next()));
                        List<Configuration> AC1 = edge2ACminus.get(edge1);
                        BitSet edge2 = new BitSet();
                        edge2.set(_netNode2id.get(node));
                        edge2.set(_netNode2id.get(childNode.next()));
                        List<Configuration> AC2 = edge2ACminus.get(edge2);
                        int minimum = Integer.MAX_VALUE;
                        boolean totalCover = _totalCoverNode.get(_netNode2id.get(node));
                        BitSet[][] netNodeLineages = null;
                        if(totalCover){
                            netNodeLineages = new BitSet[_netNodeNum][2];
                            for(int i=0; i< _netNodeNum; i++){
                                for(int j=0; j<2; j++){
                                    netNodeLineages[i][j] = new BitSet(numLineages);
                                }
                            }
                        }
                        if(_printDetail){
                            System.out.print("AC1: {");
                            for(Configuration config: AC1){
                                System.out.print(config.toString(gtClusters)+"  ");
                            };
                            System.out.println("}");
                            System.out.print("AC2: {");
                            for(Configuration config: AC2){
                                System.out.print(config.toString(gtClusters)+"  ");
                            };
                            System.out.println("}");
                        }
                        for(Configuration config1: AC1){
                            for(Configuration config2: AC2){
                                //System.out.println("here2");
                                if(config1.isCompatible(config2)){
                                    Configuration mergedConfig = new Configuration(config1, config2, numLineages);

                                    if(totalCover){

                                        if(minimum > mergedConfig._xl){
                                            minimum = mergedConfig._xl;
                                            CACs.clear();
                                            CACs.add(mergedConfig);
                                            for(int i=0; i< _netNodeNum; i++){
                                                for(int j=0; j<2; j++){
                                                    netNodeLineages[i][j] = (BitSet)(mergedConfig._netNodeLineages[i][j].clone());
                                                }
                                            }
                                        }
                                        else if(minimum == mergedConfig._xl){
                                            for(int i=0; i< _netNodeNum; i++){
                                                for(int j=0; j<2; j++){
                                                    netNodeLineages[i][j].and(mergedConfig._netNodeLineages[i][j]);
                                                    //System.out.print(netNodeLineages[i][0].cardinality() + "/" + netNodeLineages[i][1].cardinality() + "   ");
                                                }
                                            }
                                        }

                                        /*
                                                  else if(minimum == mergedConfig._xl){
                                                      CACs.get(0).addNetNodeLineageNum(mergedConfig);
                                                  }
                                                  */
                                    }
                                    else{

                                        CACs.add(mergedConfig);
                                        /*
                                                  int index = CACs.indexOf(mergedConfig);
                                                  if(index == -1){
                                                      CACs.add(mergedConfig);
                                                  }
                                                  else{

                                                      Configuration exist = CACs.get(index);
                                                      exist.setExtraLineage(Math.min(mergedConfig._xl, exist._xl));
                                                  }
                                                  */
                                    }
                                }
                                
                            }
                        }
                        if(totalCover){
                            //System.out.println(CACs.size());
                            CACs.get(0).setNetNodeLineageNum(netNodeLineages);
                        }
                    }
                }
                if(_printDetail){
                    System.out.print("AC: {");
                    for(Configuration config: CACs){
                        System.out.print(config.toString(gtClusters)+"  ");
                    };
                    System.out.println("}");
                }
                //set AC- for a node
                if(node.isRoot()){

                    if(CACs.size()!=1){
                        System.err.println("Error");
                    }
                    Configuration optimalConfig = CACs.get(0);
                    xl = optimalConfig._xl;
                    xlList.add(xl);

                    for(int i=0; i< _netNodeNum; i++){
                        //System.out.print(optimalConfig._netNodeLineages[i][0].cardinality() + "/" + optimalConfig._netNodeLineages[i][1].cardinality() + "   ");
                        _totalNetNodeLinNum[i][0] += optimalConfig._netNodeLineages[i][0].cardinality();
                        _totalNetNodeLinNum[i][1] += optimalConfig._netNodeLineages[i][1].cardinality();
                    }
                    //System.out.println();
                }
                else if(node.isTreeNode()){
                    List<Configuration> ACminus = new ArrayList<Configuration>();
                    for(Configuration config: CACs){
                        for(int i=0; i< gtclConstitution.length; i++){
                            config.mergeCluster(i, gtclConstitution[i]);
                        }
                        config.addExtraLineage(Math.max(0, config.getLineageCount()-1));
                        ACminus.add(config);
                        /*
                              int index = ACminus.indexOf(config);
                              if(index == -1){
                                  ACminus.add(config);
                              }
                              else{
                                  //System.out.println(mergedConfig);
                                  System.out.println("here");
                                  Configuration exist = ACminus.get(index);
                                  exist.setExtraLineage(Math.min(config._xl,exist._xl));
                              }
                              */
                    }

                    BitSet newEdge = new BitSet();
                    newEdge.set(_netNode2id.get(node));
                    newEdge.set(_netNode2id.get(node.getParents().iterator().next()));
                    edge2ACminus.put(newEdge, ACminus);

                    if(_printDetail){
                        System.out.print("ACminus: {");
                        for(Configuration config: ACminus){
                            System.out.print(config.toString(gtClusters)+"  ");
                        };
                        System.out.println("}");
                    }
                }
                else {
                    List<Configuration> ACminus1 = new ArrayList<Configuration>();
                    List<Configuration> ACminus2 = new ArrayList<Configuration>();
                    int configIndex = 1;
                    for(Configuration config: CACs){
                        int numLineage = config.getLineageCount();
                        for(int i=0; i<=numLineage; i++){
                            for(BitSet selectedLineages: getSelected(numLineage,i)){
                                Configuration newConfig1 = new Configuration(gtTaxa, numLineages);
                                Configuration newConfig2 = new Configuration(gtTaxa, numLineages);
                                newConfig1.setNetNodeLineageNum(config._netNodeLineages);
                                newConfig2.setNetNodeLineageNum(config._netNodeLineages);

                                int index = 0;
                                for (int k = config._lineages.nextSetBit(0); k >= 0; k = config._lineages.nextSetBit(k+1)) {
                                    if(selectedLineages.get(index)){
                                        newConfig1.addLineage(gtClusters.get(k), k);
                                        newConfig1.addNetNodeLineageNum(netNodeIndex, 0, k);
                                    }
                                    else{
                                        newConfig2.addLineage(gtClusters.get(k), k);
                                        newConfig2.addNetNodeLineageNum(netNodeIndex, 1, k);
                                    }
                                    index ++;

                                }

                                newConfig1.setNetNodeChoice(config._netNodeIndex);
                                //newConfig1.setNetNodeLineageNum(config._netNodeLineages);
                                newConfig1.setExtraLineage(config._xl);
                                newConfig1.addNetNodeChoice(netNodeIndex, configIndex);
                                //newConfig1.addNetNodeLineageNum(netNodeIndex, i, 0);
                                ACminus1.add(newConfig1);

                                newConfig2.setNetNodeChoice(config._netNodeIndex);
                                //newConfig2.setNetNodeLineageNum(config._netNodeLinNum);
                                newConfig2.setExtraLineage(0);
                                newConfig2.addNetNodeChoice(netNodeIndex, configIndex);
                                //newConfig2.addNetNodeLineageNum(netNodeIndex, 0, numLineage-i);
                                ACminus2.add(newConfig2);
                                configIndex ++;
                            }
                        }
                    }
                    /*
                         System.out.println(ACminus1);
                         System.out.println(ACminus2);
                         System.out.println();
                         */
                    Iterator<NetNode<List<Configuration>>> it = node.getParents().iterator();
                    for(int i=0; i<2; i++){
                        List<Configuration> ACminus;
                        if(i==0){
                            ACminus = ACminus1;
                        }
                        else{
                            ACminus = ACminus2;
                        }
                        NetNode<List<Configuration>> parentNode = it.next();
                        if(node.getParentDistance(parentNode) != 0){
                            for(Configuration config: ACminus){
                                config.addExtraLineage(Math.max(0, config.getLineageCount()-1));
                            }
                        }
                        BitSet newEdge = new BitSet();
                        newEdge.set(_netNode2id.get(node));
                        newEdge.set(_netNode2id.get(parentNode));
                        edge2ACminus.put(newEdge, ACminus);

                        if(_printDetail){
                            System.out.print("ACminus to " + parentNode.getName()+ ": {");
                            for(Configuration config: ACminus){
                                System.out.print(config.toString(gtClusters)+"  ");
                            };
                            System.out.println("}");

                        }
                    }
                    netNodeIndex ++;
                }

            }

        }
        
        return xlList;
    }

    public void setPrint(boolean toPrint){
        _printDetail = toPrint;
    }

    private boolean checkLeafAgreement(Map<String, List<String>> species2alleles, String[] gtTaxa){
        if(species2alleles==null){
            for(String alleleG: gtTaxa){
                boolean found = false;
                for(String alleleS: _netTaxa){
                    if(alleleG.equals(alleleS)){
                        found = true;
                        break;
                    }
                }
                if(!found){
                    return false;
                }
            }
        }
        else{
            //HashSet<String> speciesSet = new HashSet<String>();
            for(String alleleG: gtTaxa){
                boolean found = false;
                for(Map.Entry<String, List<String>> entry: species2alleles.entrySet()){
                    for(String alleleS: entry.getValue()){
                        if(alleleS.equals(alleleG)){
                            found = true;
                            //speciesSet.add(entry.getKey());
                            break;
                        }
                    }
                    if(found){
                        break;
                    }
                }
                if(!found){
                    return false;
                }
            }
            /*
               if(speciesSet.size()!=species2alleles.size()){
                   return false;
               }
               */
        }
        return true;
    }


    private void processNetwork(Network<List<Configuration>> net){
        removeBinaryNodes(net);
        _netNodeNum = 0;
        _totalNodeNum = 0;
        _netNode2id = new HashMap<NetNode, Integer>();
        List<String> taxa = new ArrayList<String>();
        for(NetNode<List<Configuration>> node: net.dfs()){
            _netNode2id.put(node, _totalNodeNum++);
            if(node.isLeaf()){
                taxa.add(node.getName());
            }else if(node.isNetworkNode()){
                _netNodeNum++;
            }
        }
        _netTaxa = (String[]) taxa.toArray(new String[0]);
        _totalNetNodeLinNum = new int[_netNodeNum][2];
        computeNodeCoverage(net);
    }


    private void processGT(Tree gt, String[] gtTaxa, int[][] gtclConstitution,List<STITreeCluster> gtClusters){
        Map<TNode, STITreeCluster> map = new HashMap<TNode, STITreeCluster>();
        for (TNode node : gt.postTraverse()) {
            STITreeCluster cl = new STITreeCluster(gtTaxa);
            if (node.isLeaf()) {
                cl.addLeaf(node.getName());
                //_gtClusters.add(cl);
            }
            else {
                for(TNode child : node.getChildren()) {
                    cl = cl.merge(map.get(child));
                }

                int i = 0;
                for(; i< gtClusters.size(); i++){
                    if(gtClusters.get(i).getClusterSize() > cl.getClusterSize()){
                        break;
                    }
                }
                gtClusters.add(i, cl);
            }
            map.put(node, cl);
        }

        for(int i=0; i< gtTaxa.length; i++){
            STITreeCluster cl = new STITreeCluster(gtTaxa);
            cl.addLeaf(gtTaxa[i]);
            gtClusters.add(cl);
        }

        for (Map.Entry<TNode, STITreeCluster> entry: map.entrySet()) {
            TNode pNode = entry.getKey();
            if(!entry.getKey().isLeaf()){
                int parent = gtClusters.indexOf(map.get(pNode));
                int index = 0;
                for(TNode child : pNode.getChildren()) {
                    gtclConstitution[parent][index++] = gtClusters.indexOf(map.get(child));
                }
            }
        }

    }


    private void removeBinaryNodes(Network<List<Configuration>> net)
    {
        // Find all binary nodes.
        List<NetNode<List<Configuration>>> binaryNodes = new LinkedList<NetNode<List<Configuration>>>();
        for (NetNode<List<Configuration>> node : net.bfs()) {
            if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                binaryNodes.add(node);
            }
        }

        // Remove them.
        for (NetNode<List<Configuration>> node : binaryNodes) {
            NetNode<List<Configuration>> child = node.getChildren().iterator().next();	// Node's only child.
            if(child.getIndeg() != 1){
                continue;
            }
            NetNode<List<Configuration>> parent = node.getParents().iterator().next();	// Node's only parent.
            double distance = node.getParentDistance(parent) + child.getParentDistance(node);
            double gamma = node.getParentProbability(parent) * child.getParentProbability(node);
            parent.removeChild(node);
            node.removeChild(child);
            parent.adoptChild(child, distance);
            child.setParentProbability(parent, gamma);
        }
    }


    private List<NetNode> walkNetwork(Network net){
        Stack<NetNode> stack = new Stack<NetNode>();
        List<NetNode> searchedNodes = new ArrayList<NetNode>();
        stack.push(net.getRoot());
        Map<NetNode, Integer> node2index = new HashMap<NetNode, Integer>();
        node2index.put(net.getRoot(), 0);

        while(!stack.isEmpty()){

            NetNode topNode = stack.peek();
            int index = node2index.get(topNode);
            if(index == topNode.getOutdeg()){
                searchedNodes.add(stack.pop());
            }
            else{
                Iterator<NetNode> it = topNode.getChildren().iterator();
                for(int i=0; i<index; i++){
                    it.next();
                }
                NetNode<List<Configuration>> child = it.next();
                if(searchedNodes.contains(child)){
                    node2index.put(topNode, index + 1);
                }
                else{
                    stack.push(child);
                    node2index.put(child, 0);
                }
            }
        }

        return searchedNodes;
    }


    private void cleanNetwork(Network<List<Configuration>> net){
        for(NetNode<List<Configuration>> node: net.bfs()){
            node.setData(null);
        }
    }

    private void computeNodeCoverage(Network<List<Configuration>> net){
        _totalCoverNode = new BitSet(_totalNodeNum);
        for(NetNode<List<Configuration>> trNode: net.getTreeNodes()){
            if(trNode.isRoot()){
                //System.out.println(_netNode2id.get(trNode)==null);
                _totalCoverNode.set(_netNode2id.get(trNode), true);
            }
            else if(!trNode.isLeaf()){
                NetNode parent = trNode.getParents().iterator().next();
                double distance = trNode.getParentDistance(parent);
                parent.removeChild(trNode);
                boolean disconnect = isValidNetwork(net);
                parent.adoptChild(trNode, distance);
                if(disconnect){
                    //System.out.println(_netNode2id.get(trNode)==null);
                    _totalCoverNode.set(_netNode2id.get(trNode), true);
                }
            }
        }
    }

    private List<BitSet> getSelected(int n, int m){
        List<BitSet> selectedList = new ArrayList<BitSet>();
        int[] order = new int[m+1];
        for(int i=0; i<=m; i++){
            order[i] = i-1;
        }
        int k = m;
        boolean flag = true;
        while(order[0] == -1){
            if(flag){
                BitSet bs = new BitSet(n);
                for(int i=1; i<=m; i++){
                    bs.set(order[i]);
                }
                selectedList.add(bs);
                flag = false;
            }

            order[k]++;
            if(order[k] == n){
                order[k--] = 0;
                continue;
            }

            if(k < m){
                order[++k] = order[k-1];
                continue;
            }

            if(k == m)
                flag = true;
        }

        return selectedList;
    }


    private boolean isValidNetwork(Network<List<Configuration>> net){
        BitSet visited = new BitSet();
        BitSet seen = new BitSet();
        for(NetNode<List<Configuration>> node: net.bfs()){
            visited.set(_netNode2id.get(node), true);
            for(NetNode parent: node.getParents()){
                seen.set(_netNode2id.get(parent), true);
            }
            for(NetNode child: node.getChildren()){
                seen.set(_netNode2id.get(child), true);
            }
        }
        return visited.cardinality()==seen.cardinality();
    }

    private class Configuration{
        private STITreeCluster _coverage;
        private BitSet _lineages;
        private int _xl;
        int[] _netNodeIndex;
        //int[][] _netNodeLinNum;
        BitSet[][] _netNodeLineages;

        public Configuration(String[] gtTaxa, int numLineages){
            _lineages = new BitSet(numLineages);
            _netNodeIndex = new int[_netNodeNum];
            //_netNodeLinNum = new int[_netNodeNum][2];
            _netNodeLineages = new BitSet[_netNodeNum][2];
            for(int i=0; i< _netNodeNum; i++){
                for(int j=0; j<2; j++){
                    _netNodeLineages[i][j] = new BitSet(numLineages);
                }
            }
            _coverage = new STITreeCluster(gtTaxa);
            Arrays.fill(_netNodeIndex, 0);
        }


        public Configuration(Configuration config1, Configuration config2, int numLineages){
            _lineages = (BitSet)config1._lineages.clone();
            _lineages.or(config2._lineages);
            _xl = config1._xl + config2._xl;
            _coverage = new STITreeCluster(config1._coverage);
            _coverage = _coverage.merge(config2._coverage);
            _netNodeIndex = new int[_netNodeNum];
            //_netNodeLinNum = new int[_netNodeNum][2];
            _netNodeLineages = new BitSet[_netNodeNum][2];
            for(int i=0; i< _netNodeNum; i++){
                for(int j=0; j<2; j++){
                    _netNodeLineages[i][j] = new BitSet(numLineages);
                }
            }
            for(int i=0; i< _netNodeNum; i++){
                if(config1._netNodeIndex[i] == config2._netNodeIndex[i]){
                    _netNodeIndex[i] = config1._netNodeIndex[i];
                }
                else{
                    _netNodeIndex[i] = Math.max(config1._netNodeIndex[i], config2._netNodeIndex[i]);
                }
                //_netNodeLinNum[i][0] = config1._netNodeLinNum[i][0] + config2._netNodeLinNum[i][0];
                //_netNodeLinNum[i][1] = config1._netNodeLinNum[i][1] + config2._netNodeLinNum[i][1];
                _netNodeLineages[i][0] = (BitSet)(config1._netNodeLineages[i][0].clone());
                _netNodeLineages[i][0].or(config2._netNodeLineages[i][0]);
                _netNodeLineages[i][1] = (BitSet)(config1._netNodeLineages[i][1].clone());
                _netNodeLineages[i][1].or(config2._netNodeLineages[i][1]);
            }
        }


        public boolean isCompatible(Configuration config){
            boolean compatible = true;
            if(!_coverage.isDisjoint(config._coverage)){
                compatible = false;
            }
            for(int i=0; i< _netNodeNum; i++){
                if(_netNodeIndex[i] != config._netNodeIndex[i] && _netNodeIndex[i]!=0 && config._netNodeIndex[i]!=0){
                    compatible = false;
                }
            }
            return compatible;
        }

        public void addLineage(STITreeCluster cl, int index){
            _lineages.set(index);
            if(_coverage == null){
                _coverage = new STITreeCluster(cl);
            }
            else{
                _coverage = _coverage.merge(cl);
            }
        }


        public void mergeCluster(int coalTo, int[] coalFrom){
            if(_lineages.get(coalFrom[0]) && _lineages.get(coalFrom[1])){
                _lineages.set(coalFrom[0], false);
                _lineages.set(coalFrom[1], false);
                _lineages.set(coalTo, true);
            }
        }


        public void setExtraLineage(int xl){
            _xl = xl;
        }

        public void addExtraLineage(int xl){
            _xl += xl;
        }

        public int getLineageCount(){
            return _lineages.cardinality();
        }


        public String toString(List<STITreeCluster> gtClusters){
            String exp = "";
            for (int i = _lineages.nextSetBit(0); i >= 0; i = _lineages.nextSetBit(i+1)) {
                exp = exp + gtClusters.get(i);
            }
            exp = exp + "/[";
            /*
               for(int i=0; i<_netNodeLineages.length; i++){
                   exp = exp + _netNodeLineages[i][0].cardinality() + "/" + _netNodeLineages[i][1].cardinality();
                   if(i!=_netNodeLineages.length-1){
                       exp = exp + ",";
                   }
               }
               */
            for(int i=0; i<_netNodeIndex.length; i++){
                exp = exp + _netNodeIndex[i];
                if(i!=_netNodeLineages.length-1){
                    exp = exp + ",";
                }
            }
            exp = exp + "]:" + _xl;
            return exp;
        }

        public void addNetNodeChoice(int net, int index){
            _netNodeIndex[net] = index;
        }

        public void setNetNodeChoice(int[] choice){
            _netNodeIndex = choice.clone();
        }

        public void addNetNodeLineageNum(int net, int index, int lineage){
            _netNodeLineages[net][index].set(lineage);
        }

        /*
          public void addNetNodeLineageNum1(Configuration config){
              for(int i=0; i<_netNodeNum; i++){
                  _netNodeLinNum[i][0] += config._netNodeLinNum[i][0];
                  _netNodeLinNum[i][1] += config._netNodeLinNum[i][1];
              }
          }
          */
        public void setNetNodeLineageNum(BitSet[][] lineageNum){
            for(int i=0; i< _netNodeNum; i++){
                for(int j=0; j<2; j++){
                    _netNodeLineages[i][j] = (BitSet)(lineageNum[i][j].clone());
                }
            }
            /*
               int index = 0;
               for(int[] num: linNum){
                   _netNodeLinNum[index++] = num.clone();
               }
               */
        }

        public boolean equals(Object o) {
            if(!(o instanceof Configuration)){
                return false;
            }

            Configuration config = (Configuration) o;
            if (config._coverage.equals(_coverage) && config._lineages.equals(_lineages) && Arrays.equals(config._netNodeIndex,_netNodeIndex)) {
                return true;
            }
            else {
                return false;

            }
        }

        public int hashCode(){
            return _coverage.hashCode()+_lineages.hashCode();
        }

    }

}
