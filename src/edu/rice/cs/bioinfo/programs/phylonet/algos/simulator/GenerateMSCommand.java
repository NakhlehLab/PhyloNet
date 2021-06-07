package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;

import java.io.ByteArrayInputStream;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 9/13/12
 * Time: 3:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class GenerateMSCommand {
    //int _numPopulations;
    private List<Tuple<NetNode<Integer>,Double>> _bundleList = new ArrayList<Tuple<NetNode<Integer>,Double>>();

    public static void main(String[] args) {
        try{
            GenerateMSCommand gen = new GenerateMSCommand();
            String st = "(a:2,(b:1,c:1):1);";
            RichNewickReaderAST reader = new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
            NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
            RichNewickReadResult<Networks> readResult = reader.read(new ByteArrayInputStream(st.getBytes()));
            Network net = transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());
            Map<String, List<String>> species2alleles = new Hashtable<String, List<String>>();
            List<String> temp1 = new ArrayList<String>();
            temp1.add("a1");
            species2alleles.put("a",temp1);
            List<String> temp2 = new ArrayList<String>();
            temp2.add("b1");
            temp2.add("b2");
            species2alleles.put("b",temp2);
            List<String> temp3 = new ArrayList<String>();
            temp3.add("c1");
            species2alleles.put("c",temp3);
            Map<String,String> MSName2AlleleName = new Hashtable<String, String>();
            int numGTs = 1000;
            double epsilon = 0.001;
            String output = gen.generateMS(net, numGTs, epsilon, species2alleles, MSName2AlleleName);
            System.out.println(output);
            System.out.println(MSName2AlleleName);
        }catch (Exception e){

        }
        //gen.generateMS(args[0],Integer.parseInt(args[1]),Double.parseDouble(args[2]));

    }

    public String generateMS(Network net, int numGTs, double epsilon, Map<String, List<String>> species2alleles, Map<String,String> MSName2AlleleName){

            //netString = "((_15#H2:1.023643045::0.66273154158198799290602210021461360156536102294921875,((_17#H3:0.3000685575::0.73019072997321987639196549935149960219860076904296875,((3:0.1000228525)_17#H3:0.1000228525::0.26980927002678012360803450064850039780139923095703125,2:0.200045705)_16:0.200045705)_3:0.63992774,1:1.04001915)_2:0.22965631)_14:0.22965631,((12:0.91582847,11:0.91582847)_11:0.26578128,((_13#H1:0.644354045::0.10751748781201087012959760613739490509033203125,((((((7:0.16402161)_13#H1:0.082010805::0.89248251218798912987040239386260509490966796875)_15#H2:0.082010805::0.33726845841801200709397789978538639843463897705078125,(9:0.02216053,8:0.02216053)_10:0.30588269)_9:0.14550185,6:0.47354507)_8:0.16241837,5:0.63596344)_7:0.04045963,4:0.67642307)_6:0.131952585)_12:0.131952585,10:0.94032824)_5:0.24128151)_4:0.31772202)_1;";
            //RichNewickReaderAST reader = new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
            //NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
            //RichNewickReadResult<Networks> readResult = reader.read(new ByteArrayInputStream(netString.getBytes()));
            //Network net = transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());
            Map<String, Integer> speciesName2MSName = new Hashtable<String, Integer>();
            List<Integer> numAllelesList = new ArrayList<Integer>();
            processNetwork(net, epsilon, species2alleles, speciesName2MSName, numAllelesList, MSName2AlleleName);
            int numPopulations = speciesName2MSName.size();
            int totalAlleles = 0;
            String outputCommand = "";
            for(int numAllele: numAllelesList){
                outputCommand += " "+numAllele;
                totalAlleles += numAllele;
            }
            outputCommand = "ms " + totalAlleles + " " + numGTs + " -T -I " + numPopulations + outputCommand;

            Map<BitSet, Integer> edge2population = new HashMap<BitSet, Integer>();

            for(Tuple<NetNode<Integer>,Double> bundle: _bundleList){
                NetNode<Integer> node = bundle.Item1;
                double height = bundle.Item2;
                if(node.isLeaf()){
                    NetNode<Integer> parent = node.getParents().iterator().next();
                    BitSet bs = new BitSet();
                    bs.set(parent.getData());
                    bs.set(node.getData());
                    edge2population.put(bs, speciesName2MSName.get(node.getName()));
                }
                else if(node.isTreeNode()){
                    Iterator<NetNode<Integer>> children = node.getChildren().iterator();
                    NetNode<Integer> child1 = children.next();
                    NetNode<Integer> child2 = children.next();
                    BitSet edge1 = new BitSet();
                    edge1.set(child1.getData());
                    edge1.set(node.getData());
                    int population1 = edge2population.get(edge1);
                    BitSet edge2 = new BitSet();
                    edge2.set(node.getData());
                    edge2.set(child2.getData());
                    int population2 = edge2population.get(edge2);
                    outputCommand += " -ej " + height + " " + population1 + " " + population2;

                    if(node.isRoot()){
                        break;
                    }
                    //System.out.println(node.getName()+":"+height1);
                    NetNode<Integer> parent = node.getParents().iterator().next();
                    BitSet bs = new BitSet();
                    bs.set(parent.getData());
                    bs.set(node.getData());
                    edge2population.put(bs, population2);
                }

                else if(node.isNetworkNode()){
                    NetNode<Integer> child = node.getChildren().iterator().next();

                    //System.out.println(node.getName()+":"+height);
                    BitSet childEdge = new BitSet();
                    childEdge.set(child.getData());
                    childEdge.set(node.getData());
                    int population = edge2population.get(childEdge);

                    Iterator<NetNode<Integer>> parents = node.getParents().iterator();
                    NetNode<Integer> parent1 = parents.next();
                    double hybridProb = node.getParentProbability(parent1);
                    outputCommand += " -es " + height + " " + population + " " + hybridProb;
                    BitSet parentEdge1 = new BitSet();
                    parentEdge1.set(parent1.getData());
                    parentEdge1.set(node.getData());
                    edge2population.put(parentEdge1, population);
                    NetNode<Integer> parent2 = parents.next();
                    BitSet parentEdge2 = new BitSet();
                    parentEdge2.set(parent2.getData());
                    parentEdge2.set(node.getData());
                    numPopulations++;
                    edge2population.put(parentEdge2, numPopulations);
                }
            }
    return outputCommand;

    }

    private void processNetwork(Network<Integer> net, double epsilon, Map<String, List<String>> species2alleles, Map<String, Integer> speciesName2MSName, List<Integer> numAllelesList, Map<String,String> MSName2AlleleName){
        int nodeIndex = 0;
        int leafIndex = 1;
        int populationIndex = 1;
        if(species2alleles == null){
            species2alleles = new Hashtable<String, List<String>>();
            for(NetNode<Integer> node: net.dfs()){
                node.setData(nodeIndex++);
                if(node.isLeaf()){
                    List<String> alleles = new ArrayList<String>();
                    alleles.add(node.getName());
                    species2alleles.put(node.getName(), alleles);
                    speciesName2MSName.put(node.getName(), leafIndex++);
                    numAllelesList.add(1);
                    MSName2AlleleName.put((populationIndex++)+"", node.getName());
                }
            }
        }
        else{
            for(NetNode<Integer> node: net.dfs()){
                node.setData(nodeIndex++);
                if(node.isLeaf()){
                    if(!species2alleles.containsKey(node.getName())){
                        throw new RuntimeException("The species " + node.getName() + " doesn't specify the number of alleles sampled.");
                    }
                    speciesName2MSName.put(node.getName(), leafIndex++);
                    List<String> alleles = species2alleles.get(node.getName());
                    int numAlleles = alleles.size();
                    numAllelesList.add(numAlleles);
                    for(String allele: alleles){
                        MSName2AlleleName.put((populationIndex++)+"", allele);
                    }
                    //System.out.println(node.getName());
                    //System.out.println(MSName2AlleleName);
                }
            }
        }
        if(species2alleles.size()!=leafIndex-1){
            throw new RuntimeException("The mapping for alleles doesn't match the network.");
        }
        Map<Integer, Double> node2height = new HashMap<Integer, Double>();
        for(NetNode<Integer> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(net)){
            double height;
            if(node.isLeaf()){
                node2height.put(node.getData(), 0.0);
                height = 0;
            }
            else if(node.isTreeNode()){
                Iterator<NetNode<Integer>> children = node.getChildren().iterator();
                NetNode<Integer> child1 = children.next();
                NetNode<Integer> child2 = children.next();
                Double height1 = node2height.get(child1.getData()) + child1.getParentDistance(node)*0.5;
                Double height2 = node2height.get(child2.getData()) + child2.getParentDistance(node)*0.5;
                if(Math.abs(height1-height2)>epsilon){
                    throw new RuntimeException("Not ultrametric");
                }
                height = Math.max(height1,height2);
                node2height.put(node.getData(), height);
            }
            else{
                NetNode<Integer> child = node.getChildren().iterator().next();
                height = node2height.get(child.getData()) + child.getParentDistance(node)*1;
                node2height.put(node.getData(),height);
            }
            Tuple<NetNode<Integer>,Double> bundle = new Tuple<NetNode<Integer>,Double>(node, height);
            int i = 0;
            for(; i< _bundleList.size(); i++){
                if(_bundleList.get(i).Item2>height){
                    break;
                }
            }
            _bundleList.add(i, bundle);
        }
    }


/*
    private class Bundle{
        NetNode<Integer> _node;
        double _time;

        public Bundle(NetNode<Integer> node, double t){
            _time = t;
            _node = node;
        }
    }
*/
}
