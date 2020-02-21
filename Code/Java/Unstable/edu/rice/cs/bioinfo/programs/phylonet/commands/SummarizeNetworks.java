package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.summarize.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/*
 *@ClassName: summarizeNetworks
 *@Description
 *@Author: Zhen Cao
 *@Date:  2019-09-11 00:38
 *@Version: 1.0
 */
@CommandName("SummarizeNetworks")
public class SummarizeNetworks extends CommandBaseFileOut {
    private int _mode;
    private String _outputFile;
    private Map<Network, Double> _networkSupportMap;
    private List<NetworkNonEmpty> _networkNonEmptyList;
    private boolean _normalize = false;
    private String _OUT_DIRECTORY = System.getProperty("user.dir");


    public SummarizeNetworks(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                             Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);

    }

    private class Out {
        public Out(){
            try {
                PrintStream print = new PrintStream(_OUT_DIRECTORY+"/SummarizeNetworks.txt");
                System.setOut(print);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }

        public Out(String filename){
            try {
                PrintStream print = new PrintStream(filename);
                System.setOut(print);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }

    @Override
    protected int getMinNumParams(){
        return 2;
    }

    @Override
    protected int getMaxNumParams(){
        return 40;
    }

    @Override
    protected boolean checkParamsForCommand() {
        boolean noError = true;

        //mode
        ParameterIdent number = this.assertParameterIdent(1);
        try
        {
            _mode = Integer.parseInt(number.Content);
        }
        catch(NumberFormatException e)
        {
            errorDetected.execute("Mode must be specified. ", number.getLine(), number.getColumn());
            noError = false;

        }

        //networks
        _networkNonEmptyList = new ArrayList<>();
        Parameter networksParam = this.assertParameterIdentListOrSetList(0);
        noError = noError && networksParam != null;
        if(networksParam instanceof ParameterIdentList) {
            ParameterIdentList networkList = (ParameterIdentList)networksParam;
            for (String ident : networkList.Elements) {
                noError = noError && this.assertNetworkExists(ident, networksParam.getLine(), networksParam.getColumn());
                if (noError) {
                    _networkNonEmptyList.add(this.sourceIdentToNetwork.get(ident));
//                    _networkSupportMap.put(this.sourceIdentToNetwork.get(ident), );
//                    _geneTrees.add(Arrays.asList(this.sourceIdentToNetwork.get(ident)));
                }
            }
        }

        // output file
        ParamExtractor outfileParam = new ParamExtractor("o", this.params, this.errorDetected);
        if(outfileParam.ContainsSwitch){
            if(outfileParam.PostSwitchParam != null) {
                if (noError) {
                    _outputFile = outfileParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -o.",
                        outfileParam.SwitchParam.getLine(), outfileParam.SwitchParam.getColumn());
            }
        }

        //normalize
        ParamExtractor normalizeParam = new ParamExtractor("n", this.params, this.errorDetected);
        if(normalizeParam.ContainsSwitch){
            if(normalizeParam.PostSwitchParam != null && !normalizeParam.PostSwitchValue.startsWith("-"))
            {
                errorDetected.execute("No value expected after switch -n.", normalizeParam.SwitchParam.getLine(), normalizeParam.SwitchParam.getColumn());
            }
            else
            {
                _normalize = true;
            }
        }

        noError = noError && checkForUnknownSwitches(
                "o", "n"
        );
        checkAndSetOutFile(outfileParam,normalizeParam);

        return  noError;
    }



    /**
     * @Description:
     * @Param: method
     * @Param: networkMap
     * @return:
     * @Author: Zhen Cao
     * @Date: 2019-09-12
     */
    private Map feature(int method, Map<Network, Double> networkMap){
        switch (method){
            case 1:
                Map<Tree, Double> maxtreemap = maxTree.summarize(networkMap);
                System.out.println("---------------------max tree---------------------");
//                for(Tree t: maxtreemap.keySet()){
//                    System.out.println(t.toString()+"\tcount:"+maxtreemap.get(t).toString());
//                }
//                System.out.println("-------------------------End------------------------");
                return maxtreemap;
            case 2:
                Map<Network, Double> backboneMap = backboneNetwork.summarize(networkMap);
                System.out.println("---------------------backbone network---------------------");
//                for(Network backnet: backboneMap.keySet()){
//                    System.out.println(backnet.toString()+"\t\t"+backboneMap.get(backnet).toString());
//                }
//                System.out.println("-------------------------End------------------------");
                return backboneMap;
            case 3:
                Map<Network, Double> decomposedtreemap = decomposedTrees.summarize(networkMap);
                System.out.println("---------------------decomposed tree---------------------");
//                for(Network t: decomposedtreemap.keySet()){
//                    System.out.println(t.toString()+"\t\tcount:"+decomposedtreemap.get(t).toString());
//                }
//                System.out.println("-------------------------End------------------------");
                return decomposedtreemap;
            case 4:
                Map<netNodeTuple, Double> tripart = tripartition2.summarize(networkMap);
                System.out.println("---------------------tripartition---------------------");
//                for(netNodeTuple t: tripart.keySet()){
//                    System.out.println(t.toString()+"\t\tcount:"+tripart.get(t).toString());
//                }
//                System.out.println("-------------------------End------------------------");
                return tripart;
            case 5:
                Map<Tree, Double> majortreemap = majorTree.summarize(networkMap);
                System.out.println("---------------------major tree---------------------");
//                for(Tree t: majortreemap.keySet()){
//                    System.out.println(t.toString()+"\t\tcount:"+majortreemap.get(t).toString());
//                }
//                System.out.println("-------------------------End------------------------");
                return majortreemap;

            default:
                System.out.println("No such summarizing method!");
                return null;
        }
    }

    public static class VComparator implements Comparator<Map.Entry<Object, Double>>
    {
        public int compare(Map.Entry<Object, Double> mp1, Map.Entry<Object, Double> mp2)
        {
            if (mp2.getValue() - mp1.getValue() > 0.00001){
                return 1;
            }
            else if (mp2.getValue() - mp1.getValue() < -0.00001){
                return -1;

            }
            else{
                return 0;
            }
        }


    }

    @Override
    protected String produceResult() {
        long startTime = System.currentTimeMillis();
        System.out.println();
        Out o;
        if (_outputFile == null){
            o = new Out();
        }
        else{
            o = new Out(_outputFile);
        }
        System.out.println("SummarizeNetworks");
        _networkSupportMap = new HashMap<Network, Double>();
        for (NetworkNonEmpty network: _networkNonEmptyList){
            double weight = network.TreeProbability.execute(new TreeProbabilityAlgo<Double, RuntimeException>() {
                @Override
                public Double forEmpty(TreeProbabilityEmpty empty) {
                    return 1.0;
                }

                @Override
                public Double forNonEmpty(TreeProbabilityNonEmpty nonEmpty) {
                    return Double.parseDouble(nonEmpty.ProbString);
                }
            });

            NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
            Network netcandidate = transformer.makeNetwork(network);
            _networkSupportMap.put(netcandidate, weight);
        }

        if (_normalize){
            double total = 0;
            for (Map.Entry<Network, Double> entry: _networkSupportMap.entrySet()){
                total += entry.getValue();
            }
            for (Map.Entry<Network, Double> entry: _networkSupportMap.entrySet()){
                _networkSupportMap.put(entry.getKey(), entry.getValue()/total);
            }
        }

        Map<Object, Double> commonMap = feature(_mode, _networkSupportMap);
        VComparator vc=new VComparator();
        List<Map.Entry<Object, Double>> MapList = new ArrayList<>();
        MapList.addAll(commonMap.entrySet());
        Collections.sort(MapList, vc);

        for (Map.Entry<Object, Double> entry: MapList){
            System.out.println(entry.getKey().toString()+"\tsupport:"+entry.getValue());
        }
        System.out.println("---------------------END---------------------");

        System.out.printf("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0);

        return "Done";
    }
}
