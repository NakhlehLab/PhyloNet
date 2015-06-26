package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.prototypes;


import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ApplyNetwork
{
public static void main(String[] args)
{
//    double[] lengths = {      1.1229271169652477,
//            0.8086422644900289,
//            0.26934975645322495,
//            3.9516145290261964,
//            0.20335872513372044};

    Network net = HmmNetworkUtils.fromENewickString("((BACK:2.0,ANC#H1:1.0):3.0,(ANC2#H2:3.0,(((GERMAN:1)ANC#H1:0)ANC2#H2:1.0,S:2):2):1);");
    Map<String,List<String>> alleleMap = new HashMap<String,List<String>>();

    alleleMap.put("BACK", Arrays.asList("Back"));
    alleleMap.put("S", Arrays.asList("S"));
    alleleMap.put("GERMAN", Arrays.asList("V"));
    System.out.println(HmmNetworkUtils.generateParentalTrees(net, alleleMap));
//
//    new NetworkLengthApplier(net,lengths).apply();

    System.out.println(net.toString());

}



}
