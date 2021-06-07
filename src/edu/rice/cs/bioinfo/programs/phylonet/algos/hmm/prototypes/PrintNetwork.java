package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.prototypes;


import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util.NetworkLengthApplier;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class PrintNetwork
{

    public static void main(String[] args)
    {
        String networkString = "((BACK:2.0,ANC#H1:1.0):1.0,(GERMAN:1.0)ANC#H1:2.0);";
        Network net = HmmNetworkUtils.fromENewickString(networkString);
        Map<String, List<String>> species2alleles = new HashMap<>();
        for(Object leaf: net.getLeaves()){
            String name = ((NetNode)leaf).getName();
            species2alleles.put(name, Arrays.asList(name));
        }
        new NetworkLengthApplier(net,new double[]{
                0.4184180302587213,
                0.19053087474067623,
                2.473062480279139,
                4.372790327458064
        }, species2alleles).apply();

        System.out.println(net.toString());

        new NetworkLengthApplier(net,new double[]{
                0.2605287847619249,
                0.5839093337276636,
                0.6180823761171341,
                4.324511317655939
        }, species2alleles).apply();

        System.out.println(net.toString());
    }


}
