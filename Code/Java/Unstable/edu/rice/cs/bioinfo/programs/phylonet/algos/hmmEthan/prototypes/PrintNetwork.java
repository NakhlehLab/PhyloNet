package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.prototypes;


import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util.NetworkLengthApplier;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

public class PrintNetwork
{

    public static void main(String[] args)
    {
        String networkString = "((BACK:2.0,ANC#H1:1.0):1.0,(GERMAN:1.0)ANC#H1:2.0);";
        Network net = HmmNetworkUtils.fromENewickString(networkString);
        new NetworkLengthApplier(net,new double[]{
                0.4184180302587213,
                0.19053087474067623,
                2.473062480279139,
                4.372790327458064
        }).apply();

        System.out.println(net.toString());

        new NetworkLengthApplier(net,new double[]{
                0.2605287847619249,
                0.5839093337276636,
                0.6180823761171341,
                4.324511317655939
        }).apply();

        System.out.println(net.toString());
    }


}
