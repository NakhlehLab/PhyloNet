package edu.rice.cs.bioinfo.programs.phylonet.algos.butterfly;
/*
 * @ClassName:   DivideTrinets
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        3/6/22 8:40 PM
 */

import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.NetworkUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

public class DivideTrinets {
    /* Constructor */
    public DivideTrinets() {


    }

    public static void getTrinets(String snet){
        Network net = Networks.readNetwork(snet);
        List<Network> subNetworks = NetworkUtils.genAllSubNetworks(net, 3);
        for (Network n: subNetworks){
            System.out.println(n.toString());
        }

    }

    public static void main(String[] args) {
        String snet = "(((((((((HmelRef:6.0026877126907414E-6,HmelDisco:6.0026877126907414E-6:0.03319278608168805))#H1,Hpar:8.186901461975774E-4):0.01053669398658851,Hnum:0.011355384132786087):0.003722008703914167,Hbes:0.015077392836700254:0.04535729589631552):0.00550754040785891,((#H1,Hcyd:0.0013730107404778516):6.570567307460737E-4,Htim:0.0020300674712239253):0.01855486577333524):0.01836866117775692,Ldor:0.038953594422316086:0.12359338863788343):0.02922494602992582,(((Hhsa:0.0140330968390893,Htel:0.0140330968390893:0.12473817442178044):0.010713036854189815,((HeraRef:5.24482622065399E-5,HeraDisco:5.24482622065399E-5:0.025365861013938236):0.008250669663750003,Hhim:0.008303117925956543):0.016443015767322572):0.0073007414560471405,(Hsar:0.021688834546025864,Hdem:0.021688834546025864:0.05148570594947243):0.010358040603300391):0.03613166530291565):0.0929227321519819:0.027556360454750732,Avan:0.1611012726042238:0.030306833628154083);";
        getTrinets(snet);
    }
}
