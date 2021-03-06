package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.List;
import java.util.Map;

/**
 * Created by dw20 on 5/11/17.
 */
public interface GTTLikelihood {

    double computeProbability(Network<Object> speciesNetwork, List list1,
                              Map<String,List<String>> species2alleles, List list2);


    void summarizeData(List originalGTs, Map<String,String> allele2species,
                       List dataForStartingNetwork, List list1, List list2);

    void setParallel(int numThreads);

}
