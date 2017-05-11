package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.List;
import java.util.Map;

/**
 * Created by dw20 on 5/11/17.
 */
public interface GTTLikelihood {

    double computeProbability(Network<Object> speciesNetwork, List distinctTrees,
                              Map<String,List<String>> species2alleles, List gtCorrespondences);


    void summarizeData(List originalGTs, Map<String,String> allele2species,
                       List dataForStartingNetwork, List dataForInferNetwork, List treeCorrespondences);

    void setParallel(int numThreads);

}
