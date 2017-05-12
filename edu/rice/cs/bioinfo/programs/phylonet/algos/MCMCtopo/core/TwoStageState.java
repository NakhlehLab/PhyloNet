package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.List;

/**
 * Created by wendingqiao on 11/2/14.
 */
public interface TwoStageState extends State {

    public double calculatePseudoLikelihood();

}
