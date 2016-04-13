package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.List;

/**
 * Created by wendingqiao on 11/2/14.
 */
public interface State {

    public double propose();

    public void undo();

    public double calculateLikelihood();

    public double calculatePrior(double k);

    public String getOperation();

    public List<Tuple<String,Double>> storeState(double posterior);

    public String toString();

    public int numOfReticulation();

    public void setNetwork(Network s);

    public Network getNetwork();

}
