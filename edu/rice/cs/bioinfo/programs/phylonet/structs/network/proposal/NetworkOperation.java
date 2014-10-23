package edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal;

import java.util.Random;

/**
 * Created by wendingqiao on 10/20/14.
 */
public interface NetworkOperation<T> {

    public double operate();

    public void undo();
}
