package edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal;

/**
 * Created by wendingqiao on 10/20/14.
 */
public interface NetworkOperation {

    public double operate();

    public void undo();
}
