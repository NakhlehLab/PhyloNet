package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.topo;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 1/19/18
 * Time: 3:47 PM
 * To change this template use File | Settings | File Templates.
 */

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.dimension.AddReticulation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.dimension.DeleteReticulation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

import java.util.HashSet;
import java.util.Set;

/**
 * This is NOT a MCMC move! Only work with MLE!
 */

public class ReplaceReticulation extends NetworkOperator {

    private AddReticulation addReticulation;
    private DeleteReticulation deleteReticulation;

    public ReplaceReticulation(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        if(_network.getNetwork().getReticulationCount() >= 0) {
            return Utils.INVALID_MOVE;
        }

        deleteReticulation = new DeleteReticulation(_network);
        deleteReticulation.propose();

        addReticulation = new AddReticulation(_network);
        addReticulation.propose();

        return 0.0;
    }

    @Override
    public void undo() {
        addReticulation.undo();
        deleteReticulation.undo();
        addReticulation = null;
        deleteReticulation = null;
    }

    @Override
    public String getName() {
        return "Replace-Reticulation-NotMCMC";
    }

}
