package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.optimizer;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.MC3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.State;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.summary.SampleSummary;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 5/5/18
 * Time: 9:48 AM
 * To change this template use File | Settings | File Templates.
 */
public class MLOptimizer {
    private State _state;
    public MLOptimizer(List<Alignment> alignments) {
        try {
            _state = new State(
                    (Utils._START_NET == null || Utils._START_NET.size() == 0) ? null :
                            Utils._START_NET.get(0),
                    (Utils._START_GT_LIST == null || Utils._START_GT_LIST.size() == 0) ? null :
                            Utils._START_GT_LIST.get(0),
                    alignments,
                    Utils._POISSON_PARAM,
                    Utils._TAXON_MAP
            );

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void run() {
        _state.calculateLikelihood();
        OptimizeAll optimizer = new OptimizeAll(_state);
        double newScore = optimizer.propose();
        System.out.println(newScore);
        System.out.println(Networks.getFullString(_state.getNetworkObject()));
    }
}
