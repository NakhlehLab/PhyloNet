package gridSearch;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

import phylogeny.EvoTree;
import runHmm.AllInformation;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;
import edu.rice.cs.bioinfo.library.programming.BijectiveHashtable;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

/**
 * This knob is used for the hybridization frequency parameter from
 * the TransitionProbabilityParameters class.
 * @author k3kathy
 *
 */
public class HybridizationFreqNob extends Nob {

    private Hmm thisHmm;

    private double[][] backupMatrix;

    private TransitionProbabilityParameters probsParam;
    private ArrayList<HiddenState> trees_states;
    private Map<Network<Double>,Set<HiddenState>> parentalTreeClasses;

    public HybridizationFreqNob(int gIn, double minIn, double maxIn,
            Hmm hmmIn, TransitionProbabilityParameters probsParamIn,
            ArrayList<HiddenState> treeStatesIn,
            Map<Network<Double>, Set<HiddenState>> parentalTreeIn) {
        super(gIn, minIn, maxIn);
        this.thisHmm = hmmIn;
        this.probsParam = probsParamIn;
        this.trees_states = treeStatesIn;
        this.parentalTreeClasses = parentalTreeIn;
    }

    public void set_param(double value) {
        backupParam = probsParam.getHybridizationFrequency();
        backupMatrix = thisHmm.getTransitionMatrix();
        probsParam.setHybridizationFrequency(value);
        double[][] newTransition = AllInformation.calculateAij(trees_states,
                probsParam.getRecombinationFrequency(), value,
                parentalTreeClasses);
        thisHmm.setTransitionMatrix(newTransition);
    }


    public double get_param() {
        return probsParam.getHybridizationFrequency();
    }

    @Override
    public void restoreParameterValue() {
        // TODO Auto-generated method stub
        probsParam.setHybridizationFrequency(backupParam);
        thisHmm.setTransitionMatrix(backupMatrix);
    }

}
