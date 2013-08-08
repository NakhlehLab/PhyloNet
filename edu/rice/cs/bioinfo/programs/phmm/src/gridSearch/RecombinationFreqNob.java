package gridSearch;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

import runHmm.AllInformation;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

/**
 * This knob is used for the recombination frequency parameter from
 * the TransitionProbabilityParameters class.
 * @author k3kathy
 *
 */
public class RecombinationFreqNob extends Nob {

    private Hmm<ObservationMap> thisHmm;

    private double[][] backupMatrix;

    private TransitionProbabilityParameters probsParam;
    private ArrayList<HiddenState> trees_states;
    private Map<Network<Double>,Set<HiddenState>> parentalTreeClasses;

    public RecombinationFreqNob(int gIn, double minIn, double maxIn,
            Hmm<ObservationMap> hmmIn, TransitionProbabilityParameters probsParamIn, ArrayList<HiddenState> treeStatesIn,
            Map<Network<Double>, Set<HiddenState>> parentalTreeIn) {
        super(gIn, minIn, maxIn);
        this.thisHmm = hmmIn;
        this.probsParam = probsParamIn;
        this.trees_states = treeStatesIn;
        this.parentalTreeClasses = parentalTreeIn;
    }

    public void set_param(double value) {
        backupParam = probsParam.getRecombinationFrequency();
        backupMatrix = thisHmm.getTransitionMatrix();
        probsParam.setRecombinationFrequency(value);
        double[][] newTransition = AllInformation.calculateAij(trees_states, value, probsParam.getHybridizationFrequency(),
                parentalTreeClasses);
        thisHmm.setTransitionMatrix(newTransition);
    }


    public double get_param() {
        return probsParam.getRecombinationFrequency();
    }

    public void restoreParameterValue() {
        probsParam.setRecombinationFrequency(backupParam);
        thisHmm.setTransitionMatrix(backupMatrix);
    }

}
