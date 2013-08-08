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
 * This knob is used for the hybridization frequency parameter from
 * the TransitionProbabilityParameters class.
 * @author k3kathy
 *
 */
public class HybridizationFreqNob extends Nob {

    private Hmm<ObservationMap> thisHmm;

    private double[][] backupMatrix;

    private TransitionProbabilityParameters probsParam;
    private ArrayList<HiddenState> trees_states;
    private Map<Network<Double>,Set<HiddenState>> parentalTreeClasses;

    public HybridizationFreqNob(int gIn, double minIn, double maxIn,
            Hmm<ObservationMap> hmmIn, TransitionProbabilityParameters probsParamIn,
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

    public void restoreParameterValue() {
        probsParam.setHybridizationFrequency(backupParam);
        thisHmm.setTransitionMatrix(backupMatrix);
    }

}
