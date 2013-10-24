package gridSearch;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;


import runHmm.runHmm;
//import runHmm.AllInformation;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF.CoalescePattern;

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
    private Map<Network<CoalescePattern[]>,Set<HiddenState>> parentalTreeClasses;

    protected runHmm runHmmObject;

    public HybridizationFreqNob(int gIn, double minIn, double maxIn,
            Hmm<ObservationMap> hmmIn, TransitionProbabilityParameters probsParamIn,
            ArrayList<HiddenState> treeStatesIn,
				Map<Network<CoalescePattern[]>, Set<HiddenState>> parentalTreeIn, runHmm inRunHmmObject) {
        super(gIn, minIn, maxIn);
        this.thisHmm = hmmIn;
        this.probsParam = probsParamIn;
        this.trees_states = treeStatesIn;
        this.parentalTreeClasses = parentalTreeIn;
	runHmmObject = inRunHmmObject;
    }

    public void set_param(double value) {
        backupParam = probsParam.getHybridizationFrequency();
        backupMatrix = thisHmm.getTransitionMatrix();
        probsParam.setHybridizationFrequency(value);
        // double[][] newTransition = AllInformation.calculateAij(trees_states,
        //         probsParam.getRecombinationFrequency(), value,
        //         parentalTreeClasses);
        // thisHmm.setTransitionMatrix(newTransition);
	runHmmObject.updateTransitionProbabilities();
    }


    public double get_param() {
        return probsParam.getHybridizationFrequency();
    }

    public void restoreParameterValue() {
        probsParam.setHybridizationFrequency(backupParam);
        thisHmm.setTransitionMatrix(backupMatrix);
    }

}
