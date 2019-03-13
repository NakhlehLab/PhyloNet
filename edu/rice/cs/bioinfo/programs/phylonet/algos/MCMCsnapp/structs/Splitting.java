package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.LineageConfiguration;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.RPattern;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimLCInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;

import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 12/2/18
 * Time: 8:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class Splitting extends StateNode {
    public static final boolean CANCELED_LIKELIHOOD = true;
    public LineageConfiguration lineageConfigurations;
    public double logWeight = 0.0;

    public LineageConfiguration lineageConfigurationsPrev = null;

    private UltrametricNetwork _network = null;
    private RPattern _RPattern;
    private SimLCInNetwork _simlc = null;
    private BiAllelicGTR _model = null;
    private Map<String, List<String>> _species2alleles = null;

    public Splitting(UltrametricNetwork network, Map<String, List<String>> species2alleles, RPattern RPattern, BiAllelicGTR model) {
        _network = network;
        _model = model;
        _species2alleles = species2alleles;
        _RPattern = RPattern;
        _simlc = new SimLCInNetwork(model, Randomizer.getRandom().nextLong());
    }

    public Splitting() {
    }

    public RPattern getRPattern() {
        return _RPattern;
    }

    @Override
    public double propose() {
        lineageConfigurationsPrev = lineageConfigurations;
        lineageConfigurations = _simlc.generateLC(_network.getNetwork(), _species2alleles);
        setDirty(true);
        return 0.0;
    }

    @Override
    public boolean isValid() {
        return true;
    }

    @Override
    public boolean mayViolate() {
        return _operator.mayViolate();
    }

    @Override
    public double logDensity() {


        return Utils.INVALID_MOVE;
    }

    @Override
    public void accept() {
        _dirty = false;
        lineageConfigurationsPrev = null;
    }

    @Override
    public void reject() {
        _dirty = false;
    }

    @Override
    public void undo() {
        lineageConfigurations = lineageConfigurationsPrev;
        lineageConfigurationsPrev = null;
    }
}
