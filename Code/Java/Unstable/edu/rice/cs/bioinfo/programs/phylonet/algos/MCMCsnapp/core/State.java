package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SpeciesNetPriorDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.PopulationSize;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 2/2/16.
 */
public class State {

    // inferred or fixed topologies/parameters
    private List<UltrametricTree> _geneTrees;       // (gene tree and heights)s
    private UltrametricNetwork _speciesNet;         // species network, heights, population sizes
    private double _gtOpWeight;
    private double _popSizeParamWeight = 0.005;
    private PopulationSize _populationSize;

    private StateNode _moveNode;

    // prior distribution
    private SpeciesNetPriorDistribution _priorDistribution;

    public State (String network,
                  List<String> trees,
                  List<Alignment> alignments,
                  double poissonParameter,
                  Map<String, List<String>> species2alleles,
                  BiAllelicGTR BAGTRModel
    ) {

        this._geneTrees = new ArrayList<>(alignments.size());
        for(int i = 0; i < alignments.size(); i++) {
            if(trees == null) {
                this._geneTrees.add(new UltrametricTree(alignments.get(i)));
            } else {
                this._geneTrees.add(new UltrametricTree(trees.get(i), alignments.get(i)));
            }
        }
        this._speciesNet = new UltrametricNetwork(network, _geneTrees, alignments, species2alleles, BAGTRModel);
        this._populationSize = new PopulationSize();
        this._priorDistribution = new SpeciesNetPriorDistribution(poissonParameter, _populationSize);
        //_gtOpWeight = 1.0 - Math.min(0.3, Math.max(0.1, 8.0 / (this._geneTrees.size() + 8.0)));
    }

    /**
     * Serialize state
     * @return    serialization string
     */
    public String toString() {
        return "[" + _speciesNet.getNetwork().getRoot().getRootPopSize() + "]"
                + _speciesNet.getNetwork().toString() + "\nTopo: " + Networks.getTopologyString(_speciesNet.getNetwork());
    }

    /**
     * Serialize state to list
     * @return    serialization string
     */
    public List<String> toList() {
        List<String> list = new ArrayList<>();
        list.add(_speciesNet.getNetwork().toString());
        /*for(UltrametricTree t : this._geneTrees) {
            list.add(t.toString());
        }*/
        if(Utils._ESTIMATE_POP_SIZE) {
            list.add(_populationSize.toString());
        }
        return list;
    }

    /**
     * @return    the hastings ratio
     */
    public double propose() {
        double rand = Randomizer.getRandomDouble();

        if(Utils._ESTIMATE_POP_SIZE && rand < _popSizeParamWeight) {
            _moveNode = _populationSize;
        }
        else {
            // always do network operation (jiafan)
            _moveNode = _speciesNet;
        }

        double logHR = _moveNode.propose();
        _moveNode.setDirty(true);
        if(getOperation().getName().contains("Scale-All")) {
            /*for(UltrametricTree ut : _geneTrees) {
                ut.setDirty(true);
            }*/
            //logHR = Utils.INVALID_MOVE;
        }

        if(getOperation().getName().contains("Add-Reticulation") &&
                _speciesNet.getNetwork().getReticulationCount() > Utils._NET_MAX_RETI) {
            logHR = Utils.INVALID_MOVE;
        } else if(_moveNode.mayViolate()){
            if(_speciesNet.isDirty() && !_priorDistribution.isValid(_speciesNet.getNetwork())) {
                if(Utils.DEBUG_MODE) System.err.println(getOperation());
                logHR = Utils.INVALID_MOVE;
            } else if(!_speciesNet.isValid()) {
                logHR = Utils.INVALID_MOVE;
            }
        }
        return logHR;
    }

    /**
     * Undo the proposal, restore the last state
     */
    public void undo(double logAlpha) {
        _moveNode.undo();
        _speciesNet.reject();
        //for (UltrametricTree gt : _geneTrees) gt.reject();

        _moveNode._operator._rejCounter++;
        _moveNode._operator._rejCorrectionCounter++;
        if (logAlpha != Utils.INVALID_MOVE)
            _moveNode._operator.optimize(logAlpha);
        _moveNode = null;
    }

    /**
     * accept the proposal
     */
    public void accept(double logAlpha) {
        _speciesNet.accept();
        //for (UltrametricTree gt : _geneTrees) gt.accept();

        _moveNode._operator._acCounter++;
        _moveNode._operator._acCorrectionCounter++;
        if (logAlpha != Utils.INVALID_MOVE)
            _moveNode._operator.optimize(logAlpha);
        _moveNode = null;
    }

    /**
     * @return   the last proposal name
     */
    public Operator getOperation() {
        return _moveNode.getOperation();
    }

    /**
     * @return    prior value
     */
    public double calculatePrior() {
        return _priorDistribution.logPrior(_speciesNet.getNetwork());
    }

    /**
     * @return    likelihood value
     */
    public double calculateLikelihood() {
        double logL = _speciesNet.logDensity();
        /*for(UltrametricTree gt : _geneTrees) {
            logL += gt.logDensity();
        }*/
        return logL;
    }

    public int numOfReticulation() {
        return _speciesNet.getNetwork().getReticulationCount();
    }

    // this function should only be called under DEBUG_MODE
    public boolean isValidState() {
        if(!Utils.DEBUG_MODE) return true;

        if(!_priorDistribution.isValid(_speciesNet.getNetwork())) {
            System.err.println("Invalid network");
            return false;
        } else if (!_speciesNet.isValid()) {
            System.err.println("Invalid temporal constraints");
            return false;
        } else if(!_speciesNet.isUltrametric()) {
            System.err.println("Invalid network - not ultrametric");
            return false;
        } else {
            /*for(UltrametricTree gt : _geneTrees) {
                if(gt.isValid()) continue;
                System.err.println("Invalid gene tree - not ultrametric");
                return false;
            }*/
            return true;
        }
    }

    public String getNetwork() {
        return _speciesNet.getNetwork().toString();
    }

}
