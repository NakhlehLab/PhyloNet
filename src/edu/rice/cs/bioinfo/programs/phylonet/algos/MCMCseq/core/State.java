package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.distribution.SpeciesNetPriorDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 2/2/16.
 */
public class State {

    // inferred or fixed topologies/parameters
    List<UltrametricTree> _geneTrees;       // (gene tree and heights)s
    private UltrametricNetwork _speciesNet;         // species network, heights, population sizes
    private double _gtOpWeight;
    private double _popSizeParamWeight = 0.005;
    private PopulationSize _populationSize;
    private int count = 0;
    private int gtMoveNum = -1;

    private StateNode _moveNode;

    // prior distribution
    private SpeciesNetPriorDistribution _priorDistribution;

    public State (String network,
                  Map<String, String> trees,
                  List<Alignment> alignments,
                  double poissonParameter,
                  Map<String, List<String>> species2alleles
    ) {

        this._geneTrees = new ArrayList<>(alignments.size());
        for(int i = 0; i < alignments.size(); i++) {
            if(trees == null) {
                this._geneTrees.add(new UltrametricTree(alignments.get(i)));
            } else {
                this._geneTrees.add(new UltrametricTree(trees.get(alignments.get(i).getName()), alignments.get(i)));
            }
        }
        Utils.NUM_LOCI = alignments.size();
        this._speciesNet = new UltrametricNetwork(network, this._geneTrees, species2alleles);
        this._populationSize = new PopulationSize();
        this._priorDistribution = new SpeciesNetPriorDistribution(poissonParameter, _populationSize);
        _gtOpWeight = 1.0 - Math.min(0.3, Math.max(0.1, 8.0 / (this._geneTrees.size() + 8.0)));
        if(Utils.ONLY_BACKBONE_OP) _gtOpWeight = 0.5;
    }

    /**
     * Serialize state
     * @return    serialization string
     */
    public String toString() {
        return "[" + _speciesNet.getNetwork().getRoot().getRootPopSize() + "]"
                + _speciesNet.getNetwork().toString();
    }

    /**
     * Serialize state to list
     * @return    serialization string
     */
    public List<String> toList() {
        List<String> list = new ArrayList<>();
        list.add(_speciesNet.getNetwork().toString());
        for(UltrametricTree t : this._geneTrees) {
            list.add(t.toString());
        }
        if(Utils._ESTIMATE_POP_SIZE) {
            list.add(_populationSize.toString());
        }
        if (Utils.SAMPLE_MUTATION_RATE) {
            for (UltrametricTree t : this._geneTrees) {
                list.add(Double.toString(t.get_mutationRate()));
            }
        }
        return list;
    }

    /**
     * @return    the hastings ratio
     */
    public double propose() {
        double rand = Randomizer.getRandomDouble();

        if(Utils._ESTIMATE_POP_SIZE && Utils._ESTIMATE_POP_PARAM && rand < _popSizeParamWeight) {
            _moveNode = _populationSize;
        }
        // tree operation
        else if((rand < _gtOpWeight || Utils._FIX_NET) && !Utils._FIX_GENE_TREES) {
            _moveNode = _geneTrees.get(gtMoveNum = Randomizer.getRandomInt(_geneTrees.size()));
            if(Utils.DEBUG_MODE) {
                System.out.println("Proposing gene tree operation.");
            }
        }
        // network operation
        else {
            gtMoveNum = -1;
            _moveNode = _speciesNet;
            if(Utils.DEBUG_MODE) {
                System.out.println("Proposing network operation.");
            }
        }

        if(count == 10000) {
            System.out.print("");
        }

        double logHR = _moveNode.propose();
        _moveNode.setDirty(true);
        if(getOperation().getName().contains("Scale-All")) {
            if(Utils._FIX_GENE_TREES) {
                return Utils.INVALID_MOVE;
            }

            for(UltrametricTree ut : _geneTrees) {
                ut.setDirty(true);
            }
        }

        count++;
        if(Utils.DEBUG_MODE) {
            System.out.println("Count: " + count + " Purposed: " + getOperation().getName() + " mayViolate: " + _moveNode.mayViolate() + " : " + _speciesNet.toString());
            if (getOperation().getName().contains("DeltaExchange")){
                for(UltrametricTree ut : _geneTrees) {
                    System.out.print(ut.get_mutationRate()+"\t");
                }
                System.out.println();
            }
        }
//        if (getOperation().getName().contains("DeltaExchange")){
//            for(UltrametricTree ut : _geneTrees) {
//                System.out.print(ut.get_mutationRate()+"\t");
//            }
//            System.out.println();
//        }

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
        for (UltrametricTree gt : _geneTrees) gt.reject();

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
        for (UltrametricTree gt : _geneTrees) gt.accept();

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
        return _priorDistribution.logPrior(_speciesNet.getNetwork(), _geneTrees);
    }

    /**
     * @return    likelihood value
     */
    public double calculateLikelihood() {
        double logL = _speciesNet.logDensity();
//        System.out.println("species log density:"+logL);
        for(UltrametricTree gt : _geneTrees) {
            logL += gt.logDensity();
//            System.out.println("gt log density:"+gt.logDensity());
        }
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
            for(UltrametricTree gt : _geneTrees) {
                if(gt.isValid()) continue;
                System.err.println("Invalid gene tree - not ultrametric");
                return false;
            }
            return true;
        }
    }

    public String getNetwork() {
        return _speciesNet.getNetwork().toString();
    }

    public Network getNetworkObject() {
        return _speciesNet.getNetwork();
    }

    public UltrametricNetwork getUltrametricNetworkObject() {
        return _speciesNet;
    }

    protected void setPreBurnInParams() {
        if(Utils._ESTIMATE_POP_SIZE) {
            double rootPS = _speciesNet.getNetwork().getRoot().getRootPopSize();
            Utils._POP_SIZE_MEAN = rootPS;
            Utils._POP_SIZE_WINDOW_SIZE = 0.05 * rootPS;
            System.out.println("SET _POP_SIZE_MEAN = " + Utils._POP_SIZE_MEAN);
            System.out.println("SET _POP_SIZE_WINDOW_SIZE = " + Utils._POP_SIZE_WINDOW_SIZE);
        }
        List<Double> list = new ArrayList<>();
        for(NetNode node : Networks.postTraversal( _speciesNet.getNetwork())) {
            for(Object c : node.getChildren()) {
                NetNode child = (NetNode) c;
                list.add(child.getParentDistance(node));
            }
        }
        Collections.sort(list);
        Utils._TIME_WINDOW_SIZE = list.get(list.size() / 2) * 0.05;
        System.out.println("SET _TIME_WINDOW_SIZE = " + Utils._TIME_WINDOW_SIZE);
    }

}
