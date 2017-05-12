/*
 * Copyright (c) 2013 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.MC3Organizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.NormalMCMC;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.TwoStageMCMC;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary.Convergence;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.*;


/**
 * Created by IntelliJ IDEA.
 * User: dw20
 * Date: 3/29/16
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("MCMC_GT")
public class MCMC_GT extends CommandBaseFileOut{

    // for inference
    private List<List<NetworkNonEmpty>> _geneTrees;
    private List<NetworkNonEmpty> _startingNets = new ArrayList<NetworkNonEmpty>();
    private Map<String, List<String>> _taxonMap = null;

    private int _parallel = 1;
    private boolean _oneGTPerLocus = true;
    private double[] _weights = new double[] {
            0.01, 0.50, // kappa (dimension changing), kappa_1 (add reticulation)
            0.60, 0.80, // omega (non-topology changing), omega_1 (change-length) (note: if no reticulation edge, omega_1 wouldn't be used)
            0.33, 0.66, // move-tail, move-head
    };
    private long _chainLength = 1100000;
    private long _burnInLength = 100000;
    private long _sampleFrequency = 1000;

    private List<Double> _temperatures = new ArrayList<>();

    private double _poissonParam = 1.0;
    private long _seed = 1234567;
    private int _maxReti = 10;
    private boolean _usePseudo = false;
    private boolean _twoStage = false;
    private Class _stateClass;
    private Class _mcmcClass;

    // for summarization
    private List<String> _files = new ArrayList<>();

    public MCMC_GT(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                   Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                   Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader, Random rand){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams(){
        return 1;
    }

    @Override
    protected int getMaxNumParams(){
        return 40;
    }

    @Override
    protected boolean checkParamsForCommand(){
        boolean noError = true;

        // summary
        ParamExtractor sumParam = new ParamExtractor("sum", this.params, this.errorDetected);
        if(sumParam.ContainsSwitch) {

            if (sumParam.PostSwitchParam != null) {
                try {
                    if (!(sumParam.PostSwitchParam instanceof ParameterIdentList)) {
                        throw new RuntimeException();
                    }
                    ParameterIdentList temps = (ParameterIdentList) sumParam.PostSwitchParam;
                    for (String tExp : temps.Elements) {
                        System.out.println(tExp);
                        _files.add(tExp);
                    }
                } catch (NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -sum.",
                            sumParam.PostSwitchParam.getLine(), sumParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -sum.",
                        sumParam.SwitchParam.getLine(), sumParam.SwitchParam.getColumn());
            }

            noError = noError && checkForUnknownSwitches("sum");
            checkAndSetOutFile(sumParam);

        }  else {

            Parameter geneTreeParam = this.assertParameterIdentListOrSetList(0);
            noError = noError & geneTreeParam != null;

            // gene trees
            _geneTrees = new ArrayList<List<NetworkNonEmpty>>();
            if(geneTreeParam instanceof ParameterIdentList) {
                ParameterIdentList geneTreeList = (ParameterIdentList)geneTreeParam;
                for (String ident : geneTreeList.Elements) {
                    noError = noError && this.assertNetworkExists(ident, geneTreeParam.getLine(), geneTreeParam.getColumn());
                    if (noError) {
                        _geneTrees.add(Arrays.asList(this.sourceIdentToNetwork.get(ident)));
                    }
                }
            }
            else{
                _oneGTPerLocus = false;
                ParameterTaxonSetList geneTreeSetList = (ParameterTaxonSetList)geneTreeParam;
                for(Iterable<String> gtSet : geneTreeSetList.TaxonSetList)
                {
                    List<NetworkNonEmpty> geneTreesForOneLocus = new ArrayList<NetworkNonEmpty>();
                    for (String ident : gtSet) {
                        noError = noError && this.assertNetworkExists(ident, geneTreeParam.getLine(), geneTreeParam.getColumn());
                        if (noError) {
                            geneTreesForOneLocus.add(this.sourceIdentToNetwork.get(ident));
                        }
                    }
                    _geneTrees.add(geneTreesForOneLocus);
                }
            }

            // chain length
            ParamExtractor clParam = new ParamExtractor("cl", this.params, this.errorDetected);
            if(clParam.ContainsSwitch){
                if(clParam.PostSwitchParam != null)
                {
                    try
                    {
                        _chainLength = Long.parseLong(clParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized chain length " + clParam.PostSwitchValue,
                                clParam.PostSwitchParam.getLine(), clParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -cl.",
                            clParam.SwitchParam.getLine(), clParam.SwitchParam.getColumn());
                }
            }

            // burn-in length
            ParamExtractor blParam = new ParamExtractor("bl", this.params, this.errorDetected);
            if(blParam.ContainsSwitch){
                if(blParam.PostSwitchParam != null)
                {
                    try
                    {
                        _burnInLength = Long.parseLong(blParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized chain length " + blParam.PostSwitchValue,
                                blParam.PostSwitchParam.getLine(), blParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -bl.",
                            blParam.SwitchParam.getLine(), blParam.SwitchParam.getColumn());
                }
            }

            // sample frequency
            ParamExtractor sfParam = new ParamExtractor("sf", this.params, this.errorDetected);
            if(sfParam.ContainsSwitch){
                if(sfParam.PostSwitchParam != null)
                {
                    try
                    {
                        _sampleFrequency = Long.parseLong(sfParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized burn-in length " + sfParam.PostSwitchValue,
                                sfParam.PostSwitchParam.getLine(), sfParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -bl.",
                            sfParam.SwitchParam.getLine(), sfParam.SwitchParam.getColumn());
                }
            }

            // seed
            ParamExtractor sdParam = new ParamExtractor("sd", this.params, this.errorDetected);
            if(sdParam.ContainsSwitch){
                if(sdParam.PostSwitchParam != null)
                {
                    try
                    {
                        _seed = Long.parseLong(sdParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized random seed " + sdParam.PostSwitchValue,
                                sdParam.PostSwitchParam.getLine(), sdParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -sd.",
                            sdParam.SwitchParam.getLine(), sdParam.SwitchParam.getColumn());
                }
            }

            // parallel
            ParamExtractor plParam = new ParamExtractor("pl", this.params, this.errorDetected);
            if(plParam.ContainsSwitch){
                if(plParam.PostSwitchParam != null)
                {
                    try
                    {
                        _parallel = Integer.parseInt(plParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized number of thread " + plParam.PostSwitchValue,
                                plParam.PostSwitchParam.getLine(), plParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -pl.",
                            plParam.SwitchParam.getLine(), plParam.SwitchParam.getColumn());
                }
            }

            // maximum reticulation
            ParamExtractor mrParam = new ParamExtractor("mr", this.params, this.errorDetected);
            if(mrParam.ContainsSwitch)
            {
                if(mrParam.PostSwitchParam != null)
                {
                    try
                    {
                        _maxReti = Integer.parseInt(mrParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value of maximum number of reticulation nodes " + mrParam.PostSwitchValue,
                                mrParam.PostSwitchParam.getLine(), mrParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -mr.",
                            mrParam.SwitchParam.getLine(), mrParam.SwitchParam.getColumn());
                }
            }

            // Poisson parameter
            ParamExtractor ppParam = new ParamExtractor("pp", this.params, this.errorDetected);
            if(ppParam.ContainsSwitch)
            {
                if(ppParam.PostSwitchParam != null)
                {
                    try
                    {
                        _poissonParam = Double.parseDouble(ppParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized Poisson parameter " + ppParam.PostSwitchValue,
                                ppParam.PostSwitchParam.getLine(), ppParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -pp.",
                            ppParam.SwitchParam.getLine(), ppParam.SwitchParam.getColumn());
                }
            }

            // operation weights
            ParamExtractor wParam = new ParamExtractor("w", this.params, this.errorDetected);
            if(wParam.ContainsSwitch)
            {
                if(wParam.PostSwitchParam != null)
                {
                    try
                    {
                        if(!(wParam.PostSwitchParam instanceof ParameterIdentList)){
                            throw new RuntimeException();
                        }
                        ParameterIdentList weights = (ParameterIdentList)wParam.PostSwitchParam;
                        int index = 0;
                        for(String wExp: weights.Elements){
                            if(index < _weights.length){
                                double w = Double.parseDouble(wExp.trim());
                                if(w < 0){
                                    throw new RuntimeException();
                                }
                                _weights[index++] = w;
                            }
                            else{
                                throw new RuntimeException();
                            }
                        }
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Invalid value after switch -w.", wParam.PostSwitchParam.getLine(), wParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -w.", wParam.SwitchParam.getLine(), wParam.SwitchParam.getColumn());
                }
            }

            // temperatures
            ParamExtractor tpParam = new ParamExtractor("tp", this.params, this.errorDetected);
            if(tpParam.ContainsSwitch)
            {
                if(tpParam.PostSwitchParam != null)
                {
                    try
                    {
                        if(!(tpParam.PostSwitchParam instanceof ParameterIdentList)){
                            throw new RuntimeException();
                        }
                        ParameterIdentList temps = (ParameterIdentList) tpParam.PostSwitchParam;
                        for(String tExp: temps.Elements){
                            double t = Double.parseDouble(tExp.trim());
                            if(t < 0) throw new RuntimeException();
                            _temperatures.add(t);
                        }
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Invalid value after switch -tp.",
                                tpParam.PostSwitchParam.getLine(), tpParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -tp.",
                            tpParam.SwitchParam.getLine(), tpParam.SwitchParam.getColumn());
                }
            }

            // starting networks
            ParamExtractor snParam = new ParamExtractor("sn", this.params, this.errorDetected);
            if(snParam.ContainsSwitch){
                if(snParam.PostSwitchParam != null)
                {
                    try
                    {
                        if(!(snParam.PostSwitchParam instanceof ParameterIdentList)){
                            throw new RuntimeException();
                        }
                        ParameterIdentList networks = (ParameterIdentList) snParam.PostSwitchParam;
                        for(String ident: networks.Elements){
                            noError = noError && this.assertNetworkExists(ident,
                                    snParam.PostSwitchParam.getLine(), snParam.PostSwitchParam.getColumn());
                            if (noError) {
                                _startingNets.add(this.sourceIdentToNetwork.get(ident));
                            }
                        }
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Invalid value after switch -sn.",
                                snParam.PostSwitchParam.getLine(), snParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -sn.",
                            snParam.SwitchParam.getLine(), snParam.SwitchParam.getColumn());
                }
            }

            // taxon map
            ParamExtractor tmParam = new ParamExtractor("tm", this.params, this.errorDetected);
            if(tmParam.ContainsSwitch){
                ParamExtractorAllelListMap aaParam = new ParamExtractorAllelListMap("tm", this.params, this.errorDetected);
                noError = noError && aaParam.IsValidMap;
                if(aaParam.IsValidMap){
                    _taxonMap = aaParam.ValueMap;
                }
            }

            // pseudo likelihood
            ParamExtractor pseudoParam = new ParamExtractor("pseudo", this.params, this.errorDetected);
            if(pseudoParam.ContainsSwitch) {
                _usePseudo = true;
            }

            // two stage MCMC
            ParamExtractor twoStageParam = new ParamExtractor("-ts", this.params, this.errorDetected);
            if(twoStageParam.ContainsSwitch) {
                _twoStage = true;
            }

            noError = noError && checkForUnknownSwitches(
                    "cl", "bl", "sf", "sd", "pp", "mr", "pl", "tp", "sn", "tm", "pseudo", "ts");
            checkAndSetOutFile(clParam, blParam, sfParam, sdParam, ppParam,
                    mrParam, plParam, tpParam, snParam, tmParam, pseudoParam);
        }

        return  noError;
    }

    @Override
    protected String produceResult() {

        StringBuffer result = new StringBuffer("\n");
        System.out.println();

        if(_files.size() > 0) {
            Convergence conv = new Convergence(_files);
            conv.summarizeTopo(false);
            conv.computePSRF();
            conv.generateSRQ();
            conv.generateTracePlot();
            return result.toString();
        }

        List<List<MutableTuple<Tree,Double>>> gts = new ArrayList<List<MutableTuple<Tree,Double>>>();
        for(List<NetworkNonEmpty> geneTrees : _geneTrees){
            List<MutableTuple<Tree,Double>> gtsForOneLocus = new ArrayList<MutableTuple<Tree,Double>>();
            for(NetworkNonEmpty geneTree : geneTrees) {
                double weight = geneTree.TreeProbability.execute(new TreeProbabilityAlgo<Double, RuntimeException>() {
                    @Override
                    public Double forEmpty(TreeProbabilityEmpty empty) {
                        return 1.0;
                    }

                    @Override
                    public Double forNonEmpty(TreeProbabilityNonEmpty nonEmpty) {
                        return Double.parseDouble(nonEmpty.ProbString);
                    }
                });

                String phylonetGeneTree = NetworkTransformer.toENewickTree(geneTree);
                NewickReader nr = new NewickReader(new StringReader(phylonetGeneTree));
                STITree<Double> newtr = new STITree<Double>(true);
                try {
                    nr.readTree(newtr);
                } catch (Exception e) {
                    errorDetected.execute(e.getMessage(),
                            this._motivatingCommand.getLine(), this._motivatingCommand.getColumn());
                }
                Trees.removeBinaryNodes(newtr);
                gtsForOneLocus.add(new MutableTuple<Tree, Double>(newtr, weight));
            }
            gts.add(gtsForOneLocus);
        }
        List<Network> networks = new ArrayList<>();
        for(NetworkNonEmpty snet : _startingNets) {
            NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
            Network sNetwork = transformer.makeNetwork(snet);
            networks.add(sNetwork);
        }

        if(_usePseudo) {
            _stateClass = _oneGTPerLocus ? NetworkFromGTTPseudoSinglePerLocusState.class :
                    NetworkFromGTTPseudoMultiPerLocusState.class;
        } else {
            _stateClass = _oneGTPerLocus ? NetworkFromGTTSinglePerLocusState.class :
                    NetworkFromGTTMultiPerLocusState.class;
        }
        if(_twoStage) {
            _mcmcClass = TwoStageMCMC.class;
            _stateClass = _oneGTPerLocus ? NetworkFromGTTTwoStageSinglePerLocusState.class :
                    NetworkFromGTTTwoStageMultiPerLocusState.class;
        } else {
            _mcmcClass = NormalMCMC.class;
        }

        if(_temperatures.size() == 0) {
            _temperatures.add(1.0);
        }

        MC3Organizer mc3 = new MC3Organizer(
                _stateClass,
                _mcmcClass,
                networks,
                gts,
                _taxonMap,
                _chainLength,
                _burnInLength,
                _sampleFrequency,
                _seed,
                _temperatures,
                _poissonParam,
                _maxReti,
                _parallel,
                _weights
                );

        long startTime = System.currentTimeMillis();
        mc3.run();
        result.append(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));

        return result.toString();
    }

}
