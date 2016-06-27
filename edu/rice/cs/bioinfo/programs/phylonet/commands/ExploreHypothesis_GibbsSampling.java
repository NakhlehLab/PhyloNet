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
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferNetworkMP;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.*;


/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/27/16
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("ExploreHypothesis_GibbsSampling")
public class ExploreHypothesis_GibbsSampling extends CommandBaseFileOut{
    private List<List<NetworkNonEmpty>> _geneTrees = null;
    private NetworkNonEmpty _speciesNetwork = null;
    private HashMap<String, List<String>> _species2alleles = null;
    private double _maxBranchLength = 6;
    private int _numIterations = 11000;
    private int _burnIn = 1000;
    private int _sampleInterval = 100;
    private double _pruneThreshold = 0.01;
    private int _numProcessors = 1;
    private boolean _dentroscropeOutput = false;
    private boolean _pseudoLikelihood = false;
    private boolean _oneGTPerLocus = true;

    public ExploreHypothesis_GibbsSampling(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                           Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                           Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams(){
        return 2;
    }

    @Override
    protected int getMaxNumParams(){
        return 11;
    }

    @Override
    protected boolean checkParamsForCommand(){
        boolean noError = true;

        _speciesNetwork = this.assertAndGetNetwork(0);
        noError = noError && _speciesNetwork != null;

        Parameter geneTreeParam = this.assertParameterIdentListOrSetList(1);
        noError = noError && geneTreeParam != null;


        if(noError)
        {
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

            ParamExtractor aParam = new ParamExtractor("a", this.params, this.errorDetected);
            if(aParam.ContainsSwitch){
                ParamExtractorAllelListMap alParam = new ParamExtractorAllelListMap("a", this.params, this.errorDetected);
                noError = noError && alParam.IsValidMap;
                if(alParam.IsValidMap){
                    _species2alleles = alParam.ValueMap;
                }
            }

            ParamExtractor clParam = new ParamExtractor("cl", this.params, this.errorDetected);
            if(clParam.ContainsSwitch){
                if(clParam.PostSwitchParam != null)
                {
                    try
                    {
                        _numIterations = Integer.parseInt(clParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value of chain length " + clParam.PostSwitchValue, clParam.PostSwitchParam.getLine(), clParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -cl.", clParam.SwitchParam.getLine(), clParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor biParam = new ParamExtractor("bi", this.params, this.errorDetected);
            if(biParam.ContainsSwitch){
                if(biParam.PostSwitchParam != null)
                {
                    try
                    {
                        _burnIn = Integer.parseInt(biParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value of burn-in length " + biParam.PostSwitchValue, biParam.PostSwitchParam.getLine(), biParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -bi.", biParam.SwitchParam.getLine(), biParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor sfParam = new ParamExtractor("sf", this.params, this.errorDetected);
            if(sfParam.ContainsSwitch)
            {
                if(sfParam.PostSwitchParam != null)
                {
                    try
                    {
                        _sampleInterval = Integer.parseInt(sfParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value of sample frequency " + sfParam.PostSwitchValue, sfParam.PostSwitchParam.getLine(), sfParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -sf.", sfParam.SwitchParam.getLine(), sfParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor ptParam = new ParamExtractor("pt", this.params, this.errorDetected);
            if(ptParam.ContainsSwitch){
                if(ptParam.PostSwitchParam != null)
                {
                    try
                    {
                        _pruneThreshold = Double.parseDouble(ptParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value of threshold to prune networks " + ptParam.PostSwitchValue, ptParam.PostSwitchParam.getLine(), ptParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -pt.", ptParam.SwitchParam.getLine(), ptParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor mbParam = new ParamExtractor("mb", this.params, this.errorDetected);
            if(mbParam.ContainsSwitch){
                if(mbParam.PostSwitchParam != null)
                {
                    noError = noError && this.assertNetworkExists(mbParam.PostSwitchValue, mbParam.PostSwitchParam.getLine(), mbParam.PostSwitchParam.getColumn());
                    if(noError)
                    {
                        _maxBranchLength = Double.parseDouble(mbParam.PostSwitchValue);
                    }

                }
                else
                {
                    errorDetected.execute("Expected value after switch -mb.", mbParam.SwitchParam.getLine(), mbParam.SwitchParam.getColumn());
                }
            }


            ParamExtractor plParam = new ParamExtractor("pl", this.params, this.errorDetected);
            if(plParam.ContainsSwitch){
                if(plParam.PostSwitchParam != null)
                {
                    try
                    {
                        _numProcessors = Integer.parseInt(plParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value of number of processors " + plParam.PostSwitchValue, plParam.PostSwitchParam.getLine(), plParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -pl.", plParam.SwitchParam.getLine(), plParam.SwitchParam.getColumn());
                }
            }


            ParamExtractor psdParam = new ParamExtractor("psd", this.params, this.errorDetected);
            if(psdParam.ContainsSwitch){
                _pseudoLikelihood = true;

            }


            ParamExtractor diParam = new ParamExtractor("di", this.params, this.errorDetected);
            if(diParam.ContainsSwitch)
            {
                _dentroscropeOutput = true;
            }

            noError = noError && checkForUnknownSwitches("a", "cl", "bi", "sf", "pt", "mb", "psd", "pl", "di");
            checkAndSetOutFile(aParam, clParam, biParam, sfParam, ptParam, mbParam, psdParam, plParam, diParam);
        }



        return  noError;
    }

    @Override
    protected String produceResult() {
        StringWriter result = new StringWriter();
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

        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        Network speciesNetwork = transformer.makeNetwork(_speciesNetwork);

        GibbsSamplingForPruningNetworks gs;
        if(!_pseudoLikelihood){
            if(_oneGTPerLocus) {
                gs = new GibbsSamplingForPruningNetworksFromGTT_SingleTreePerLocus();
            }else{
                gs = new GibbsSamplingForPruningNetworksFromGTT_MultiTreesPerLocus();
            }
        }else{
            if(_oneGTPerLocus) {
                gs = new GibbsSamplingForPruningNetworksFromGTTPseudo_SingleTreePerLocus();
            }else{
                gs = new GibbsSamplingForPruningNetworksFromGTTPseudo_MultiTreesPerLocus();
            }
        }
        gs.setParameters(_maxBranchLength, _numIterations, _burnIn, _sampleInterval, _pruneThreshold, _numProcessors);
        gs.setOutputString(result);

        List<MutableTuple<Network,Integer>> samples = gs.sample(speciesNetwork, gts, _species2alleles);

        if(_dentroscropeOutput){
            result.append("\n");
            result.append("Visualize in Dendroscope: \n");
            for (MutableTuple<Network, Integer> tuple : samples) {
                edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.removeAllParameters(tuple.Item1);
                result.append(tuple.Item2 + ": " + tuple.Item1.toString() + "\n");
            }
        }

        return result.toString();

    }
}
