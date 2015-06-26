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
 * Date: 3/13/12
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("infernetwork_mp")
public class InferNetwork_MP extends CommandBaseFileOut{
    private HashMap<String, List<String>> _species2alleles = null;
    private HashMap<String, String> _allele2species = null;
    private List<List<NetworkNonEmpty>> _geneTrees;
    private double _bootstrap = 100;
    private NetworkNonEmpty _startSpeciesNetwork = null;
    private int _maxReticulations;
    private long _maxExaminations = -1;
    private int _moveDiameter = -1;
    private int _reticulationDiameter = -1;
    private int _maxFailure = 100;
    private int _returnNetworks = 1;
    private int _numProcessors = 1;
    private boolean _dentroscropeOutput = false;
    private int _numRuns = 5;
    private double[] _operationWeight = {0.1,0.1,0.15,0.55,0.15,0.15};
    private Long _seed = null;

    private Set<String> _fixedHybrid = new HashSet<String>();


    public InferNetwork_MP(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
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
        return 25;
    }

    @Override
    protected boolean checkParamsForCommand(){
        boolean noError = true;

        Parameter geneTreeParam = this.assertParameterIdentListOrSetList(0);
        noError = noError && geneTreeParam != null;
        _geneTrees = new ArrayList<List<NetworkNonEmpty>>();

        ParameterIdent number = this.assertParameterIdent(1);
        try
        {
            _maxReticulations = Integer.parseInt(number.Content);
        }
        catch(NumberFormatException e)
        {
            errorDetected.execute("Maximum number of reticulation nodes must be specified. ", number.getLine(), number.getColumn());
            noError = false;

        }


        if(noError)
        {

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

                ParamExtractorAllelMap aaParam = new ParamExtractorAllelMap("a", this.params, this.errorDetected);
                if(aaParam.IsValidMap){
                    _allele2species = aaParam.ValueMap;
                }
            }

            ParamExtractor mdParam = new ParamExtractor("md", this.params, this.errorDetected);
            if(mdParam.ContainsSwitch){
                if(mdParam.PostSwitchParam != null)
                {
                    try
                    {
                        _moveDiameter = Integer.parseInt(mdParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized maximum diameter for network search " + mdParam.PostSwitchValue, mdParam.PostSwitchParam.getLine(), mdParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -md.", mdParam.SwitchParam.getLine(), mdParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor rdParam = new ParamExtractor("rd", this.params, this.errorDetected);
            if(rdParam.ContainsSwitch){
                if(rdParam.PostSwitchParam != null)
                {
                    try
                    {
                        _reticulationDiameter = Integer.parseInt(rdParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized maximum diameter for network search " + rdParam.PostSwitchValue, rdParam.PostSwitchParam.getLine(), rdParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -rd.", rdParam.SwitchParam.getLine(), rdParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor nParam = new ParamExtractor("n", this.params, this.errorDetected);
            if(nParam.ContainsSwitch)
            {
                if(nParam.PostSwitchParam != null)
                {
                    try
                    {
                        _returnNetworks = Integer.parseInt(nParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value of returned networks " + nParam.PostSwitchValue, nParam.PostSwitchParam.getLine(), nParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -n.", nParam.SwitchParam.getLine(), nParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor mParam = new ParamExtractor("m", this.params, this.errorDetected);
            if(mParam.ContainsSwitch){
                if(mParam.PostSwitchParam != null)
                {
                    try
                    {
                        _maxExaminations = new Long(Integer.parseInt(mParam.PostSwitchValue));
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized maximum number of networks for search " + mParam.PostSwitchValue, mParam.PostSwitchParam.getLine(), mParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -m.", mParam.SwitchParam.getLine(), mParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor sParam = new ParamExtractor("s", this.params, this.errorDetected);
            if(sParam.ContainsSwitch){
                if(sParam.PostSwitchParam != null)
                {
                    noError = noError && this.assertNetworkExists(sParam.PostSwitchValue, sParam.PostSwitchParam.getLine(), sParam.PostSwitchParam.getColumn());
                    if(noError)
                    {
                        _startSpeciesNetwork = this.sourceIdentToNetwork.get(sParam.PostSwitchValue);
                        //_geneTrees.add(this.sourceIdentToNetwork.get(ident));
                    }

                }
                else
                {
                    errorDetected.execute("Expected value after switch -s.", sParam.SwitchParam.getLine(), sParam.SwitchParam.getColumn());
                }
            }




            ParamExtractor bParam = new ParamExtractor("b", this.params, this.errorDetected);
            if(bParam.ContainsSwitch)
            {
                if(bParam.PostSwitchParam != null)
                {
                    try
                    {
                        _bootstrap = Double.parseDouble(bParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized bootstrap value " + bParam.PostSwitchValue, bParam.PostSwitchParam.getLine(), bParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -b.", bParam.SwitchParam.getLine(), bParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor hParam = new ParamExtractor("h", this.params, this.errorDetected);
            if(hParam.ContainsSwitch)
            {
                if(hParam.PostSwitchParam != null)
                {
                    try
                    {
                        if(!(hParam.PostSwitchParam instanceof ParameterIdentSet)){
                            throw new RuntimeException();
                        }
                        ParameterIdentSet hybrids = (ParameterIdentSet)hParam.PostSwitchParam;
                        for(String value: hybrids.Elements){
                            _fixedHybrid.add(value);
                        }

                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value after switch -h.", hParam.PostSwitchParam.getLine(), hParam.PostSwitchParam.getColumn());
                    }

                }
                else
                {
                    errorDetected.execute("Expected value after switch -h.", hParam.SwitchParam.getLine(), hParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor fParam = new ParamExtractor("f", this.params, this.errorDetected);
            if(fParam.ContainsSwitch)
            {
                if(fParam.PostSwitchParam != null)
                {
                    try
                    {
                        _maxFailure = Integer.parseInt(fParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value of maximum consecutive failure " + fParam.PostSwitchValue, fParam.PostSwitchParam.getLine(), fParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -n.", fParam.SwitchParam.getLine(), fParam.SwitchParam.getColumn());
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
                        errorDetected.execute("Unrecognized maximum number of networks for search " + plParam.PostSwitchValue, plParam.PostSwitchParam.getLine(), plParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -pl.", plParam.SwitchParam.getLine(), plParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor xParam = new ParamExtractor("x", this.params, this.errorDetected);
            if(xParam.ContainsSwitch)
            {
                if(xParam.PostSwitchParam != null)
                {
                    try
                    {
                        _numRuns = Integer.parseInt(xParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value after -x " + xParam.PostSwitchValue, xParam.PostSwitchParam.getLine(), xParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -x.", xParam.SwitchParam.getLine(), xParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor rsParam = new ParamExtractor("rs", this.params, this.errorDetected);
            if(rsParam.ContainsSwitch){
                if(rsParam.PostSwitchParam != null)
                {
                    try
                    {
                        _seed = Long.parseLong(rsParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized seed for network search " + rsParam.PostSwitchValue, rsParam.PostSwitchParam.getLine(), rsParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -rs.", rsParam.SwitchParam.getLine(), rsParam.SwitchParam.getColumn());
                }
            }

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
                        double total = 0;
                        for(String wExp: weights.Elements){
                            if(index < _operationWeight.length){
                                double w = Double.parseDouble(wExp.trim());
                                if(w<0){
                                    throw new RuntimeException();
                                }
                                _operationWeight[index++] = w;
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


            ParamExtractor diParam = new ParamExtractor("di", this.params, this.errorDetected);
            if(diParam.ContainsSwitch)
            {
                _dentroscropeOutput = true;
            }

            noError = noError && checkForUnknownSwitches("a","b","s","n", "m", "md", "rd", "di", "h","f","pl","x","rs","w");
            checkAndSetOutFile(aParam, bParam, sParam, nParam, mParam, mdParam, rdParam, diParam, hParam, fParam, plParam,xParam,rsParam,wParam);
        }



        return  noError;
    }

    @Override
    protected String produceResult() {
        StringWriter result = new StringWriter();

        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        for(List<NetworkNonEmpty> geneTreesForOneLocus : _geneTrees) {
            final int size = geneTreesForOneLocus.size();
            for (NetworkNonEmpty geneTree : geneTreesForOneLocus) {
                {
                    double weight = geneTree.TreeProbability.execute(new TreeProbabilityAlgo<Double, RuntimeException>() {
                        @Override
                        public Double forEmpty(TreeProbabilityEmpty empty) {
                            return 1.0/size;
                        }

                        @Override
                        public Double forNonEmpty(TreeProbabilityNonEmpty nonEmpty) {
                            return Double.parseDouble(nonEmpty.ProbString)/size;
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
                    if (_bootstrap < 100) {
                        if (Trees.handleBootStrapInTree(newtr, _bootstrap) == -1) {
                            throw new IllegalArgumentException("Input gene tree " + newtr + " have nodes that don't have bootstrap value");
                        }
                    }
                    for (TNode node : newtr.getNodes()) {
                        node.setParentDistance(TNode.NO_DISTANCE);
                    }

                    String exp = Trees.getLexicographicNewickString(newtr, _allele2species);
                    MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
                    if (existingTuple == null) {
                        existingTuple = new MutableTuple<Tree, Double>(newtr, weight);
                        exp2tree.put(exp, existingTuple);
                    } else {
                        existingTuple.Item2 += weight;
                    }
                }
            }
        }

        List<MutableTuple<Tree,Double>> tuples = new ArrayList<MutableTuple<Tree, Double>>();
        tuples.addAll(exp2tree.values());

        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        Network speciesNetwork = null;
        if(_startSpeciesNetwork!=null){
            speciesNetwork = transformer.makeNetwork(_startSpeciesNetwork);
        }


        //long start = System.currentTimeMillis();
        //InferILSNetworkParsimoniouslyParallelBackup inference = new InferILSNetworkParsimoniouslyParallelBackup();
        InferNetworkMP inference = new InferNetworkMP();
        inference.setSearchParameter(_maxExaminations, _maxFailure, _moveDiameter, _reticulationDiameter, speciesNetwork, _fixedHybrid, _numProcessors, _operationWeight, _numRuns, _seed);
        List<Tuple<Network, Double>> resultTuples = inference.inferNetwork(tuples,_species2alleles,_maxReticulations, _returnNetworks);
        //System.out.print(System.currentTimeMillis()-start);
        int index = 1;
        for(Tuple<Network, Double> tuple: resultTuples){
            result.append("\nInferred Network #" + index++ + ":");

            Network n = tuple.Item1;

            for(Object node : n.bfs())
            {
                NetNode netNode = (NetNode)node;
                if(!netNode.isLeaf())
                {
                    netNode.setName(NetNode.NO_NAME);
                }
            }

            result.append("\n" + n.toString());
            result.append("\n" + "Total number of extra lineages: " + tuple.Item2);

            if(_dentroscropeOutput){
                edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.removeAllParameters(n);
                result.append("\nVisualize in Dendroscope : " + n.toString());
            }
        }

        return result.toString();

    }
}
