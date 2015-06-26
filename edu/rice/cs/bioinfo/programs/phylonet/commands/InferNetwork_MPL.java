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
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
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
 * User: yy9
 * Date: 3/13/12
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("infernetwork_MPL")
public class InferNetwork_MPL extends CommandBaseFileOut{
    private HashMap<String, List<String>> _taxonMap = null;
    private List<List<NetworkNonEmpty>> _geneTrees;
    private double _bootstrap = 100;
    private NetworkNonEmpty _startSpeciesNetwork = null;
    private int _maxReticulations;
    private long _maxExaminations = -1;
    private int _maxFailure = 100;
    private int _moveDiameter = -1;
    private int _reticulationDiameter = -1;
    private int _returnNetworks = 5;
    private int _maxRounds = 100;
    private int _maxTryPerBranch = 100;
    private double _maxBranchLength = 6;
    private double _improvementThreshold = 0.01;
    private double _Brent1 = 0.01;
    private double _Brent2 = 0.001;
    private boolean _dentroscropeOutput = false;
    private int _parallel = 1;
    private Set<String> _fixedHybrid = new HashSet<String>();
    private double[] _operationWeight = {0.1,0.1,0.15,0.55,0.15,0.15,2.8};
    private boolean _postOptimization = false;
    private boolean _optimizeBL = false;
    private int _numRuns = 10;
    private Long _seed = null;
    private boolean _oneGTPerLocus = true;



    public InferNetwork_MPL(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                            Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                            Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader, Random rand){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams(){
        return 2;
    }

    @Override
    protected int getMaxNumParams(){
        return 40;
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
                ParamExtractorAllelListMap aaParam = new ParamExtractorAllelListMap("a", this.params, this.errorDetected);
                noError = noError && aaParam.IsValidMap;
                if(aaParam.IsValidMap){
                    _taxonMap = aaParam.ValueMap;
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
                        errorDetected.execute("Unrecognized maximum number of networks for search " + plParam.PostSwitchValue, plParam.PostSwitchParam.getLine(), plParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -pl.", plParam.SwitchParam.getLine(), plParam.SwitchParam.getColumn());
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
            //It's for hill climbing
/*
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
*/

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

            ParamExtractor iParam = new ParamExtractor("i", this.params, this.errorDetected);
            if(iParam.ContainsSwitch)
            {
                if(iParam.PostSwitchParam != null)
                {
                    try
                    {
                        _improvementThreshold = Double.parseDouble(iParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value after -i " + iParam.PostSwitchValue, iParam.PostSwitchParam.getLine(), iParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -i.", iParam.SwitchParam.getLine(), iParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor rParam = new ParamExtractor("r", this.params, this.errorDetected);
            if(rParam.ContainsSwitch)
            {
                if(rParam.PostSwitchParam != null)
                {
                    try
                    {
                        _maxRounds = Integer.parseInt(rParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value after -r " + rParam.PostSwitchValue, rParam.PostSwitchParam.getLine(), rParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -r.", rParam.SwitchParam.getLine(), rParam.SwitchParam.getColumn());
                }
            }

            /*
            ParamExtractor htParam = new ParamExtractor("ht", this.params, this.errorDetected);
            if(htParam.ContainsSwitch)
            {
                if(htParam.PostSwitchParam != null)
                {
                    try
                    {
                        _hasTried = Integer.parseInt(htParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value of maximum branch length " + htParam.PostSwitchValue, htParam.PostSwitchParam.getLine(), htParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -l.", htParam.SwitchParam.getLine(), htParam.SwitchParam.getColumn());
                }
            }
            */

            ParamExtractor lParam = new ParamExtractor("l", this.params, this.errorDetected);
            if(lParam.ContainsSwitch)
            {
                if(lParam.PostSwitchParam != null)
                {
                    try
                    {
                        _maxBranchLength = Double.parseDouble(lParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value of maximum branch length " + lParam.PostSwitchValue, lParam.PostSwitchParam.getLine(), lParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -l.", lParam.SwitchParam.getLine(), lParam.SwitchParam.getColumn());
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

            ParamExtractor poParam = new ParamExtractor("po", this.params, this.errorDetected);
            if(poParam.ContainsSwitch){
                if(poParam.PostSwitchParam != null && !poParam.PostSwitchValue.startsWith("-"))
                {
                    errorDetected.execute("No value expected after switch -po.", poParam.SwitchParam.getLine(), poParam.SwitchParam.getColumn());
                }
                else
                {
                    _postOptimization = true;
                }
            }


            ParamExtractor oParam = new ParamExtractor("o", this.params, this.errorDetected);
            if(oParam.ContainsSwitch){
                if(oParam.PostSwitchParam != null && !oParam.PostSwitchValue.startsWith("-"))
                {
                    errorDetected.execute("No value expected after switch -o.", oParam.SwitchParam.getLine(), oParam.SwitchParam.getColumn());
                }
                else
                {
                    _optimizeBL = true;
                }
            }



            ParamExtractor pParam = new ParamExtractor("p", this.params, this.errorDetected);
            if(pParam.ContainsSwitch)
            {
                if(pParam.PostSwitchParam != null)
                {
                    try
                    {
                        if(!(pParam.PostSwitchParam instanceof ParameterIdentList)){
                            throw new RuntimeException();
                        }
                        ParameterIdentList brent = (ParameterIdentList)pParam.PostSwitchParam;
                        int index = 0;
                        for(String value: brent.Elements){
                            if(index==0){
                                _Brent1 = Double.parseDouble(value.trim());
                            }else if(index==1){
                                _Brent2 = Double.parseDouble(value.trim());
                            }
                            else{
                                throw new RuntimeException();
                            }
                            index++;
                        }
                        if(_Brent1==0 || _Brent2==0){
                            throw new RuntimeException();
                        }
                        //ParameterIdentList brent = this.assertParameterIdentList(pParam.PostSwitchValue);
                        /*
                        String brentPrecision = pParam.PostSwitchValue;
                        if(!brentPrecision.contains("~")){
                            throw new NumberFormatException();
                        }
                        String[] splitValues = brentPrecision.split("~");

                        _Brent1 = Double.parseDouble(splitValues[0].trim());
                        _Brent2 = Double.parseDouble(splitValues[1].trim());
                        assert(_Brent1 > _Brent2);
                        */
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value after switch -p.", pParam.PostSwitchParam.getLine(), pParam.PostSwitchParam.getColumn());
                    }

                }
                else
                {
                    errorDetected.execute("Expected value after switch -p.", pParam.SwitchParam.getLine(), pParam.SwitchParam.getColumn());
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
                        for(String wExp: weights.Elements){
                            if(index < _operationWeight.length){
                                double w = Double.parseDouble(wExp.trim());
                                if(w < 0){
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

            ParamExtractor tParam = new ParamExtractor("t", this.params, this.errorDetected);
            if(tParam.ContainsSwitch)
            {
                if(tParam.PostSwitchParam != null)
                {
                    try
                    {
                        _maxTryPerBranch = Integer.parseInt(tParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized value after -t " + tParam.PostSwitchValue, tParam.PostSwitchParam.getLine(), tParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -t.", tParam.SwitchParam.getLine(), tParam.SwitchParam.getColumn());
                }
            }

            ParamExtractor diParam = new ParamExtractor("di", this.params, this.errorDetected);
            if(diParam.ContainsSwitch)
            {
                _dentroscropeOutput = true;
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

            noError = noError && checkForUnknownSwitches("a","b","s","m","n","md","rd","p","l","r","i","t","di","pl","h","w","x","rs","o","po");
            checkAndSetOutFile(aParam, bParam, sParam, mParam, nParam, mdParam, rdParam, pParam, lParam, rParam, iParam,tParam, diParam,plParam,hParam,wParam,xParam,rsParam,oParam,poParam);

        }



        return  noError;
    }

    @Override
    protected String produceResult() {
        StringBuffer result = new StringBuffer();
        List<List<MutableTuple<Tree,Double>>> gts = new ArrayList<List<MutableTuple<Tree,Double>>>();
        //List<Tree> gts2 = new ArrayList<Tree>();
        //List<Integer> counter = new ArrayList<Integer>();
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
                if (_bootstrap < 100) {
                    if (Trees.handleBootStrapInTree(newtr, _bootstrap) == -1) {
                        throw new IllegalArgumentException("Input gene tree " + newtr + " have nodes that don't have bootstrap value");
                    }

                }

                gtsForOneLocus.add(new MutableTuple<Tree, Double>(newtr, weight));
            }
            gts.add(gtsForOneLocus);
            //((STINode<Double>)newtr.getRoot()).setData(weight);
            //gts2.add(newtr);
        }

        /*
        if(_fixedHybrid.size()!=0){
            for(int i=0; i<_operationWeight.length-1; i++){
                if(i == 3){
                    _operationWeight[i] = 1.0;
                }
                else{
                    _operationWeight[i] = 0;
                }
            }
            _operationWeight[_operationWeight.length-1] = 1/0.3*0.7;
        }
        */

        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        Network speciesNetwork = null;
        if(_startSpeciesNetwork!=null){
            speciesNetwork = transformer.makeNetwork(_startSpeciesNetwork);
        }


        InferNetworkMLFromGT inference;
        if(_oneGTPerLocus)
            inference = new InferNetworkPseudoMLFromGTT_SingleTreePerLocus();
        else
            inference = new InferNetworkPseudoMLFromGTT_MultiTreesPerLocus();

        inference.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength, _Brent1, _Brent2, _maxExaminations, _maxFailure, _moveDiameter, _reticulationDiameter, _parallel, speciesNetwork, _fixedHybrid, _operationWeight, _numRuns, _optimizeBL, _seed);
        LinkedList<Tuple<Network, Double>> resultTuples = new LinkedList<>();
        inference.inferNetwork(gts,_taxonMap, _maxReticulations, _returnNetworks, _postOptimization, resultTuples);


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
            result.append("\n" + "Total log probability: " + tuple.Item2);

            if(_dentroscropeOutput){
                edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.removeAllParameters(n);
                result.append("\nVisualize in Dendroscope : " + n.toString());
            }
        }

        return result.toString();

    }





}
