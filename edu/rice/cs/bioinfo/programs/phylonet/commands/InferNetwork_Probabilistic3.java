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

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferILSNetworkProbabilistically;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferILSNetworkProbabilistically3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferILSNetworkUsingBLProbabilistically;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferILSNetworkUsingBLProbabilistically3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
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
 * Created with IntelliJ IDEA.
 * User: jianrongdong
 * Date: 3/18/13
 * Time: 3:41 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("infernetwork_ML3")
public class InferNetwork_Probabilistic3 extends CommandBaseFileOut{
    private HashMap<String, List<String>> _taxonMap = null;
    private List<NetworkNonEmpty> _geneTrees;
    private double _bootstrap = 100;
    private NetworkNonEmpty _startSpeciesNetwork = null;
    private int _maxReticulations;
    private Long _maxExaminations = null;
    private int _maxDiameter = 0;
    private int _returnNetworks = 1;
    private int _maxRounds = 100;
    private int _maxTryPerBranch = 100;
    private double _maxBranchLength = 6;
    private double _improvementThreshold = 0.001;

    private double _Brent1 = 0.01;
    private double _Brent2 = 0.001;
    private boolean  _dentroscropeOutput = false;
    private boolean _usingBL = false;

    public InferNetwork_Probabilistic3(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
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
        return 28;
    }

    @Override
    protected boolean checkParamsForCommand(){
        boolean noError = true;

        ParameterIdentList geneTreeParam = this.assertParameterIdentList(0);
        noError = noError && geneTreeParam != null;
        _geneTrees = new LinkedList<NetworkNonEmpty>();

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

            for(String ident : geneTreeParam.Elements)
            {
                noError = noError && this.assertNetworkExists(ident, geneTreeParam.getLine(), geneTreeParam.getColumn());
                if(noError)
                {
                    _geneTrees.add(this.sourceIdentToNetwork.get(ident));
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

            ParamExtractor dParam = new ParamExtractor("d", this.params, this.errorDetected);
            if(dParam.ContainsSwitch){
                if(dParam.PostSwitchParam != null)
                {
                    try
                    {
                        _maxDiameter = Integer.parseInt(dParam.PostSwitchValue);
                    }
                    catch(NumberFormatException e)
                    {
                        errorDetected.execute("Unrecognized maximum diameter for network search " + dParam.PostSwitchValue, dParam.PostSwitchParam.getLine(), dParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -d.", dParam.SwitchParam.getLine(), dParam.SwitchParam.getColumn());
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
                    catch(RuntimeException e)
                    {
                        errorDetected.execute("Unexpected value after switch -p.", pParam.PostSwitchParam.getLine(), pParam.PostSwitchParam.getColumn());
                    }
                }
                else
                {
                    errorDetected.execute("Expected value after switch -p.", pParam.SwitchParam.getLine(), pParam.SwitchParam.getColumn());
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

            ParamExtractor blParam = new ParamExtractor("bl", this.params, this.errorDetected);
            if(blParam.ContainsSwitch)
            {
                _usingBL = true;
            }

            noError = noError && checkForUnknownSwitches("a","b","s","m","n","d","p","l","r","i","t","di","bl");
            checkAndSetOutFile(aParam, bParam, sParam, mParam, nParam, dParam, pParam, lParam, rParam, iParam,tParam,diParam,blParam);
        }

        return  noError;
    }

    @Override
    protected String produceResult() {
        StringBuffer result = new StringBuffer();

        List<Tree> gts = new ArrayList<Tree>();
        //List<Integer> counter = new ArrayList<Integer>();
        for(NetworkNonEmpty geneTree : _geneTrees){

            double prob = geneTree.TreeProbability.execute(new TreeProbabilityAlgo<Double, RuntimeException>() {
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
            try
            {
                nr.readTree(newtr);
            }
            catch(Exception e)
            {
                errorDetected.execute(e.getMessage(),
                        this._motivatingCommand.getLine(), this._motivatingCommand.getColumn());
            }
            Trees.removeBinaryNodes(newtr);
            if(_bootstrap<100){
                if(Trees.handleBootStrapInTree(newtr, _bootstrap)==-1){
                    throw new IllegalArgumentException("Input gene tree " + newtr + " have nodes that don't have bootstrap value");
                }

            }
            newtr.getRoot().setData(prob);
            gts.add(newtr);
        }

        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        Network speciesNetwork = null;
        if(_startSpeciesNetwork!=null){
            speciesNetwork = transformer.makeNetwork(_startSpeciesNetwork);
        }

        //long start = System.currentTimeMillis();
        if (!_usingBL){   // the old no branch length route
            InferILSNetworkProbabilistically3 inference = new InferILSNetworkProbabilistically3();
            inference.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength,
                                         _Brent1, _Brent2, _maxExaminations, _maxDiameter, speciesNetwork);
            int rightK = inference.CV(gts,_taxonMap,_maxReticulations,_returnNetworks);
            return Integer.toString(rightK);
        }
        else{ // with new branch length route
            InferILSNetworkUsingBLProbabilistically3 inference = new InferILSNetworkUsingBLProbabilistically3();
            inference.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength,
                                         _Brent1, _Brent2, _maxExaminations, _maxDiameter, speciesNetwork);
            int rightK = inference.CV(gts,_taxonMap,_maxReticulations, _returnNetworks);
            return Integer.toString(rightK);
        }
        //System.out.print(System.currentTimeMillis()-start);
    }
}

