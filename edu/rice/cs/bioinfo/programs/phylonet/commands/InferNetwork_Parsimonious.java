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
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferILSNetworkParsimoniously;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
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
@CommandName("infernetwork_parsimony")
public class InferNetwork_Parsimonious extends CommandBaseFileOut{
    private HashMap<String, List<String>> _taxonMap = null;
    private List<NetworkNonEmpty> _geneTrees;
    private double _bootstrap = 100;
    private NetworkNonEmpty _startSpeciesNetwork = null;
    private Long _maxReticulations;
    private Long _maxExaminations = null;
    private int _maxDiameter = 0;
    private int _returnNetworks = 1;
    private  boolean _dentroscropeOutput = false;


    public InferNetwork_Parsimonious(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                                     Map<String,NetworkNonEmpty>  sourceIdentToNetwork,
                                     Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams(){
        return 2;
    }

    @Override
    protected int getMaxNumParams(){
        return 13;
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
            _maxReticulations = new Long(Integer.parseInt(number.Content));
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

            ParamExtractor diParam = new ParamExtractor("di", this.params, this.errorDetected);
            if(diParam.ContainsSwitch)
            {
                _dentroscropeOutput = true;
            }

            noError = noError && checkForUnknownSwitches("a","b","s","n", "m", "d", "di");
            checkAndSetOutFile(aParam, bParam, sParam, nParam, mParam, dParam, diParam);
        }



        return  noError;
    }

    @Override
    protected String produceResult() {
        StringWriter result = new StringWriter();

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
        InferILSNetworkParsimoniously inference = new InferILSNetworkParsimoniously();
        List<Tuple<Network, Double>> resultTuples = inference.inferNetwork(gts,_taxonMap,_maxExaminations,_maxReticulations,_maxDiameter,speciesNetwork, _returnNetworks);
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

            StringWriter writer = new StringWriter();
            RnNewickPrinter printer = new RnNewickPrinter();
            printer.print(tuple.Item1, writer);
            result.append("\n" + writer.toString());
            result.append("\n" + "Total number of extra lineages: " + tuple.Item2);

            if(_dentroscropeOutput){
                for(Object node : n.getNetworkNodes())
                {
                    NetNode netNode = (NetNode)node;
                    for(Object parent: netNode.getParents())
                    {
                        NetNode parentNode = (NetNode)parent;
                        netNode.setParentProbability(parentNode, Double.NaN);
                    }
                }
                writer = new StringWriter();
                printer = new RnNewickPrinter();
                printer.print(tuple.Item1, writer);
                result.append("\nVisualize in Dendroscope : " + writer.toString());
            }
        }

        return result.toString();

    }
}
