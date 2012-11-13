/*
 * Copyright (c) 2012 Rice University.
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
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentSet;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.DeepCoalescencesCounter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.StringReader;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 4:45 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("deepcoalcount_tree")
public class DeepCoalCount extends CommandBaseFileOut {

    private HashMap<String,String> _taxonMap = null;

    private boolean  _treatGeneTreesAsRooted = true;

    private double _bootstrap = 1.0;

    private List<NetworkNonEmpty> _speciesTrees;

    private List<NetworkNonEmpty> _geneTrees;

    public DeepCoalCount(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    public int getMinNumParams()
    {
        return 2;
    }

    public int getMaxNumParams()
    {
        return 8;
    }

    public boolean checkParamsForCommand() {

        boolean noError = true;

        ParameterIdentSet speciesTreeParam = this.assertParameterIdentSet(0);
        noError = noError && speciesTreeParam != null;

        ParameterIdentSet geneTreeParam    = this.assertParameterIdentSet(1);
        noError = noError && geneTreeParam != null;

        ParamExtractor uParam = new ParamExtractor("u", this.params, this.errorDetected);
        if(uParam.ContainsSwitch)
        {
            _treatGeneTreesAsRooted = false;
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


        ParamExtractorAllelMap aParam = new ParamExtractorAllelMap("a", this.params, this.errorDetected);

        if(aParam.ContainsSwitch)
        {
            noError = noError && aParam.IsValidMap;

            if(aParam.IsValidMap)
            {
                _taxonMap = aParam.ValueMap;
            }
        }

        noError = noError && checkForUnknownSwitches("u", "b", "a");
        checkAndSetOutFile(bParam, aParam);

        if(noError)
        {
            return checkContext(speciesTreeParam, geneTreeParam);
        }
        else
        {
            return noError;
        }

    }

    private boolean checkContext(ParameterIdentSet speciesTreeParam, ParameterIdentSet geneTreeParam ) {

        boolean noError = true;

        _speciesTrees = new LinkedList<NetworkNonEmpty>();
        for(String ident : speciesTreeParam.Elements)
        {
            noError = noError && this.assertNetworkExists(ident, speciesTreeParam.getLine(), speciesTreeParam.getColumn());
            if(noError)
            {
                _speciesTrees.add(this.sourceIdentToNetwork.get(ident));
            }
        }

        _geneTrees = new LinkedList<NetworkNonEmpty>();
        for(String ident : geneTreeParam.Elements)
        {
            noError = noError && this.assertNetworkExists(ident, geneTreeParam.getLine(), geneTreeParam.getColumn());
            if(noError)
            {
                _geneTrees.add(this.sourceIdentToNetwork.get(ident));
            }
        }

        if(!noError)
        {
            _speciesTrees = null;
            _geneTrees = null;

        }

        return noError;

    }

    @Override
    protected String produceResult() {

        StringBuffer result = new StringBuffer();


           List<Tree> geneTrees = new LinkedList<Tree>();

           for(NetworkNonEmpty geneTree : _geneTrees)
           {
               String phylonetGeneTree = NetworkTransformer.toENewickTree(geneTree);
               NewickReader nr = new NewickReader(new StringReader(phylonetGeneTree));
			   STITree<Double> gt = new STITree<Double>(true);
			   try
               {
                   nr.readTree(gt);
                   geneTrees.add(gt);
               }
               catch(Exception e)
               {
                    errorDetected.execute(e.getMessage(), -1, -1);
               }

           }

           int index = 1;
           for(NetworkNonEmpty st : _speciesTrees)
           {
               String phylonetSpeciesTree = NetworkTransformer.toENewickTree(st);
               NewickReader nr = new NewickReader(new StringReader(phylonetSpeciesTree));

               try
               {
                    Tree speciesTree = nr.readTree();
                    int coalNum;

                    if(_taxonMap == null)
                    {
                       coalNum = DeepCoalescencesCounter.countExtraCoal(geneTrees, speciesTree, _treatGeneTreesAsRooted, _bootstrap);
                    }
                    else
                    {
                       coalNum = DeepCoalescencesCounter.countExtraCoal(geneTrees, speciesTree, _taxonMap, _treatGeneTreesAsRooted, _bootstrap);
                    }
                    String speciesTreeString = speciesTree.toStringWD();
                    result.append("\nSpecies_Tree#" + (index++ ) + " = " + speciesTreeString + "\n");
			        result.append("Total number of extra lineages: " + coalNum);
                    this.treeGenerated(speciesTreeString);
               }
               catch(Exception e)
               {
                   errorDetected.execute(e.getMessage(), -1, -1);
               }
           }


        return result.toString();

    }
}
