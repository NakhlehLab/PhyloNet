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
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.GeneTreeRefinement;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.StringReader;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/26/11
 * Time: 3:10 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("processgt")
public class ProcessGT extends CommandBaseFileOut
{
    LinkedList<NetworkNonEmpty> _speciesTrees = new LinkedList<NetworkNonEmpty>();

    LinkedList<NetworkNonEmpty> _geneTrees = new LinkedList<NetworkNonEmpty>();

    private  double _bootstrap = 100;

    private HashMap<String, List<String>> _taxonMap = null;

    private boolean _geneTreesRooted = true;

    public ProcessGT(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
              Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    protected int getMinNumParams()
    {
        return 2;
    }

    protected int getMaxNumParams()
    {
        return 8;
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = true;


        ParameterIdentSet speciesTreesIdents = super.assertParameterIdentSet(0);
        if(speciesTreesIdents != null)
        {
            _speciesTrees = assertNetworksExist(speciesTreesIdents);
            noError = _speciesTrees != null && noError;
        }

        ParameterIdentSet geneTreesIdents = super.assertParameterIdentSet(1);
        if(geneTreesIdents != null)
        {
            _geneTrees = assertNetworksExist(geneTreesIdents);
            noError = _geneTrees != null && noError;
        }

         ParamExtractor uParam = new ParamExtractor("u", this.params, this.errorDetected);
        if(uParam.ContainsSwitch)
        {
            _geneTreesRooted = false;
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


        ParamExtractorAllelListMap aParam = new ParamExtractorAllelListMap("a", this.params, this.errorDetected);

        if(aParam.ContainsSwitch)
        {
            noError = noError && aParam.IsValidMap;

            if(aParam.IsValidMap)
            {
                _taxonMap = aParam.ValueMap;
            }
        }

        checkAndSetOutFile(bParam, aParam);

        return noError;
    }

    @Override
    protected String produceResult()
    {
        List<Tree> speciesTrees;
		List<Tree> geneTrees = new ArrayList<Tree>();
		//Map<String,List<String>> taxonMap = null;

        for(NetworkNonEmpty geneTree : _geneTrees)
        {

            String eNewickString = NetworkTransformer.toENewick(geneTree);
		    NewickReader nr = new NewickReader(new StringReader(eNewickString));
            STITree<Double> gt = new STITree<Double>(true);
            try
            {
			    nr.readTree(gt);
            }
            catch(Exception e)
            {
                throw new RuntimeException(e);
            }
            geneTrees.add(gt);
		}

        speciesTrees = new ArrayList<Tree>();
        for(NetworkNonEmpty speciesTree : _speciesTrees)
        {
            String eNewickString = NetworkTransformer.toENewick(speciesTree);
            NewickReader nr = new NewickReader(new StringReader(eNewickString));
            try
            {
			    speciesTrees.add(nr.readTree());
            }
            catch(Exception e)
            {
                throw new RuntimeException(e);
            }
		}

        StringBuffer result = new StringBuffer();
        int index = 1;

		for(Tree st: speciesTrees){
			GeneTreeRefinement.processGeneTrees(geneTrees, st, _taxonMap, _geneTreesRooted, _bootstrap);
            String stString = st.toString();
            this.richNewickGenerated(stString);
			result.append("\nSpecies_Tree#" + (index++ ) + " = " + stString + "\n");
			result.append("Resulting gene trees:");
			for(Tree gt: geneTrees){
                String geneTreeString = gt.toString();
				result.append("\n"+geneTreeString);
                this.richNewickGenerated(geneTreeString);
			}
		}

        return result.toString();


    }

}
