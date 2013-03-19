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
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/9/11
 * Time: 4:06 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class InferSTBase extends CommandBaseFileOut
{
    class ThresholdResult
    {
        public final boolean NoError;

        public final double Threshold;

        public final ParamExtractor Extractor;

        public ThresholdResult(boolean  noError, double threshold, ParamExtractor extractor)
        {
            NoError = noError;
            Threshold = threshold;
            Extractor = extractor;
        }
    }

    class TaxonMapResult
    {
        public final boolean NoError;

        public final Map<String,String> TaxonMap;

        public final ParamExtractor Extractor;

        public TaxonMapResult(boolean  noError, Map<String,String> taxonMap, ParamExtractor extractor)
        {
            NoError = noError;
            TaxonMap = taxonMap;
            Extractor = extractor;
        }
    }

    protected Iterable<NetworkNonEmpty> _geneTrees;



    InferSTBase(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected abstract int getMinNumParams();

    @Override
    protected abstract int getMaxNumParams();

    @Override
    protected boolean checkParamsForCommand()
    {
       boolean noError = true;

        ParameterIdentList geneTreesParam = this.assertParameterIdentList(0);
        noError = noError && geneTreesParam != null;

        if(noError)
        {
            _geneTrees = this.assertNetworksExist(geneTreesParam);
        }

        return noError;
    }

     @Override
    protected abstract String produceResult();

    protected List<Tree> GetGeneTreesAsSTIDoubleTreeList()
    {
        List<Tree> trees = new ArrayList<Tree>();

        for(NetworkNonEmpty geneTree : _geneTrees)
        {
            try
            {
                NewickReader nr = new NewickReader(new StringReader(NetworkTransformer.toENewick(geneTree)));
			    STITree<Double> gt = new STITree<Double>(true);
		        nr.readTree(gt);
			    trees.add(gt);
            }
            catch(Exception e)
            {
                throw new RuntimeException(e);
            }
        }

        return trees;
    }

     protected List<Tree> GetGeneTreesAsTreeList()
    {
        List<Tree> trees = new ArrayList<Tree>();

        for(NetworkNonEmpty geneTree : _geneTrees)
        {
            try
            {
                NewickReader nr = new NewickReader(new StringReader(NetworkTransformer.toENewick(geneTree)));
			    Tree tr = nr.readTree();
				trees.add(tr);
            }
            catch(Exception e)
            {
                throw new RuntimeException(e);
            }
        }

        return trees;
    }

    protected ThresholdResult assignThreshold(double defaultBootstrap)
    {
        double bootstrap = defaultBootstrap;
        boolean noError = true;

        ParamExtractor bParam = new ParamExtractor("b", this.params, this.errorDetected);
        if(bParam.ContainsSwitch)
        {
            if(bParam.PostSwitchParam != null)
            {
                try
                {
                    bootstrap = Double.parseDouble(bParam.PostSwitchValue);
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unknown bootstrap '" + bParam.PostSwitchValue + "'.", bParam.PostSwitchParam.getLine(), bParam.PostSwitchParam.getColumn());
                    noError = false;
                }
            }
            else
            {
                errorDetected.execute("Expected value after bootstrap switch.", bParam.SwitchParam.getLine(), bParam.SwitchParam.getColumn());
                noError = false;
            }
        }

        if(noError && bParam.ContainsSwitch && _geneTrees != null)
        {
            for(NetworkNonEmpty gt : _geneTrees)
            {
                noError = noError && assignThresholdHelp(gt.PrincipleInfo, gt.PrincipleDescendants, true);
            }
        }

        return new ThresholdResult(noError, bootstrap, bParam);
    }

    private boolean assignThresholdHelp(final NetworkInfo info, final DescendantList descendantList, final boolean isRoot)
    {
        if(info.Support.execute(new SupportAlgo<Boolean, Object, RuntimeException>()
        {
            public Boolean forSupportNonEmpty(SupportNonEmpty supportNonEmpty, Object o) throws RuntimeException {
                return false;
            }

            public Boolean forSupportEmpty(SupportEmpty supportEmpty, Object o) throws RuntimeException {
                return true && !isRoot && descendantList.Subtrees.iterator().hasNext();  //To change body of implemented methods use File | Settings | File Templates.
            }
        }, null))
        {
            errorDetected.execute(
                    String.format("If bootstrap switch is specified for '%s', all internal gene tree nodes must have a support value.", _motivatingCommand.getName()),
                                  _motivatingCommand.getLine(), _motivatingCommand.getColumn());
            return false;
        }

        for(Subtree childTree : descendantList.Subtrees)
        {
            if(!assignThresholdHelp(childTree.NetworkInfo, childTree.Descendants, false))
            {
                return false;
            }
        }

        return true;
    }

    protected TaxonMapResult assignTaxonMap()
    {
        ParamExtractorAllelMap aParam = new ParamExtractorAllelMap("a", this.params, this.errorDetected);
        boolean noError = true;
        HashMap<String,String> taxonMap = null;

        if(aParam.ContainsSwitch)
        {
            noError = noError && aParam.IsValidMap;

            if(aParam.IsValidMap)
            {
                taxonMap = aParam.ValueMap;
            }
        }

        return new TaxonMapResult(noError, taxonMap, aParam);
    }




}
