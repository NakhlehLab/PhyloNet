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
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterTaxaMap;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.GLASSInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.TaxaDistanceMatrix;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/14/11
 * Time: 1:08 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("Infer_ST_GLASS")
public class InfterST_MDC_GLASS extends InferSTBase
{
    private Map<String,String> _taxonMap;

    private TaxaDistanceMatrix _distanceMatrix;

    public InfterST_MDC_GLASS(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                              Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 6;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand()
    {
       boolean noError = true;

       Parameter firstParam = this.params.get(0);
       if(firstParam instanceof ParameterIdentList)
       {
           _geneTrees = this.assertNetworksExist((ParameterIdentList)firstParam);
       }
       else if(firstParam instanceof ParameterTaxaMap)
       {
           ParameterTaxaMap distanceMap = (ParameterTaxaMap)firstParam;

           HashSet<String> taxaInMap = new HashSet<String>();

           for(Map.Entry<String, List<String>> entry : distanceMap._mappings)
           {
               taxaInMap.add(entry.getKey());
               List<String> taxonDistancePairs = entry.getValue();
               for(int i = 0; i<taxonDistancePairs.size(); i+=2)
               {
                   taxaInMap.add(taxonDistancePairs.get(i));
               }
           }
           String[] taxa = taxaInMap.toArray(new String[0]);
           _distanceMatrix = new TaxaDistanceMatrix(taxa);


           for(Map.Entry<String, List<String>> entry : distanceMap._mappings)
           {
               String taxon1 = entry.getKey();
               List<String> taxonDistancePairs = entry.getValue();

               if(taxonDistancePairs.size() % 2 != 0)
               {
                    noError = false;
                    this.errorDetected.execute("Odd number of taxon/distance entries found in the entries for taxon '" + taxon1 + "'.",
                                               distanceMap.getLine(), distanceMap.getColumn());
               }

               HashSet<String> seenMapKeys = new HashSet<String>();
               for(int i = 0; i<taxonDistancePairs.size(); i+=2)
               {
                   String taxon2 = taxonDistancePairs.get(i);
                   String stringDistance = taxonDistancePairs.get(i+1);

                   try
                   {
                        double distance = Double.parseDouble(stringDistance);
                        String mapKey1 = taxon1.toLowerCase() + taxon2.toLowerCase();
                        String mapKey2 = taxon2.toLowerCase() + taxon1.toLowerCase();

                       if(!seenMapKeys.contains(mapKey1) && !seenMapKeys.contains(mapKey2))
                       {
                           seenMapKeys.add(mapKey1);
                           seenMapKeys.add(mapKey2);

                           STITreeCluster c = new STITreeCluster(taxa);
				           c.addLeaf(taxon1);
				           c.addLeaf(taxon2);
				           _distanceMatrix.put(c, distance);
                       }
                       else
                       {
                           this.errorDetected.execute(
                                   String.format("Duplicate mapping between taxa '%s' and '%s'.", taxon1, taxon2),
                                   distanceMap.getLine(), distanceMap.getColumn());
                           noError = false;
                       }


                   }
                   catch(NumberFormatException e)
                   {
                       noError = false;
                       this.errorDetected.execute("Unknown number '" + stringDistance + "'.", distanceMap.getLine(), distanceMap.getColumn());
                   }
               }



           }


       }
       else
       {
           noError = false;
           this.errorDetected.execute("Expected first parameter to command InfterST_MDC_GLASS to be an identifier list or distance map.",
                                       firstParam.getLine(), firstParam.getColumn());
       }

        TaxonMapResult result = assignTaxonMap();
        noError = noError && result.NoError;
        _taxonMap = result.TaxonMap;

         noError = noError && checkForUnknownSwitches("a");

        if(!noError)
        {
            _geneTrees = null;
            _distanceMatrix = null;
        }

        return noError;
    }

    @Override
    protected String produceResult() {
        if(_geneTrees == null && _distanceMatrix == null)
        {
            throw new IllegalStateException();
        }

        StringBuffer result = new StringBuffer();



        GLASSInference inference = new GLASSInference();
		Tree inferredTree;

        if(_distanceMatrix != null)
        {
            inferredTree = inference.inferSpeciesTreeFromTaxa(_distanceMatrix);
        }
        else if(_taxonMap == null)
        {
                List<Tree> trees = GetGeneTreesAsTreeList();
				inferredTree = inference.inferSpeciesTree(trees);
	    }
        else
        {
                List<Tree> trees = GetGeneTreesAsTreeList();
				inferredTree = inference.inferSpeciesTree(trees, _taxonMap);
		}

        String tree = inferredTree.toString();
        this.richNewickGenerated(tree);
        result.append("\n" + tree);

        return result.toString();
    }
}
