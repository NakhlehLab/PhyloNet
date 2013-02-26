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

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.ilp.IlpGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/8/11
 * Time: 2:02 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("gencplex")
public class GenCPLEX extends CommandBase  {

    private Iterable<NetworkNonEmpty> _syntaxTrees;

    private Iterable<NetworkNonEmpty> _geneTrees;

    private Double _deepCoalWeightParam;

    private Double _numCoalCount;

    private File _outDir;

    public GenCPLEX(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                    Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 5;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 5;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = true;

        ParameterIdentSet syntaxTreesParam = this.assertParameterIdentSet(0);
        noError = noError && syntaxTreesParam != null;

        if(noError)
        {
            _syntaxTrees = this.assertTreesExist(syntaxTreesParam);
        }

        ParameterIdentSet geneTreesParam = this.assertParameterIdentSet(1);
        noError = noError && geneTreesParam != null;

        if(noError)
        {
            _geneTrees = this.assertTreesExist(geneTreesParam);
        }

        ParameterIdent deepCoalWeightParam = this.assertParameterIdent(2);
        noError = noError && deepCoalWeightParam != null;

        if(deepCoalWeightParam != null)
        {
            try
            {
                _deepCoalWeightParam = new Double(deepCoalWeightParam.execute(GetSimpleParamValue.Singleton, null));
            }
            catch(NumberFormatException e)
            {
                noError = false;
                errorDetected.execute("Unrecognized number.", deepCoalWeightParam.getLine(), deepCoalWeightParam.getColumn());

            }
        }

        ParameterIdent numCoalCountParam = this.assertParameterIdent(3);
        noError = noError && numCoalCountParam != null;

        if(numCoalCountParam != null)
        {
            try
            {
                _numCoalCount = new Double(numCoalCountParam.execute(GetSimpleParamValue.Singleton, null));
            }
            catch(NumberFormatException e)
            {
                noError = false;
                errorDetected.execute("Unrecognized number.", numCoalCountParam.getLine(), numCoalCountParam.getColumn());

            }
        }

        ParameterQuote outDirectoryParam = this.assertQuotedParameter(4);
        noError = noError && outDirectoryParam != null;

        if(outDirectoryParam != null)
        {
            String outDirValue = outDirectoryParam.execute(GetSimpleParamValue.Singleton, null);
            _outDir = new File(outDirValue);

            if(!_outDir.isDirectory() || !_outDir.exists())
            {
                noError = false;
                this.errorDetected.execute("Invalid directory: '" + _outDir.getAbsolutePath() + "'.", outDirectoryParam.getLine(), outDirectoryParam.getColumn());
                _outDir = null;
            }

        }
        return noError;
    }

    @Override
    protected void executeCommandHelp(Proc1<String> displayResult) throws IOException
    {
        if(_deepCoalWeightParam == null || _geneTrees == null || _numCoalCount == null || _outDir == null || _syntaxTrees == null)
        {
            throw new IllegalStateException();
        }

        List<STITree<Integer>> speciesTrees = new LinkedList<STITree<Integer>>();
		List<Tree> geneTrees = new LinkedList<Tree>();

		try {
			String line;
			for(NetworkNonEmpty syntaxTree : _syntaxTrees)
            {
				NewickReader nw = new NewickReader(new StringReader(NetworkTransformer.toENewick(syntaxTree)));
				STITree<Integer> st = new STITree<Integer>(nw.readTree());

				speciesTrees.add(st);
			}

			for(NetworkNonEmpty geneTree : _geneTrees)
            {
					NewickReader nw = new NewickReader(new StringReader(NetworkTransformer.toENewick(geneTree)));
					geneTrees.add(nw.readTree());
			}

		}
		catch (Exception e) {
			throw new RuntimeException(e);
		}

		for (int i = 0; i < speciesTrees.size(); i++) {
			String cplexFileName = _outDir.getAbsolutePath() + File.separatorChar + "input" + i;
			String variableFileName = _outDir.getAbsolutePath() + File.separatorChar + "var" + i;
			String scriptFileName = _outDir.getAbsolutePath() + File.separatorChar + "script" + i;

			IlpGenerator builder = new IlpGenerator();

			builder.setSfWeight(_deepCoalWeightParam);
			builder.setSgWeight(_numCoalCount);
			builder.generateCplexInput(speciesTrees.get(i), geneTrees, cplexFileName, variableFileName, scriptFileName);
		}

    }
}
