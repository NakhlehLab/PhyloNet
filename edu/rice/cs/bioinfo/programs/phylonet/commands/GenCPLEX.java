package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc;
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
public class GenCPLEX extends CommandBase  {

    private Iterable<NetworkNonEmpty> _syntaxTrees;

    private Iterable<NetworkNonEmpty> _geneTrees;

    private Double _deepCoalWeightParam;

    private Double _numCoalCount;

    private File _outDir;

    GenCPLEX(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
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
            _syntaxTrees = this.assertNetworksExist(syntaxTreesParam);
        }

        ParameterIdentSet geneTreesParam = this.assertParameterIdentSet(1);
        noError = noError && geneTreesParam != null;

        if(noError)
        {
            _geneTrees = this.assertNetworksExist(geneTreesParam);
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
    protected void executeCommandHelp(Proc<String> displayResult) throws IOException
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
