package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentSet;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.ArrayList;
import java.util.Map;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/29/11
 * Time: 11:08 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferST_Bootstrap extends CommandBaseFileOut {

    private int _numRepititions;

    private double _threshold = 0.5;

    private Random _rand;

    private  SyntaxCommand _exampleCommand;

    InferST_Bootstrap(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected, Random rand) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
        _rand = rand;
    }

    @Override
    protected int getMinNumParams() {
        return 3;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 16;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = true;

        ParameterIdent repititions = this.assertParameterIdent(0);
        noError = noError && repititions != null;

        if(repititions != null)
        {
            try
            {
                _numRepititions = Integer.parseInt(repititions.Content);
            }
            catch(NumberFormatException e)
            {
                this.errorDetected.execute("Unknown number of repetitions '" + repititions.Content + "'.", repititions.getLine(), repititions.getColumn());
                noError = false;
            }
        }

        ParamExtractor sParam = new ParamExtractor("s", this.params, this.errorDetected);

        if(sParam.ContainsSwitch)
        {
            if(sParam.PostSwitchParam != null)
            {
                try
                {
                    _threshold = Double.parseDouble(sParam.PostSwitchValue);
                }
                 catch(NumberFormatException e)
                {
                    this.errorDetected.execute("Unknown threshold value '" + repititions.Content + "'.", sParam.PostSwitchParam.getLine(), sParam.PostSwitchParam.getColumn());
                    noError = false;
                }
            }
            else
            {
                 this.errorDetected.execute("Expected value after switch 's'.", sParam.SwitchParam.getLine(), sParam.SwitchParam.getColumn());
                 noError = false;
            }
        }

        final int stCommandIndex =  sParam.ContainsSwitch ? 2 : 1;
        final ParameterIdent stCommand = this.assertParameterIdent(stCommandIndex);
        noError = noError && stCommand != null;

        if(stCommand != null)
        {
            _exampleCommand = new SyntaxCommand() {
                public int getLine() {
                    return stCommand.getLine();  //To change body of implemented methods use File | Settings | File Templates.
                }

                public int getColumn() {
                    return stCommand.getColumn();  //To change body of implemented methods use File | Settings | File Templates.
                }

                public String getName() {
                    return stCommand.Content;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Iterable<Parameter> getParameters() {
                    return InferST_Bootstrap.this.params.subList(stCommandIndex + 1, InferST_Bootstrap.this.params.size());
                }
            };

        InferSTBase inferCommand = (InferSTBase) CommandFactory.make(_exampleCommand, this.sourceIdentToNetwork, this.errorDetected, _rand);
        noError = noError && inferCommand.checkParams();
        }


        return noError;

    }

    @Override
      protected String produceResult() {

        for(int i=0; i<_numRepititions; i++)
        {
            final ParameterIdentSet originalSet = (ParameterIdentSet)_exampleCommand.getParameters().iterator().next();
            final Iterable<String> resampledGeneTreeIdents = resampleGeneTrees(originalSet);
            SyntaxCommand executionCommand = new SyntaxCommand() {
                public int getLine() {
                    return _exampleCommand.getLine();
                }

                public int getColumn() {
                   return  _exampleCommand.getColumn();  //To change body of implemented methods use File | Settings | File Templates.
                }

                public String getName() {
                   return _exampleCommand.getName();
                }

                public Iterable<Parameter> getParameters() {

                    Object[] params = IterableHelp.toArray(_exampleCommand.getParameters());
                    params[0] = new ParameterIdentSet(originalSet.getLine(), originalSet.getColumn(), resampledGeneTreeIdents);

                    ArrayList<Parameter> tbr = new ArrayList<Parameter>();

                    for(Object obj : params)
                    {
                        tbr.add((Parameter)obj);
                    }

                    return tbr;
                }
            };

            InferSTBase inferCommand = (InferSTBase) CommandFactory.make(executionCommand, this.sourceIdentToNetwork, this.errorDetected, _rand);
            inferCommand.checkParams();
            inferCommand.produceResult();
        }

        return "";
      }

    private Iterable<String> resampleGeneTrees(ParameterIdentSet parameter)
    {
        Object[] geneTreeIdents = IterableHelp.toArray(parameter.Elements);

        ArrayList<String> sample = new ArrayList<String>();

        for(int i = 0; i<geneTreeIdents.length; i++)
        {
            sample.add((String)geneTreeIdents[_rand.nextInt(geneTreeIdents.length)]);
        }

        return sample;
    }

}
