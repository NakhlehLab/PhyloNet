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
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.consensus.TreeConsensusCalculator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.StringReader;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/29/11
 * Time: 11:08 AM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("Infer_ST_BOOTSTRAP")
public class InferST_Bootstrap extends CommandBaseFileOut {

    private int _numRepititions;

    private double _threshold = 0.5;

    private Random _rand;

    private  SyntaxCommand _exampleCommand;

    public InferST_Bootstrap(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected,
                             RichNewickReader<Networks> rnReader, Random rand) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
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

        final int stCommandIndex =  sParam.ContainsSwitch ? 3 : 1;
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

                public AssignmentIdent getAssigment() {
                    return getDefiningSyntaxCommand().getAssigment();  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Iterable<Parameter> getParameters() {
                    return InferST_Bootstrap.this.params.subList(stCommandIndex + 1, InferST_Bootstrap.this.params.size());
                }
            };

        InferSTBase inferCommand = (InferSTBase) CommandFactory.make(_exampleCommand, this.sourceIdentToNetwork, this.errorDetected, rnReader, _rand);
        noError = noError && inferCommand.checkParams();
        }


        return noError;

    }

    @Override
      protected String produceResult() {

        Set<Tree> st_set = new LinkedHashSet<Tree>();
        TreeConsensusCalculator tcc = new TreeConsensusCalculator();
        for(int i=0; i<_numRepititions; i++)
        {
            final ParameterIdentList originalSet = (ParameterIdentList)_exampleCommand.getParameters().iterator().next();
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

                public AssignmentIdent getAssigment() {
                    return _exampleCommand.getAssigment();
                }

                public Iterable<Parameter> getParameters() {

                    Object[] params = IterableHelp.toArray(_exampleCommand.getParameters());
                    params[0] = new ParameterIdentList(originalSet.getLine(), originalSet.getColumn(), resampledGeneTreeIdents);

                    ArrayList<Parameter> tbr = new ArrayList<Parameter>();

                    for(Object obj : params)
                    {
                        tbr.add((Parameter)obj);
                    }

                    return tbr;
                }
            };

            InferSTBase inferCommand = (InferSTBase) CommandFactory.make(executionCommand, this.sourceIdentToNetwork, this.errorDetected, this.rnReader, _rand);

            final Set<Tree> sts = new LinkedHashSet<Tree>();
            inferCommand.addRichNewickGeneratedListener(new Proc1<String>() {
                public void execute(String tree) {
                    NewickReader nr = new NewickReader(new StringReader(tree));

                    try {
                        sts.add(nr.readTree());
                    } catch (Exception e) {
                        throw new RuntimeException(e);
                    }

                }
            });
            boolean noError = inferCommand.checkParams();
            if(!noError)
            {
                throw new RuntimeException("Failed execution error check: " + this.getDefiningSyntaxCommand().getName());
            }
            inferCommand.produceResult();

            if(sts.size() == 1){
					st_set.add(sts.iterator().next());
				}
				else{
					st_set.add(tcc.computeRootedConsensus(sts, 0.5));
				}
        }

        tcc.setOutputEdgeWeightsAreSupports(true);
		Tree st = tcc.computeRootedConsensus(st_set, _threshold);

        String stString = st.toNewick();

        this.richNewickGenerated(stString);
        return "\n" + stString;
      }

    private Iterable<String> resampleGeneTrees(ParameterIdentList parameter)
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
