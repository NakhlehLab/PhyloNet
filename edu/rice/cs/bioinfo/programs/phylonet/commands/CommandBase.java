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
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.ContainsHybridNode;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/19/11
 * Time: 6:35 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class CommandBase implements Command {

    protected SyntaxCommand _motivatingCommand;

    protected final ArrayList<Parameter> params;

    protected final Map<String,NetworkNonEmpty> sourceIdentToNetwork;

    protected final Proc3<String, Integer, Integer> errorDetected;

    private boolean _checkParamsCalled = false;

    protected final RichNewickReader<Networks> rnReader;

    private ArrayList<Proc1<String>> _rnGeneratedObservers = new ArrayList<Proc1<String>>();

    public void addRichNewickGeneratedListener(Proc1<String> listener)
    {
        _rnGeneratedObservers.add(listener);
    }

    protected void richNewickGenerated(String richNewick)
    {
        if(_motivatingCommand.getAssigment() != null)
        {
            Networks networks = rnReader.readAnyErrorToRuntimeException(richNewick);
            NetworkNonEmpty networkNonEmpty = networks.Networks.iterator().next();
            sourceIdentToNetwork.put(_motivatingCommand.getAssigment().Identifier, networkNonEmpty);
        }

        for(Proc1<String> observer : _rnGeneratedObservers)
        {
            observer.execute(richNewick);
        }
    }

    public CommandBase(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String,NetworkNonEmpty> sourceIdentToNetwork,
                Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader)
    {
        _motivatingCommand = motivatingCommand;
        this.params = params;
        this.sourceIdentToNetwork = sourceIdentToNetwork;
        this.errorDetected = errorDetected;
        this.rnReader = rnReader;
    }

    private boolean assertParamsCount(int minParamCount, int maxParamCount)
    {
        boolean noError = true;
        if(minParamCount > params.size())
        {
               SyntaxCommand directive = this.getDefiningSyntaxCommand();
                errorDetected.execute(String.format("Expected at least %s parameters for command %s but found %s.",
                                                    minParamCount, directive.getName(), params.size()),
                                                    directive.getLine(), directive.getColumn());
               noError = false;
        }

         if(maxParamCount < params.size())
         {
                SyntaxCommand directive = this.getDefiningSyntaxCommand();
                errorDetected.execute(String.format("Expected at most %s parameters for command %s but found %s.",
                                                    maxParamCount, directive.getName(), params.size()),
                                                    directive.getLine(), directive.getColumn());
               noError = false;
         }

        return noError;
    }

    protected NetworkNonEmpty assertAndGetNetwork(int paramIndex)
    {
        ParameterIdent ident = this.assertParameterIdent(paramIndex);

        if(ident == null)
            return null;

        if(!sourceIdentToNetwork.containsKey(ident.Content))
        {
            errorDetected.execute(String.format("Unknown identifier '%s'.", ident.Content), ident.getLine(), ident.getColumn());
            return null;
        }
        else
        {
            return sourceIdentToNetwork.get(ident.Content);
        }


    }

    protected NetworkNonEmpty assertAndGetTree(int paramIndex)
    {
        NetworkNonEmpty network = assertAndGetNetwork(paramIndex);

        if(network != null)
        {
            if(network.execute(ContainsHybridNode.Singleton, null))
            {
                ParameterIdent p = (ParameterIdent) params.get(paramIndex);
                errorDetected.execute(String.format("Expected '%s' to be a tree but found network with hybrid node.", p.Content), p.getLine(), p.getColumn());
                return null;
            }
            else
            {
                return network;
            }
        }
        else
        {
            return null;
        }

    }



    private boolean assertNetworkExists(ParameterIdent parameterIdent)
    {
        return assertNetworkExists(parameterIdent.Content, parameterIdent.getLine(), parameterIdent.getColumn());
    }

    protected boolean assertNetworkExists(String networkName, int line, int col)
    {
        if(!sourceIdentToNetwork.containsKey(networkName))
        {
            errorDetected.execute(String.format("Unknown identifier '%s'.", networkName), line, col);
            return false;
        }
        else
        {
            return true;
        }

    }

    protected boolean assertTreeExists(String treeName, int line, int col)
    {
        if(!sourceIdentToNetwork.containsKey(treeName))
        {
            errorDetected.execute(String.format("Unknown identifier '%s'.", treeName), line, col);
            return false;
        }
        else
        {
            Network network = sourceIdentToNetwork.get(treeName);
            if(network.execute(ContainsHybridNode.Singleton, null))
            {
                errorDetected.execute(String.format("Expected '%s' to be a tree but found network with hybrid node.", treeName), line,col);
                return false;
            }

            return true;
        }

    }



    protected LinkedList<NetworkNonEmpty> assertNetworksExist(ParameterIdentSet identSet)
    {
        LinkedList<NetworkNonEmpty> accum = new LinkedList<NetworkNonEmpty>();

        for(String ident : identSet.Elements)
        {
            if(assertNetworkExists(ident, identSet.getLine(), identSet.getColumn()))
            {
                if(accum != null)
                {
                    accum.add(sourceIdentToNetwork.get(ident));
                }
            }
            else
            {
                accum = null;
            }
        }

        return accum;
    }

     protected LinkedList<NetworkNonEmpty> assertTreesExist(ParameterIdentSet identSet)
    {
        LinkedList<NetworkNonEmpty> accum = new LinkedList<NetworkNonEmpty>();

        for(String ident : identSet.Elements)
        {
            if(assertTreeExists(ident, identSet.getLine(), identSet.getColumn()))
            {
                if(accum != null)
                {
                    accum.add(sourceIdentToNetwork.get(ident));
                }
            }
            else
            {
                accum = null;
            }
        }

        return accum;
    }

    protected LinkedList<NetworkNonEmpty> assertNetworksExist(ParameterIdentList identList)
    {
        LinkedList<NetworkNonEmpty> accum = new LinkedList<NetworkNonEmpty>();

        for(String ident : identList.Elements)
        {
            if(assertNetworkExists(ident, identList.getLine(), identList.getColumn()))
            {
                if(accum != null)
                {
                    accum.add(sourceIdentToNetwork.get(ident));
                }
            }
            else
            {
                accum = null;
            }
        }

        return accum;
    }

    public boolean checkParams()
    {
        _checkParamsCalled = true;

        boolean noError = assertParamsCount(getMinNumParams(), getMaxNumParams());

        if(noError)
        {
            return checkParamsForCommand();
        }
        else
        {
            return noError;
        }
    }

    protected abstract int getMinNumParams();

    protected abstract int getMaxNumParams();

    protected abstract boolean checkParamsForCommand();

    public SyntaxCommand getDefiningSyntaxCommand()
    {
        return _motivatingCommand;
    }

    public void executeCommand(Proc1<String> displayResult) throws IOException
    {
        if(!_checkParamsCalled)
        {
            if(checkParams())
            {
                executeCommandHelp(displayResult);
            }
        }
        else
        {
             executeCommandHelp(displayResult);
        }
    }

    protected abstract void executeCommandHelp(Proc1<String> displayResult) throws IOException;

    protected ParameterIdent assertParameterIdent(int paramIndent)
    {
        Parameter param = this.params.get(paramIndent);

        return param.execute(new ParameterAlgo<ParameterIdent, Object, RuntimeException>() {

            public ParameterIdent forIdentifier(ParameterIdent parameterIdent, Object o) throws RuntimeException {
                return parameterIdent;
            }

            public ParameterIdent forIdentList(ParameterIdentList parameterIdentList, Object o) throws RuntimeException {

                errorDetected.execute("Expected identifier, found parameter identifier list.", parameterIdentList.getLine(), parameterIdentList.getColumn());
                return null;
            }

            public ParameterIdent forQuote(ParameterQuote parameterQuote, Object o) throws RuntimeException {

                errorDetected.execute("Expected identifier, found quoted parameter.", parameterQuote.getLine(), parameterQuote.getColumn());
                return null;
            }

            public ParameterIdent forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Object o) throws RuntimeException {
               errorDetected.execute("Expected identifier, found taxon set list.", parameterTaxonSetList.getLine(), parameterTaxonSetList.getColumn());
                return null;
            }

            public ParameterIdent forIdentSet(ParameterIdentSet parameterIdentSet, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier, found identifier set list.", parameterIdentSet.getLine(), parameterIdentSet.getColumn());
                return null;
            }

            public ParameterIdent forTaxaMap(ParameterTaxaMap parameterTaxaMap, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier, found taxa map.", parameterTaxaMap.getLine(), parameterTaxaMap.getColumn());
                return null;
            }
        }, null);
    }

     protected ParameterQuote assertQuotedParameter(int paramIndex)
    {
        Parameter param = this.params.get(paramIndex);

        return param.execute(new ParameterAlgo<ParameterQuote, Object, RuntimeException>() {

            public ParameterQuote forIdentifier(ParameterIdent parameterIdent, Object o) throws RuntimeException {
                errorDetected.execute("Expected quoted parameter, found identifier.", parameterIdent.getLine(), parameterIdent.getColumn());
                return null;
            }

            public ParameterQuote forIdentList(ParameterIdentList parameterIdentList, Object o) throws RuntimeException {

                errorDetected.execute("Expected quoted parameter, found parameter identifier list.", parameterIdentList.getLine(), parameterIdentList.getColumn());
                return null;
            }

            public ParameterQuote forQuote(ParameterQuote parameterQuote, Object o) throws RuntimeException {

               return parameterQuote;
            }

            public ParameterQuote forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Object o) throws RuntimeException {
               errorDetected.execute("Expected quoted parameter, found taxon set list.", parameterTaxonSetList.getLine(), parameterTaxonSetList.getColumn());
                return null;
            }

            public ParameterQuote forIdentSet(ParameterIdentSet parameterIdentSet, Object o) throws RuntimeException {
                errorDetected.execute("Expected quoted parameter, found identifier set list.", parameterIdentSet.getLine(), parameterIdentSet.getColumn());
                return null;
            }

            public ParameterQuote forTaxaMap(ParameterTaxaMap parameterTaxaMap, Object o) throws RuntimeException {
                errorDetected.execute("Expected quoted parameter, found taxa map.", parameterTaxaMap.getLine(), parameterTaxaMap.getColumn());
                return null;
            }
        }, null);
    }

    protected ParameterIdentSet assertParameterIdentSet(int paramIndex)
    {
        Parameter param = this.params.get(paramIndex);
        return param.execute(new ParameterAlgo<ParameterIdentSet, Object, RuntimeException>()
        {
            public ParameterIdentSet forIdentifier(ParameterIdent parameterIdent, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier set, found identifier.", parameterIdent.getLine(), parameterIdent.getColumn());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public ParameterIdentSet forIdentList(ParameterIdentList parameterIdentList, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier set, found identifier list.", parameterIdentList.getLine(), parameterIdentList.getColumn());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public ParameterIdentSet forQuote(ParameterQuote parameterQuote, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier set, found quoted parameter.", parameterQuote.getLine(), parameterQuote.getColumn());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public ParameterIdentSet forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier set, found identifier set list.", parameterTaxonSetList.getLine(), parameterTaxonSetList.getColumn());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public ParameterIdentSet forIdentSet(ParameterIdentSet parameterIdentSet, Object o) throws RuntimeException {

                return parameterIdentSet;

            }

            public ParameterIdentSet forTaxaMap(ParameterTaxaMap parameterTaxaMap, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier set, found taxa map.", parameterTaxaMap.getColumn(), parameterTaxaMap.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }
        }, null);
    }


    protected ParameterIdentList assertParameterIdentList(int paramIndex)
    {
        Parameter param = this.params.get(paramIndex);
        return param.execute(new ParameterAlgo<ParameterIdentList, Object, RuntimeException>()
        {
            public ParameterIdentList forIdentifier(ParameterIdent parameterIdent, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier list, found identifier.", parameterIdent.getColumn(), parameterIdent.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public ParameterIdentList forIdentList(ParameterIdentList parameterIdentList, Object o) throws RuntimeException {
                  return parameterIdentList;
            }

            public ParameterIdentList forQuote(ParameterQuote parameterQuote, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier list, found quoted parameter.", parameterQuote.getColumn(), parameterQuote.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public ParameterIdentList forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier list, found identifier set list.", parameterTaxonSetList.getColumn(), parameterTaxonSetList.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public ParameterIdentList forIdentSet(ParameterIdentSet parameterIdentSet, Object o) throws RuntimeException {

                  errorDetected.execute("Expected identifier list, found identifier set.", parameterIdentSet.getColumn(), parameterIdentSet.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.

            }

            public ParameterIdentList forTaxaMap(ParameterTaxaMap parameterTaxaMap, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier list, found taxa map.", parameterTaxaMap.getColumn(), parameterTaxaMap.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }
        }, null);
    }



    protected Parameter assertParameterIdentListOrSetList(int paramIndex)
    {
        Parameter param = this.params.get(paramIndex);
        return param.execute(new ParameterAlgo<Parameter, Object, RuntimeException>()
        {
            public Parameter forIdentifier(ParameterIdent parameterIdent, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier list or set list, found identifier.", parameterIdent.getColumn(), parameterIdent.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public ParameterIdentList forIdentList(ParameterIdentList parameterIdentList, Object o) throws RuntimeException {
                return parameterIdentList;
            }

            public Parameter forQuote(ParameterQuote parameterQuote, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier list or set list, found quoted parameter.", parameterQuote.getColumn(), parameterQuote.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public Parameter forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Object o) throws RuntimeException {
                return parameterTaxonSetList;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public Parameter forIdentSet(ParameterIdentSet parameterIdentSet, Object o) throws RuntimeException {

                errorDetected.execute("Expected identifier list or set list, found identifier set.", parameterIdentSet.getColumn(), parameterIdentSet.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.

            }

            public ParameterTaxaMap forTaxaMap(ParameterTaxaMap parameterTaxaMap, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier list or set list, found identifier taxa map.", parameterTaxaMap.getColumn(), parameterTaxaMap.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }
        }, null);
    }





    protected  boolean checkForUnknownSwitches(String... expectedSwitches)
    {
        boolean  noError = true;

        HashSet<String> exepctedLower = new HashSet<String>();

        for(String s : expectedSwitches)
        {
            exepctedLower.add(s.toLowerCase());
        }

        for(Parameter p : this.params)
        {
            if(p instanceof ParameterIdent)
            {
                ParameterIdent ident = (ParameterIdent)p;
                if(ident.Content.startsWith("-"))
                {
                    if(ident.Content.length() > 1)
                    {
                        String switchKey = ident.Content.substring(1).toLowerCase();

                        if(!exepctedLower.contains(switchKey))
                        {
                            errorDetected.execute("Unknown switch -" + switchKey, p.getLine(), p.getColumn());
                            noError = false;
                        }
                    }
                    else
                    {
                        errorDetected.execute("Expected switch name after -.", p.getLine(), p.getColumn());
                        noError = false;
                    }
                }
            }
        }

        return noError;
    }


}
