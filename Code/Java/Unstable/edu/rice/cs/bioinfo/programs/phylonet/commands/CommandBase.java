package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import sun.text.normalizer.Trie;

import java.io.IOException;
import java.util.ArrayList;
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

    private SyntaxCommand _motivatingCommand;

    protected final ArrayList<Parameter> params;

    protected  final Map<String,NetworkNonEmpty> sourceIdentToNetwork;

    protected  final Proc3<String, Integer, Integer> errorDetected;

    private boolean _checkParamsCalled = false;

    CommandBase(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String,NetworkNonEmpty> sourceIdentToNetwork,
                Proc3<String, Integer, Integer> errorDetected)
    {
        _motivatingCommand = motivatingCommand;
        this.params = params;
        this.sourceIdentToNetwork = sourceIdentToNetwork;
        this.errorDetected = errorDetected;
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

        if(ident != null)
        {
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

    protected LinkedList<NetworkNonEmpty> assertNetworksExist(ParameterIdentSet identList)
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

    public void executeCommand(Proc<String> displayResult) throws IOException
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

    protected abstract void executeCommandHelp(Proc<String> displayResult) throws IOException;

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

     protected ParameterQuote assertQuotedParameter(int paramIndent)
    {
        Parameter param = this.params.get(paramIndent);

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
                errorDetected.execute("Expected identifier set, found identifier.", parameterIdent.getColumn(), parameterIdent.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public ParameterIdentSet forIdentList(ParameterIdentList parameterIdentList, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier set, found identifier list.", parameterIdentList.getColumn(), parameterIdentList.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public ParameterIdentSet forQuote(ParameterQuote parameterQuote, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier set, found quoted identifier.", parameterQuote.getColumn(), parameterQuote.getLine());
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public ParameterIdentSet forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Object o) throws RuntimeException {
                errorDetected.execute("Expected identifier set, found identifier set list.", parameterTaxonSetList.getColumn(), parameterTaxonSetList.getLine());
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


}
