package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/13/11
 * Time: 11:19 AM
 * To change this template use File | Settings | File Templates.
 */
public class GetSimpleParamValue implements ParameterAlgo<String, Object, RuntimeException>
{
    public static final GetSimpleParamValue Singleton = new GetSimpleParamValue();

    private GetSimpleParamValue()
    {

    }

    public String forIdentifier(ParameterIdent parameterIdent, Object o) throws RuntimeException {
        return parameterIdent.Content;
    }

    public String forIdentList(ParameterIdentList parameterIdentList, Object o) throws RuntimeException {
        return null;
    }

    public String forQuote(ParameterQuote parameterQuote, Object o) throws RuntimeException {
        return parameterQuote.UnquotedText;
    }

    public String forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Object o) throws RuntimeException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String forIdentSet(ParameterIdentSet parameterIdentSet, Object o) throws RuntimeException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String forTaxaMap(ParameterTaxaMap parameterTaxaMap, Object o) throws RuntimeException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
