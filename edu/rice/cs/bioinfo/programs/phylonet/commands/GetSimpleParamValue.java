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
