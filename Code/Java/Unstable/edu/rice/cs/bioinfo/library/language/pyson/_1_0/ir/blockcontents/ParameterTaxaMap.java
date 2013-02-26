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

package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.Identifier;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.TaxaMap;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.TaxaMapEntry;

import java.util.AbstractMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 2:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParameterTaxaMap extends ParameterBase
{
    public final LinkedList<Map.Entry<String, List<String>>> _mappings = new LinkedList<Map.Entry<String, List<String>>>();

    public ParameterTaxaMap(int line, int column, TaxaMap map) {
        super(line, column);

       for(TaxaMapEntry e : map.Entries)
       {
           String key = e.Key.Content;

           List<String> values = new LinkedList<String>();
           for(Identifier value : e.Values)
           {
                values.add(value.Content);
           }
           _mappings.add(new AbstractMap.SimpleEntry<String, List<String>>(key, values));
       }

    }

    public <R, T, E extends Exception> R execute(ParameterAlgo<R, T, E> algo, T input) throws E {

        return algo.forTaxaMap(this, input);

    }
}
