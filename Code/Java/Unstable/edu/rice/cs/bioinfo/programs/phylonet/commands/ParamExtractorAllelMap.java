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
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterTaxaMap;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/22/11
 * Time: 1:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParamExtractorAllelMap extends ParamExtractor
{
    public final HashMap<String,String> ValueMap = new HashMap<String, String>();

    public final boolean IsValidMap;

    public ParamExtractorAllelMap(String switchValue, List<Parameter> params, Proc3<String, Integer, Integer> errorDetected)
    {
        super(switchValue, params, errorDetected);

        boolean isValidMap = false;
        if(this.ContainsSwitch)
        {
            isValidMap = true;
            ParameterTaxaMap map = (ParameterTaxaMap)this.PostSwitchParam;

            for(Map.Entry<String, List<String>> entry : map._mappings)
            {
                for(String allele : entry.getValue())
                {
                    if(!ValueMap.containsKey(allele))
                    {
                        ValueMap.put(allele, entry.getKey());
                    }
                    else
                    {
                        errorDetected.execute("Duplicate allele '" + allele  +"'", map.getLine(), map.getColumn());
                        isValidMap = false;
                    }
                }
            }
        }

        IsValidMap = isValidMap;
    }
}
