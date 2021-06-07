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
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/31/11
 * Time: 4:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParamExtractor
{
    public final boolean ContainsSwitch;

    public final Parameter SwitchParam;

    public final String PostSwitchValue;

    public final Parameter PostSwitchParam;

    public final boolean DuplicateSwitch;

    public ParamExtractor(String switchValue, List<Parameter> params, Proc3<String, Integer, Integer> errorDetected)
    {
        for(int i = 0; i<params.size(); i++)
        {
            Parameter ithParam = params.get(i);
            String simpleValue = ithParam.execute(GetSimpleParamValue.Singleton, null);

            if(simpleValue != null && simpleValue.toLowerCase().equals("-" + switchValue))
            {

                ContainsSwitch = true;
                SwitchParam = ithParam;
                if(i<params.size()-1)
                {
                    PostSwitchParam = params.get(i + 1);
                    PostSwitchValue = PostSwitchParam.execute(GetSimpleParamValue.Singleton, null);
                }
                else
                {
                    PostSwitchValue = null;
                    PostSwitchParam = null;

                }

                for(int j = i + 1; j<params.size(); j++)
                {
                    Parameter jthParam = params.get(j);
                    simpleValue = jthParam.execute(GetSimpleParamValue.Singleton, null);

                     if(simpleValue != null && simpleValue.toLowerCase().equals("-" + switchValue))
                     {
                         errorDetected.execute("Duplicate switch -" + switchValue, jthParam.getLine(), jthParam.getColumn());
                         DuplicateSwitch = true;
                         return;
                     }
                }
                DuplicateSwitch = false;

                return;
            }
        }
        ContainsSwitch = false;
        SwitchParam = null;
        PostSwitchValue = null;
        PostSwitchParam = null;
        DuplicateSwitch = false;

    }
}


