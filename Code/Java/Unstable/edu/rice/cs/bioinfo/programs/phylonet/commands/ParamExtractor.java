package edu.rice.cs.bioinfo.programs.phylonet.commands;

import com.sun.org.apache.bcel.internal.generic.NEW;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;

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

    ParamExtractor(String switchValue, List<Parameter> params)
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
                return;
            }
        }
        ContainsSwitch = false;
        SwitchParam = null;
        PostSwitchValue = null;
        PostSwitchParam = null;

    }
}


