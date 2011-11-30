package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterTaxaMap;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/22/11
 * Time: 3:00 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParamExtractorAllelListMap extends ParamExtractor
{
    public final HashMap<String,List<String>> ValueMap = new HashMap<String, List<String>>();

    public final boolean IsValidMap;

    ParamExtractorAllelListMap(String switchValue, List<Parameter> params, Proc3<String, Integer, Integer> errorDetected)
    {
        super(switchValue, params, errorDetected);
        HashSet<String> allAlleles = new HashSet<String>();

        boolean isValidMap = false;
        if(this.ContainsSwitch)
        {
            isValidMap = true;
            ParameterTaxaMap map = (ParameterTaxaMap)this.PostSwitchParam;

            for(Map.Entry<String, List<String>> entry : map._mappings)
            {
                for(String allele : entry.getValue())
                {
                    if(allAlleles.contains(allele))
                    {
                        errorDetected.execute("Duplicate allele '" + allele  +"'", map.getLine(), map.getColumn());
                        isValidMap = false;
                    }
                }
                ValueMap.put(entry.getKey(), entry.getValue());
            }
        }

        IsValidMap = isValidMap;
    }
}
