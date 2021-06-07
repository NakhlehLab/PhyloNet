package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.commands.GetSimpleParamValue;
import edu.rice.cs.bioinfo.programs.phylonet.commands.ParamExtractor;
import edu.rice.cs.bioinfo.programs.phylonet.commands.ParamExtractorAllelMap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ParamParser
{
    List<ParamExtractor> switches;
    List<Parameter> params;
    List<String> switchParamsNames;
    Proc3<String, Integer, Integer> errorDetected;

    public ParamParser(List<Parameter> params,  Proc3<String, Integer, Integer> errorDetected)
    {
        this.params = params;
        this.switches = new ArrayList<ParamExtractor>();
        this.switchParamsNames = new ArrayList<String>();
        this.errorDetected = errorDetected;

    }

    private String getParamName(Parameter p)
    {
        return p.execute(GetSimpleParamValue.Singleton,null);
    }

    private ParamExtractor createParamExtractor(String paramName)
    {
        ParamExtractor result = new ParamExtractor(paramName, this.params, this.errorDetected);
        switchParamsNames.add(paramName);
        switches.add(result);
        return result;
    }

    private ParamExtractorAllelMap createAlleleMapParamExtractor(String paramName)
    {
        ParamExtractorAllelMap result =  new ParamExtractorAllelMap(paramName, this.params, this.errorDetected);
        switchParamsNames.add(paramName);
        switches.add(result);
        return result;
    }


    public boolean getBooleanParameter(String paramName)
    {
        ParamExtractor paramExtractor = createParamExtractor(paramName);
        return paramExtractor.ContainsSwitch;
    }

    public String getStringParameter(String paramName)
    {
        ParamExtractor paramExtractor = createParamExtractor(paramName);
        if (!paramExtractor.ContainsSwitch)
            return null;

        if (paramExtractor.PostSwitchValue == null)
            throw new RuntimeException("Parameter -" + paramName + " needs a value.");

        return paramExtractor.PostSwitchValue;
    }

    public Integer getIntegerParameter(String paramName)
    {
        String value =  getStringParameter(paramName);
        if (value == null)
            return null;
        else
        {
            try
            {
                return Integer.parseInt(value);
            }
            catch (NumberFormatException e)
            {
                throw new RuntimeException(String.format("-%s requires an integer value. Instead, it was passed \"%s\".",paramName,value));
            }
        }
    }

    public Double getDoubleParameter(String paramName)
    {
        String value = getStringParameter(paramName);
        if (value == null)
            return null;
        else
        {
            try
            {
                return Double.parseDouble(value);
            }
            catch (NumberFormatException e)
            {
                throw new RuntimeException(String.format("-%s requires a double value. Instead, it was passed \"%s\".",paramName,value));
            }
        }
    }



    public Map<String,List<String>> getSpeciesToAlleleMapParameter(String paramName)
    {
        ParamExtractorAllelMap paramExtractor = createAlleleMapParamExtractor(paramName);
        if (!paramExtractor.ContainsSwitch)
            return null;

        if (!paramExtractor.IsValidMap)
            throw new RuntimeException("Allele map parameter " + paramName + " does not have a valid map.");

        return convertParamAlleleMapToActual(paramExtractor.ValueMap);

    }

    private Map<String, List<String>> convertParamAlleleMapToActual(Map<String,String> ValueMap)
    {
        Map<String, List<String>> result = new HashMap<String, List<String>>();

        for (String allele : ValueMap.keySet())
        {
            String species = ValueMap.get(allele);
            if (!result.containsKey(species))
                result.put(species,new ArrayList<String>());

            result.get(species).add(allele);
        }

        return result;
    }


    public void checkNoInvalidParameters()
    {
        for (Parameter param : params)
        {
            String paramName = getParamName(param);
            if (paramName != null && paramName.startsWith("-"))
            {
                String paramSwitchName = paramName.substring(1);
                if (!switchParamsNames.contains(paramSwitchName))
                {
                    throw new RuntimeException(String.format("%s is not a valid parameter.",paramSwitchName));
                }
            }
        }
    }

    public ParamExtractor[] getSwitches()
    {
        return switches.toArray(new ParamExtractor[switches.size()]);
    }
}
