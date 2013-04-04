package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.programs.soranus.models.data.Sequencing;

import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/5/13
 * Time: 1:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrerGeneticOnly extends edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012.SnitkinTransMapInferrerGeneticOnly<Integer,Sequencing,Integer>
{
    public SnitkinTransMapInferrerGeneticOnly(Map<Sequencing, Integer> sequencingToPatient, Func2<Sequencing, Sequencing, Integer> getGeneticDistance, Func2<Integer, Integer, Integer> add,
                                              Func2<Integer, Integer, Integer> subtract, Func1<Integer, Integer> makeDistance, Func<Integer> getMaxDistance) {
        super(sequencingToPatient, getGeneticDistance, add, subtract, makeDistance, getMaxDistance);
    }
}
