package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012.SnitkinTransMapInferrerTemplateDoubleDistance;
import edu.rice.cs.bioinfo.programs.soranus.models.data.Sequencing;
import org.joda.time.LocalDate;

import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/3/13
 * Time: 5:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrer<E> extends SnitkinTransMapInferrerTemplateDoubleDistance<E,Sequencing>
{
    public SnitkinTransMapInferrer(Map<E, Map<LocalDate, Object>> patientTraces,
                                   Map<E, LocalDate> firstPositiveDate, Map<Sequencing, E> sequencings) {
        super(patientTraces, firstPositiveDate, sequencings);
    }

    @Override
    protected Double getMaxDistance() {
        return 10000.0;
    }

    @Override
    protected int getGeneticDistance(Sequencing sequence1, Sequencing sequence2) {
        return sequence1.getGeneticDistance(sequence2);
    }
}
