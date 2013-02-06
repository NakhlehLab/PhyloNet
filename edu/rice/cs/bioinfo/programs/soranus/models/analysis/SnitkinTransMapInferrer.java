package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012.SnitkinTransMapInferrerTemplateDoubleDistance;
import org.joda.time.LocalDate;

import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/3/13
 * Time: 5:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrer extends SnitkinTransMapInferrerTemplateDoubleDistance<Integer,Sequencing>
{
    public SnitkinTransMapInferrer(Map<Integer, Map<LocalDate, Object>> patientTraces, Map<Integer, LocalDate> patientToFirstPositiveDate, Map<Sequencing, Integer> sequencingToPatient) {
        super(patientTraces, patientToFirstPositiveDate, sequencingToPatient);
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
