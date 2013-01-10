package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import edu.rice.cs.bioinfo.library.programming.Func2;
import org.joda.time.LocalDate;

import java.math.BigDecimal;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/20/12
 * Time: 11:55 AM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrerDoubleDistanceTest extends SnitkinTransMapInferrerTestBase<Double> {


    @Override
    protected SnitkinTransMapInferrer<Integer, Genome, Double> makeInferrer(Integer rootPatient, List<Genome> genomes, HashMap<Genome, Integer> genomeToPatient, Func2<Genome, Genome, Double> genDistance, Map<Integer, Map<LocalDate, Object>> patientTraceData, Map<Integer, LocalDate> patientToFirstPositive) {
         return new SnitkinTransMapInferrerDoubleDistance<Integer,Genome>(rootPatient, genomes, genDistance, patientTraceData, patientToFirstPositive, genomeToPatient);

    }
}
