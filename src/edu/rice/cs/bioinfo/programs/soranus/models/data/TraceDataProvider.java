package edu.rice.cs.bioinfo.programs.soranus.models.data;

import org.joda.time.LocalDate;

import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/1/13
 * Time: 4:09 PM
 * To change this template use File | Settings | File Templates.
 */
public interface TraceDataProvider<R,S,E extends Throwable>
{
    Map<S,Map<LocalDate,Object>> getTraces(R traceDataRecord) throws E;
}
