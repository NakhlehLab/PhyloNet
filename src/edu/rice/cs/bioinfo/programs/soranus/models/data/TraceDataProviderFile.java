package edu.rice.cs.bioinfo.programs.soranus.models.data;

import edu.rice.cs.bioinfo.programs.soranus.models.factories.PatientTraceMapFactory;
import org.joda.time.LocalDate;

import java.io.File;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/1/13
 * Time: 4:27 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TraceDataProviderFile<S> implements TraceDataProvider<File,S,Exception>
{
    public Map<S, Map<LocalDate, Object>> getTraces(File traceDataRecord) throws Exception {

        return new PatientTraceMapFactory<S,LocalDate,Object>()
        {
            @Override
            protected Object makeLocation(String locationText) {
                return locationText;  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected LocalDate makeDate(String dateText) {
                return new LocalDate(dateText);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected S makeId(String sourceIdText) {
                return makeSourceId(sourceIdText);  //To change body of implemented methods use File | Settings | File Templates.
            }
        }.makeMapFromXMLFile(traceDataRecord);

    }

    protected abstract S makeSourceId(String sourceIdText);
}
