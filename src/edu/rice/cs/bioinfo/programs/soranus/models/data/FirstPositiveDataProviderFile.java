package edu.rice.cs.bioinfo.programs.soranus.models.data;

import edu.rice.cs.bioinfo.programs.soranus.models.factories.FirstPositiveMapFactory;
import org.joda.time.LocalDate;

import java.io.File;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/1/13
 * Time: 4:21 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class FirstPositiveDataProviderFile<S>
        implements FirstPositiveDataProvider<File,S,Exception>
{
    public Map<S, LocalDate> getFirstPositiveDates(File dataRecord) throws Exception
    {
        return new FirstPositiveMapFactory<S,LocalDate>()
        {
            @Override
            protected LocalDate makeDate(String dateText) {
                return new LocalDate(dateText);
            }

            @Override
            protected S makeId(String sourceIdText) {
                return makeSourceId(sourceIdText);
            }

        }.makeMapFromXMLFile(dataRecord);
    }

    protected abstract S makeSourceId(String sourceIdText);
}
