package edu.rice.cs.bioinfo.programs.soranus.models.data;

import edu.rice.cs.bioinfo.programs.soranus.models.factories.SequencingsMapFactory;

import java.io.File;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/1/13
 * Time: 3:42 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SequencingsDataProviderFile<S> implements SequencingsDataProvider<File,S,Exception>
{
    public Map<Sequencing, S> getSequencingToPatientMap(File sequencingsXml) throws Exception {

        return new SequencingsMapFactory<S>()
        {
            @Override
            protected S makeId(String sourceIdText) {
                return makeSourceId(sourceIdText);
            }

        }.makeMapFromXMLFile(sequencingsXml);
    };

    protected abstract S makeSourceId(String sourceIdText);
}
