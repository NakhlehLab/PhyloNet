package edu.rice.cs.bioinfo.programs.soranus.models.data;

import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/1/13
 * Time: 3:31 PM
 * To change this template use File | Settings | File Templates.
 */
public interface SequencingsDataProvider<R,P,E extends Throwable>
{
    public Map<Sequencing,P> getSequencingToPatientMap(R recordId) throws E;
}
