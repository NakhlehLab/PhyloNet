package edu.rice.cs.bioinfo.programs.soranus.views;

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.AnalysisRecord;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.DataRecord;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/6/13
 * Time: 4:02 PM
 * To change this template use File | Settings | File Templates.
 */
public interface WorkspaceView<DR extends DataRecord,
                               AR extends AnalysisRecord>
{
    public void addCreateSequencingsDataRequestListener(Proc listener);

    public Observable1<File> getAddDataRequested();

    public Observable1<DR> getNeighborJoiningAnalysisRequested();

    public Observable3<DR,DR,DR> getSnitkinTransMapAnalysisRequested();

    public Observable1<DR> getDetectRecombAnalysisRequested();

    public Observable1<AR> getAnalysisRecordSelected();

    public Observable1<DR> getMinSpanTreesMax2SnpRequested();

    public Observable1<DR> getMinSpanTreesSnpRequested();

    public Observable1<DR> getDataRecordSelected();

    public void startView();

}
