package edu.rice.cs.bioinfo.programs.soranus.views;

import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc3;
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

    public void addDataAddRequestListener(Proc1<File> listener);

    public void addSnitkinTransMapAnalysisRequestedListener(Proc3<DR,DR,DR> listener);

    public void addDataRecordSelectedListener(Proc1<DR> listener);

    public void addNeighborJoiningAnalysisRequestedListener(Proc1<DR> listener);

    public void addAnalysisRecordSelectedListener(Proc1<AR> listener);

    public void startView();

}
